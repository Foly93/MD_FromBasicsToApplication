#!/bin/bash
#SBATCH --job-name=08_PROT_PCA
#SBATCH --nodes=1
#SBATCH --gres=gpu:0
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm08.out
#SBATCH --error=error08
#SBATCH --export=ALL

ml gromacs/5.1.4

# choose AMBER99-ILDN and TIP4P-Ew and build the system; beforehand: ligands and other solutes were removed manually
printf "6\n3\n" | gmx_mpi pdb2gmx -f 5yok_clean.pdb -ignh -ter -o 5yok_clean.gro
gmx_mpi editconf -f 5yok_clean.gro -o 5yok_cube.gro -c -d 0.85 -bt cubic
gmx_mpi solvate -cp 5yok_cube.gro -cs tip4p.gro -o 5yok_solv.gro -p topol.top
gmx_mpi grompp -f em.mdp -c 5yok_solv.gro -p topol.top -o ions.tpr
printf "SOL\n" | gmx_mpi genion -s ions.tpr -o 5yok_ready.gro -p topol.top -pname NA -nname CL -conc 0.1 -neutral

# Run simulation loop; first gro file is special, since it is the output from prior simulation/setup
j='5yok_ready'
for i in em nvt npt pre prod
do
	gmx_mpi grompp -f $i.mdp -c $j.gro -p topol.top -o $i.tpr
	gmx_mpi mdrun -ntomp 16 -deffnm $i -v 
	j=$i
done

# convert the trj for better analysis with gyrate and rms
printf "Protein\nSystem\n" | gmx_mpi trjconv -s prod.tpr -f prod.xtc -o prod_center.xtc -pbc mol -center
printf "Protein\n" | gmx_mpi gyrate -f prod_center.xtc -s prod.tpr -o rgyr.xvg
printf "Protein\nProtein\n" | gmx_mpi rms -s em.tpr -f prod_center.xtc -o rmsd.xvg

# BSE analysis for both rmsd and rgyr
gmx_mpi analyze -f rmsd.xvg -ee bse_rmsd.xvg
gmx_mpi analyze -f rgyr.xvg -ee bse_rgyr.xvg

# Eigenvector analysis of the covariance matrix 
printf "Protein\nProtein\n" | gmx_mpi covar -f prod_center.xtc -s prod.tpr -xpma 
gmx_mpi xpm2ps -f covara.xpm -o covara.eps
printf "Protein\nProtein\n" | gmx_mpi anaeig -v eigenvec.trr -f prod_center.xtc -s prod.tpr -proj eigproj1_10.xvg -comp eigcomp1_10.xvg -last 10 
printf "Protein\nProtein\n" | gmx_mpi anaeig -v eigenvec.trr -f prod_center.xtc -s prod.tpr -extr  # vmd is recommended for visualisation (pymol only shows one chain smh)

# Plot the energy properties from the equilibrations
printf "Potential\n" | gmx_mpi energy -f em.edr -o epot.xvg
printf "Temperature\n" | gmx_mpi energy -f nvt.edr -o temp.xvg
printf "Pressure\n" | gmx_mpi energy -f npt.edr -o pres.xvg
printf "Density\n" | gmx_mpi energy -f npt.edr -o dens.xvg

# Cleanup
mkdir -p _output _trash
mv -v *.cpt *.edr *.gro *.log *.tpr *.trr *.xtc *.xvg _output
mv -v \#* mdout.mdp _trash

# get timings for the md simulations 
echo "total/s per.CPU/s Acc" | awk '{printf "%12s %11s %6s\n", $1,$2,$3}' > timings.txt
for log in $(ls _LOG/*)
do
	grep 'Time:' $log | awk '{printf "%12.3f %11.3f %6.1f\n", $2,$3,$4}' >> timings.txt
done

sum_total=$(awk 'NR > 1{print $1}' timings.txt | paste -sd+ | bc)
sum_pCPU=$(awk 'NR > 1{print $2}' timings.txt | paste -sd+ | bc)

echo "$sum_total $sum_pCPU" | awk '{printf "%12.3f %11.3f",$1,$2}' >> timings.txt

