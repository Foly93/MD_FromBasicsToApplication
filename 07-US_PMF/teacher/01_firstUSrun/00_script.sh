#!/bin/bash
#SBATCH --job-name=07_USPMF
#SBATCH --nodes=1
#SBATCH --gres=gpu:0
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --output=slurmsysprep.out
#SBATCH --error=errorsysprep
#SBATCH --export=ALL

ml gromacs/5.1.4

### 01_SYSPREP ###
# create topologies and coordinates for GROMACS from PDB file
printf "13\n1\n2\n0\n2\n0\n2\n0\n2\n0\n2\n0\n" | gmx_mpi pdb2gmx -f 2BEG_model1_capped.pdb -ignh -ter -o complex.gro

# this line is necessary for the pulling simulation
printf "#ifdef POSRES_B\n#include \"posre_Protein_chain_B.itp\"\n#endif\n" >> topol_Protein_chain_B.itp

# place the protofibril at the correct location and create the box, values from online tutorial.
gmx_mpi editconf -f complex.gro -o newbox.gro -center 3.280 2.181 2.4775 -box 6.560 4.362 12

# solvate the newly created box
gmx_mpi solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top

# ion generation in GROMACS requires a TPR file, therefore use grompp then genion...
gmx_mpi grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr
printf "13\n" | gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1

# run minimisation and equilibration
gmx_mpi grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx_mpi mdrun -ntomp 4 -v -deffnm em
gmx_mpi grompp -f npt.mdp -c em.gro -p topol.top -r em.gro -o npt.tpr
gmx_mpi mdrun -ntomp 4 -deffnm npt

# make an index that contains chain a and chain b explicitly
gmx_mpi make_ndx -f npt.gro<<EOF
r 1-27
name 19 Chain_A
r 28-54
name 20 Chain_B
q
EOF

### 02_SIMPREP ###
# Run the pulling simulation
gmx_mpi grompp -f md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr
gmx_mpi mdrun -ntomp 4 -deffnm pull -pf pullf.xvg -px pullx.xvg

gmx_mpi pairdist -f pull.xtc -s pull.tpr -n index.ndx -sel "com of group Chain_A" -ref "com of group Chain_B" -o dist.xvg

# delete all lines starting with @
sed '/^\@/d' dist.xvg |\
# delete all lines starting with #
sed '/^\#/d' |\
# if another spacing is wanted the i+0.2 can be changed into i+spacing
awk '$2 > i+0.2{a[k++]=$1;i=i+0.2}END{for(j=0;j<501;j++)if(a[j] != 0.000 && j > 0 || j < 1 )printf("%7.3f%1s",a[j],"\n")}' > important_frames.txt       # write out every frame that has COM dist 0.2 higher than the prior one

### 03_RUN_US ###
# generate the coordinate files from the important frames
for frame in $(cat important_frames.txt)
do
        printf "0\n" | gmx_mpi trjconv -s pull.tpr -f pull.xtc -o conf${frame%%.*}.gro -b $frame -e $frame
done

for frame in $(cat important_frames.txt)
do
        i=${frame%%.*}
        echo "running for frame $i..."
        gmx_mpi grompp -f npt_umbrella.mdp -c conf$i.gro -p topol.top -r conf$i.gro -n index.ndx -o npt$i.tpr
        gmx_mpi mdrun -ntomp 4 -v -deffnm npt$i

        gmx_mpi grompp -f md_umbrella.mdp -c npt$i.gro -t npt$i.cpt -p topol.top -r npt$i.gro -n index.ndx -o umbrella$i.tpr
        gmx_mpi mdrun -ntomp 4 -v -deffnm umbrella$i

        echo "umbrella$i.tpr" >> tpr-files.dat
        echo "umbrella${i}_pullf.xvg" >> pullf-files.dat
        echo "done"
done

gmx_mpi wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal



### 04_US_ADDONS ###
extraframe=( 50 147 186 195 257 302 340 393 ) # already added: ...

if [[ -z "$extraframe" ]]; then
   echo "no frames specified to run an umbrella simulation for." ; exit 1
else
        for i in "${extraframe[@]}"
        do
                echo "running for frame $i"
                printf "0\n" | gmx_mpi trjconv -s pull.tpr -f pull.xtc -o conf$i.gro -b $i -e $i

                gmx_mpi grompp -f npt_umbrella.mdp -c conf$i.gro -p topol.top -r conf$i.gro -n index.ndx -o npt$i.tpr
                gmx_mpi mdrun -ntomp 4 -v -deffnm npt$i

                gmx_mpi grompp -f md_umbrella.mdp -c npt$i.gro -t npt$i.cpt -p topol.top -r npt$i.gro -n index.ndx -o umbrella$i.tpr
                gmx_mpi mdrun -ntomp 4 -v -deffnm umbrella$i

                echo "umbrella$i.tpr" >> tpr-files.dat
                echo "umbrella${i}_pullf.xvg" >> pullf-files.dat
                echo "done"
        done
fi
gmx_mpi wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal


### 05_CLEAN_DIR ### 
mkdir -p _CPT  _EDR  _GRO  _LOG  _TPR  _trash  _TRR  _XTC  _XVG

mv -v *.cpt _CPT
mv -v *.edr _EDR
mv -v *.gro _GRO
mv -v *.log _LOG
mv -v *.tpr _TPR
mv -v *.trr _TRR
mv -v *.xtc _XTC
mv -v *.xvg _XVG
mv -v \#* _trash

### 06_TIMINGS ###
echo "total/s per.CPU/s Acc" | awk '{printf "%12s %11s %6s\n", $1,$2,$3}' > timings.txt
for log in $(ls _LOG/*)
do
	grep 'Time:' $log | awk '{printf "%12.3f %11.3f %6.1f\n", $2,$3,$4}' >> timings.txt
done

sum_total=$(awk 'NR > 1{print $1}' timings.txt | paste -sd+ | bc)
sum_pCPU=$(awk 'NR > 1{print $2}' timings.txt | paste -sd+ | bc)

echo "$sum_total $sum_pCPU" | awk '{printf "%12.3f %11.3f",$1,$2}' >> timings.txt

