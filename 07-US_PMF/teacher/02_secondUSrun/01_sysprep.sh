#!/bin/bash
#SBATCH --job-name=07_USPMF_sysprep \
#SBATCH --nodes=1 \
#SBATCH --gres=gpu:0 \
#SBATCH --time=1-00:00:00 \
#SBATCH --ntasks=16 \
#SBATCH --cpus-per-task=1 \
#SBATCH --output=slurmsysprep.out \
#SBATCH --error=errorsysprep \
#SBATCH --export=ALL \

ml gromacs/5.1.4

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


