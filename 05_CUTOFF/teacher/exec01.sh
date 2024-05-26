#!/bin/bash
ml gromacs/5.1.4

out="_output"
mkdir -p $out
mkdir -p _trash

# create system: solvate -> insert -> pdb2gmx
gmx_mpi solvate -cs spc216 -box 4 4 4 -o $out/spce.gro
printf "6\n" | gmx_mpi pdb2gmx -f $out/spce.gro -p $out/spce.top -water spce -o _trash/pdb2gmx.gro

# EM
gmx_mpi grompp -f em.mdp -c $out/spce.gro -p $out/spce.top -o $out/spce_em.tpr
gmx_mpi mdrun -s $out/spce_em.tpr -deffnm $out/spce_em -v

# NVT
gmx_mpi grompp -f nvt.mdp -c $out/spce_em.gro -p $out/spce.top -o $out/spce_nvt.tpr
gmx_mpi mdrun -s $out/spce_nvt.tpr -deffnm $out/spce_nvt -v 

# NPT 
gmx_mpi grompp -f npt.mdp -c $out/spce_nvt.gro -p $out/spce.top -o $out/spce_npt.tpr
gmx_mpi mdrun -s $out/spce_npt.tpr -deffnm $out/spce_npt -v

printf "filename\tCore time (s)\n" > timings_water.txt
for method in tr8 tr12 tr18 sh8 sh12 pme
do
	sbatch \
	--job-name=${method}_spce \
	--nodes=1 \
	--gres=gpu:1 \
	--time=1-00:00:00 \
	--ntasks=4 \
	--cpus-per-task=1 \
	--output=${out}/spce_${method}_slurm.out \
	--error=${out}/spce_${method}_error \
	--export=none \
	./ex01_submit.sh "spce" $method
done

