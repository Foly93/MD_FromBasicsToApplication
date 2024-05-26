#!/bin/bash

ml gromacs/5.1.4


out='_output'
tr='_trash'
mkdir -p $out
mkdir -p $tr

system_creation()
{
	printf '6\n6\n' | gmx_mpi pdb2gmx -f 10mer_${1}-DNA.pdb -o $out/${1}_dna.gro -ignh -p ${out}/${1}_dna.top
	gmx_mpi solvate -cp $out/${1}_dna.gro -cs spc216.gro -box 5 5 5 -p ${out}/${1}_dna.top -o $out/${1}_box.gro
	gmx_mpi grompp -f em.mdp -p ${out}/${1}_dna.top -c $out/${1}_box.gro -o $out/${1}_ions.tpr
	printf 'SOL\n' | gmx_mpi genion -s $out/${1}_ions.tpr -p ${out}/${1}_dna.top -o $out/${1}_ions.gro -neutral
}

system_equilibration()
{
	gmx_mpi grompp -f em.mdp -p ${out}/${1}_dna.top -c $out/${1}_ions.gro -o $out/${1}_em.tpr
	gmx_mpi mdrun -s $out/${1}_em.tpr -deffnm $out/${1}_em -ntomp 2
	
	gmx_mpi grompp -f nvt.mdp -p ${out}/${1}_dna.top -c $out/${1}_em.gro -o $out/${1}_nvt.tpr
	gmx_mpi mdrun -s $out/${1}_nvt.tpr -deffnm $out/${1}_nvt -ntomp 2
	
	gmx_mpi grompp -f npt.mdp -p ${out}/${1}_dna.top -c $out/${1}_nvt.gro -o $out/${1}_npt.tpr
	gmx_mpi mdrun -s $out/${1}_npt.tpr -deffnm $out/${1}_npt -ntomp 2
}

system_production()
{
	# Prepare files to generate solutions 
	printf "Filename\tCore Time (s)\n" > timings_dna.txt
	printf 'a 1-315\nname 9 ssDNA\nq\n' | gmx_mpi make_ndx -f $out/${1}_npt.gro -o ${out}/${1}_dna.ndx

	# send the simulations to the HEN
	for method in tr8 tr12 tr18 sh8 sh12 pme
	do
		sbatch \
		--job-name=${method}_${1} \
		--nodes=1 \
		--gres=gpu:1 \
		--time=1-00:00:00 \
		--ntasks=4 \
		--cpus-per-task=1 \
		--output=${out}/${1}_${method}_slurm.out \
		--export=none \
		./ex02_submit.sh $1 $method
	done
}

for system in b a
do
	system_creation $system
	system_equilibration $system
	system_production $system
done
