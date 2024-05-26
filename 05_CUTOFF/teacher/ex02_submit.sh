#!/bin/bash

ml gromacs/5.1.4

name=$1
method=$2

out='_output'

# simulations
gmx_mpi grompp -f ${method}_pre.mdp -p ${out}/${name}_dna.top -c $out/${name}_npt.gro -o $out/${name}_${method}_pre.tpr
gmx_mpi mdrun -s $out/${name}_${method}_pre.tpr -deffnm $out/${name}_${method}_pre -ntomp 6
gmx_mpi grompp -f ${method}_prod.mdp -p ${out}/${name}_dna.top -c $out/${name}_${method}_pre.gro -o $out/${name}_${method}_prod.tpr
gmx_mpi mdrun -s $out/${name}_${method}_prod.tpr -deffnm $out/${name}_${method}_prod -ntomp 6

# analysis 
printf "${name}_${method}_prod.log:\t" >> timings_dna.txt 

grep -A1 'Core t' ${out}/${name}_${method}_prod.log | tail -1 | awk '{print $3}' >> timings_dna.txt

printf 'ssDNA\nSystem\n' | gmx_mpi trjconv -f $out/${name}_${method}_prod.xtc \
					   -s $out/${name}_${method}_prod.tpr \
					   -o $out/${name}_${method}_whole.xtc \
					   -n ${out}/${name}_dna.ndx \
					   -pbc mol \
					   -center
printf 'ssDNA\nSystem\n' | gmx_mpi trjconv -f $out/${name}_${method}_pre.gro \
					   -s $out/${name}_${method}_pre.tpr \
					   -o $out/${name}_${method}_whole.gro \
					   -n ${out}/${name}_dna.ndx \
					   -pbc mol \
					   -center

printf 'DNA\nDNA\n' | gmx_mpi rms -f $out/${name}_${method}_whole.xtc \
				  -s $out/${name}_${method}_whole.gro \
				  -o $out/${name}_${method}_rmsd.xvg \
				  -a $out/${name}_${method}_avg_rmsd.xvg

gmx_mpi rdf -f ${out}/${name}_${method}_prod.xtc -s ${out}/${name}_${method}_prod.tpr -o ${out}/rdf_${name}_${method}_O1P_OW.xvg -ref 'name O1P' -sel 'name OW'
