#!/bin/bash

ml gromacs/5.1.4

name=$1
method=$2

out='_output'

# Simulation Runs
gmx_mpi grompp -f ${method}_pre.mdp -p $out/${name}.top -c $out/${name}_npt.gro -o $out/${name}_${method}_pre.tpr
gmx_mpi mdrun -s $out/${name}_${method}_pre.tpr -deffnm $out/${name}_${method}_pre -ntomp 6
gmx_mpi grompp -f ${method}_prod.mdp -p $out/${name}.top -c $out/${name}_${method}_pre.gro -o $out/${name}_${method}_prod.tpr
gmx_mpi mdrun -s $out/${name}_${method}_prod.tpr -deffnm $out/${name}_${method}_prod -ntomp 6

# some analysis
printf "${out}/${name}_${method}_prod.log:\t" >> timings_water.txt
grep -A1 'Core t' ${out}/${name}_${method}_prod.log | tail -1 | awk '{print $3}' >> timings_water.txt
gmx_mpi rdf -f ${out}/${name}_${method}_prod.xtc -s ${out}/${name}_${method}_prod.tpr -o ${out}/rdf_${name}_${method}_OW_OW.xvg -ref 'name OW' -sel 'name OW'
