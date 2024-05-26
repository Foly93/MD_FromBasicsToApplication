#!/bin/bash

# Formula for calculating Specific Heat Capacity C_v
#
#        ( RMSD(E_tot) * 1000 )^2 * N_mol   
# C_v = ----------------------------------
#                  R * T * T                                                         

# Version and Hardware information ### uncomment to erase results of prior runs
gmx --version | tail -24 | head -16 #> summary.txt

# polymer run with the Berendsen Thermostat to be compared with run from 02_*/
gmx grompp -f polymer.mdp -c polymer.gro -p polymer_non_bonded.top -o polymer_non_bonded.tpr
gmx mdrun -deffnm polymer_non_bonded -v 

# 3 Point water Models (SPC SPCE TIP3P)
for water in spc spce tip3p tip4p
do
	printf "\n\n\n### $water" >> summary.txt
	# EM
	gmx grompp -f em.mdp -c ${water}box.gro -p $water.top -o em_$water.tpr
	gmx mdrun -v -deffnm em_$water
	
	for tstat in nve berendsen nosehoover nhc
	do
		printf "\n##$tstat\n" >> summary.txt
		# Thermostat Simulation Equilibration
		gmx grompp -f eq_${tstat}.mdp -c em_$water.gro -p $water.top -o eq_${tstat}_$water.tpr
		gmx mdrun -v -deffnm eq_${tstat}_$water
		# Thermostat Production Run
		gmx grompp -f prod_${tstat}.mdp -c eq_${tstat}_$water.gro -p $water.top -o prod_${tstat}_$water.tpr
		gmx mdrun -v -deffnm prod_${tstat}_$water
		printf "Temperature\n7\n6\n5\n" | gmx energy -f prod_${tstat}_$water.edr -s prod_${tstat}_$water.tpr -o ${tstat}_$water.xvg -fluct_props -nmol 523 |  grep -E '\(kJ/mol|K\)' >> summary.txt
	done
done
