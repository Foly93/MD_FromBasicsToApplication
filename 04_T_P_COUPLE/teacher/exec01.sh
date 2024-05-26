#!/bin/bash
# Gromacs Formula for calculating Isothermal Compressibility k_T
#
#        RMSD(Vol) * RMSD(Vol)                                                          
# k_T = -----------------------                                                                     
#           k_B * T * <Vol>                                                               

# create system: solvate -> insert -> pdb2gmx
gmx solvate -cs spc216 -box 2.504 2.504 2.504 -maxsol 1 -o tip3p.gro
gmx insert-molecules -ci tip3p.gro -box 2.504 2.504 2.504 -nmol 523 -try 10000 -o choice_tauP.gro
printf "6\n" | gmx pdb2gmx -f choice_tauP.gro -p choice_tauP.top -water tip3p -o ../_trash/pdb2gmx.gro 

# EM
gmx grompp -f em.mdp -c choice_tauP.gro -p choice_tauP.top -o choice_tauP_em.tpr
gmx mdrun -s choice_tauP_em.tpr -deffnm choice_tauP_em

# NVT
gmx grompp -f nvt.mdp -c choice_tauP_em.gro -p choice_tauP.top -o choice_tauP_nvt.tpr
gmx mdrun -s choice_tauP_nvt.tpr -deffnm choice_tauP_nvt

# make summary file
gmx --version | tail -24 | head -16 > summary.txt
printf "\n\nExec 01 - Barostat Comparison"

# Exercise 1 main loop
for bstat in Berendsen Parrinello-Rahman
do
	printf "\n$bstat\n" >> summary.txt
	sed "s/BAROSTAT/$bstat/g" choice_npt.mdp > choice_npt_${bstat:0:3}.mdp

	for taup in 0.05 0.1 0.5 1 5 10
	do
		printf "tauP:    $taup\n" >> summary.txt
		printf "RMSD(Vol):   " >> summary.txt

		sed "s/TAUP/$taup/g" choice_npt_${bstat:0:3}.mdp > choice_npt_${bstat:0:3}_${taup}.mdp
		if (( $(echo "$taup > 0.1" | bc -l) )); then sed -i '/^nstpcouple.*/d' choice_npt_Par_${taup}.mdp; fi
		sed -i '/^nstpcouple.*/d' choice_npt_Ber_*.mdp

		gmx grompp -f choice_npt_${bstat:0:3}_${taup}.mdp -c choice_tauP_nvt.gro -p choice_tauP.top -o choice_npt_${bstat:0:3}_${taup}.tpr -maxwarn 2
		gmx mdrun -v -s choice_npt_${bstat:0:3}_${taup}.tpr -deffnm choice_npt_${bstat:0:3}_${taup}


		printf "Volume\n" |  gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o vol_choice_npt_${bstat:0:3}_${taup}.xvg | grep 'nm^3' | awk '{print $4}' >> summary.txt
		printf "Pressure\n" |  gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o pres_choice_npt_${bstat:0:3}_${taup}.xvg
		printf "Temperature\n" |  gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o temp_choice_npt_${bstat:0:3}_${taup}.xvg
		printf "Total-Energy\n" |  gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o etot_choice_npt_${bstat:0:3}_${taup}.xvg
		printf "Volume\nTemperature" | gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o ../_trash/temp.xvg -fluct_props -nmol 523 | egrep -o 'Kappa.*' >> summary.txt
	done
done
