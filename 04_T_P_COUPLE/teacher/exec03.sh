#!/bin/bash

# SOLUTION: The ensemble that is simulated is no longer the NPT ensemble but the NPH ensemble. 
# SOLUTION: Therefore the energy values drift away from the ensemble averages of the NPT ensemble.
# SOLUTION: A way to test this hypothesis is to calculate properties related to the NPH ensemble like the Joule-Thompson Coefficient.

#printf "\n\nExec 03 - Breaking the System\n" >> summary.txt

# NPH
gmx grompp -f standard.mdp -c spcebox_nvt.gro -p spcebox.top -o standard.tpr -maxwarn 2
gmx mdrun -v -s standard.tpr -deffnm standard

gmx grompp -f non_standard.mdp -c spcebox_nvt.gro -p spcebox.top -o non_standard.tpr -maxwarn 2
gmx mdrun -v -s non_standard.tpr -deffnm non_standard

printf "\nstandard\n" >> summary.txt
printf "Density\n" | gmx energy -f standard.edr -o dens_standard.xvg | egrep '(kg/m\^3)' >> summary.txt
printf "Enthalpy\n" | gmx energy -f standard.edr -o enth_standard.xvg | egrep '(kJ/mol)' >> summary.txt
printf "Potential\n" | gmx energy -f standard.edr -o epot_standard.xvg | egrep '(kJ/mol)' >> summary.txt
printf "Kinetic-En.\n" | gmx energy -f standard.edr -o ekin_standard.xvg | egrep '(kJ/mol)' >> summary.txt
printf "Total-Energy\n" | gmx energy -f standard.edr -o etot_standard.xvg | egrep '(kJ/mol)' >> summary.txt
printf "\nnon_standard\n" >> summary.txt
printf "Density\n" | gmx energy -f non_standard.edr -o dens_non_standard.xvg | egrep '(kg/m\^3)' >> summary.txt
printf "Enthalpy\n" | gmx energy -f non_standard.edr -o enth_non_standard.xvg | egrep '(kJ/mol)' >> summary.txt
printf "Potential\n" | gmx energy -f non_standard.edr -o epot_non_standard.xvg | egrep '(kJ/mol)' >> summary.txt
printf "Kinetic-En.\n" | gmx energy -f non_standard.edr -o ekin_non_standard.xvg | egrep '(kJ/mol)' >> summary.txt
printf "Total-Energy\n" | gmx energy -f non_standard.edr -o etot_non_standard.xvg | egrep '(kJ/mol)' >> summary.txt

