#!/bin/bash 
gmx grompp -f polymer.mdp -c polymer.gro -p polymer_fj.top -o polymer_fj.tpr
gmx mdrun -deffnm polymer_fj -v 
gmx grompp -f polymer.mdp -c polymer.gro -p polymer_fr.top -o polymer_fr.tpr
gmx mdrun -deffnm polymer_fr -v 
gmx grompp -f polymer.mdp -c polymer.gro -p polymer_bonded.top -o polymer_bonded.tpr
gmx mdrun -deffnm polymer_bonded -v 
gmx grompp -f polymer.mdp -c polymer.gro -p polymer_non_bonded.top -o polymer_non_bonded.tpr
gmx mdrun -deffnm polymer_non_bonded -v 

gmx mk_angndx -s polymer_non_bonded.tpr -n dihedral.ndx -type dihedral
gmx mk_angndx -s polymer_non_bonded.tpr -n angle.ndx -type angle

echo "This file summarizes the averages of the polymer analysis" > summary.txt

for i in fj fr bonded non_bonded
do
	printf "\n\nPOLYMER MODEL: $i\n" >> summary.txt
	gmx angle -f polymer_$i.xtc -n angle.ndx -od ${i}_angle_hist.xvg -type angle -all -ov ${i}_angle.xvg >> summary.txt
	gmx angle -f polymer_$i.xtc -n dihedral.ndx -od ${i}_dihedral_hist.xvg -type dihedral -all -ov ${i}_dihedral.xvg >> summary.txt
	echo "0" | gmx polystat -s polymer_$i.tpr -f polymer_$i.xtc -o ${i}_polystat.xvg >> summary.txt
done
