#!/bin/bash

weight=12
gmx grompp -f em.mdp -p topol_${weight}.top -c intro.gro -o em_${weight}.tpr
gmx mdrun -v -deffnm em_${weight}

for timestep in 0.007 0.009 0.011
do
	sed "s/TIMESTEP/$timestep/g" intro.mdp > intro_${timestep}.mdp
	gmx grompp -f intro_${timestep}.mdp -p topol_${weight}.top -c em_${weight}.gro -o intro_${timestep}_${weight}.tpr -maxwarn 3
	gmx mdrun -v -deffnm intro_${timestep}_${weight}
	printf "3\n4\n5\n" | gmx energy -f intro_${timestep}_${weight}.edr -s intro_${timestep}_${weight}.tpr -o intro_${timestep}_${weight}.xvg
done

weight=6
gmx grompp -f em.mdp -p topol_${weight}.top -c intro.gro -o em_${weight}.tpr
gmx mdrun -v -deffnm em_${weight}

for timestep in 0.005 0.006 0.007 
do
	sed "s/TIMESTEP/$timestep/g" intro.mdp > intro_${timestep}.mdp
	gmx grompp -f intro_${timestep}.mdp -p topol_${weight}.top -c em_${weight}.gro -o intro_${timestep}_${weight}.tpr -maxwarn 3
	gmx mdrun -v -deffnm intro_${timestep}_${weight}
	printf "3\n4\n5\n" | gmx energy -f intro_${timestep}_${weight}.edr -s intro_${timestep}_${weight}.tpr -o intro_${timestep}_${weight}.xvg
done

weight=4
gmx grompp -f em.mdp -p topol_${weight}.top -c intro.gro -o em_${weight}.tpr
gmx mdrun -v -deffnm em_${weight}

for timestep in 0.004 0.005 0.006
do
	sed "s/TIMESTEP/$timestep/g" intro.mdp > intro_${timestep}.mdp
	gmx grompp -f intro_${timestep}.mdp -p topol_${weight}.top -c em_${weight}.gro -o intro_${timestep}_${weight}.tpr -maxwarn 3
	gmx mdrun -v -deffnm intro_${timestep}_${weight}
	printf "3\n4\n5\n" | gmx energy -f intro_${timestep}_${weight}.edr -s intro_${timestep}_${weight}.tpr -o intro_${timestep}_${weight}.xvg
done

weight=2
gmx grompp -f em.mdp -p topol_${weight}.top -c intro.gro -o em_${weight}.tpr
gmx mdrun -v -deffnm em_${weight}

for timestep in 0.003 0.0034 0.0038
do
	sed "s/TIMESTEP/$timestep/g" intro.mdp > intro_${timestep}.mdp
	gmx grompp -f intro_${timestep}.mdp -p topol_${weight}.top -c em_${weight}.gro -o intro_${timestep}_${weight}.tpr -maxwarn 3
	gmx mdrun -v -deffnm intro_${timestep}_${weight}
	printf "3\n4\n5\n" | gmx energy -f intro_${timestep}_${weight}.edr -s intro_${timestep}_${weight}.tpr -o intro_${timestep}_${weight}.xvg
done

