#!/bin/bash

printf "\n\nExec 02 - Watermodel Comparison\n" >> summary.txt
# create system: solvate -> insert -> pdb2gmx
gmx solvate -cs tip4p -box 2.504 2.504 2.504 -maxsol 1 -o tip4p.gro
gmx solvate -cs spc216 -box 2.504 2.504 2.504 -maxsol 1 -o spce.gro
gmx insert-molecules -ci tip4p.gro -box 2.504 2.504 2.504 -nmol 523 -try 10000 -o tip4pbox.gro
gmx insert-molecules -ci spce.gro -box 2.504 2.504 2.504 -nmol 523 -try 10000 -o spcebox.gro
printf "6\n" | gmx pdb2gmx -f tip4pbox.gro -p tip4pbox.top -water tip4p -o ../_trash/pdb2gmx.gro 
printf "6\n" | gmx pdb2gmx -f spcebox.gro -p spcebox.top -water spce -o ../_trash/pdb2gmx.gro 

# EM
gmx grompp -f em.mdp -c spcebox.gro -p spcebox.top -o spcebox_em.tpr
gmx grompp -f em.mdp -c tip4pbox.gro -p tip4pbox.top -o tip4pbox_em.tpr
gmx mdrun -s spcebox_em.tpr -deffnm spcebox_em
gmx mdrun -s tip4pbox_em.tpr -deffnm tip4pbox_em

# NVT
gmx grompp -f nvt.mdp -c spcebox_em.gro -p spcebox.top -o spcebox_nvt.tpr
gmx grompp -f nvt.mdp -c tip4pbox_em.gro -p tip4pbox.top -o tip4pbox_nvt.tpr
gmx mdrun -s spcebox_nvt.tpr -deffnm spcebox_nvt
gmx mdrun -s tip4pbox_nvt.tpr -deffnm tip4pbox_nvt

# NPT
gmx grompp -f waters_npt.mdp -c spcebox_nvt.gro -p spcebox.top -o spcebox_npt.tpr
gmx grompp -f waters_npt.mdp -c tip4pbox_nvt.gro -p tip4pbox.top -o tip4pbox_npt.tpr
gmx mdrun -s spcebox_npt.tpr -deffnm spcebox_npt
gmx mdrun -s tip4pbox_npt.tpr -deffnm tip4pbox_npt

# PRE
gmx grompp -f waters_pre.mdp -c spcebox_npt.gro -p spcebox.top -o spcebox_pre.tpr
gmx grompp -f waters_pre.mdp -c tip4pbox_npt.gro -p tip4pbox.top -o tip4pbox_pre.tpr
gmx mdrun -s spcebox_pre.tpr -deffnm spcebox_pre
gmx mdrun -s tip4pbox_pre.tpr -deffnm tip4pbox_pre

# PROD
gmx grompp -f waters_prod.mdp -c spcebox_pre.gro -p spcebox.top -o spcebox_prod.tpr
gmx grompp -f waters_prod.mdp -c tip4pbox_pre.gro -p tip4pbox.top -o tip4pbox_prod.tpr
gmx mdrun -v -s spcebox_prod.tpr -deffnm spcebox_prod
gmx mdrun -v -s tip4pbox_prod.tpr -deffnm tip4pbox_prod

# DATA ANALYSIS
printf "TIP4P " >> summary.txt
printf "Volume\nTemperature\n" | gmx energy -f tip4pbox_prod.edr -o tip4pbox_prod.xvg -fluct_props -nmol 523 | egrep -o 'Kappa.*' >> summary.txt
printf "SPC/E " >> summary.txt
printf "Volume\nTemperature\n" | gmx energy -f spcebox_prod.edr -o spcebox_prod.xvg -fluct_props -nmol 523 | egrep -o 'Kappa.*' >> summary.txt
