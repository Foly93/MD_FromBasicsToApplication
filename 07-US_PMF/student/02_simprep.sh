#!/bin/bash

gmx grompp -f md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr
gmx mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg

gmx pairdist -f pull.xtc -s pull.tpr -n index.ndx -sel "com of group Chain_A" -ref "com of group Chain_B" -o dist.xvg

sed '/^\@/d' dist.xvg |\
sed '/^\#/d' |\
awk '$2 > i+0.2{a[k++]=$1;i=i+0.2}END{for(j=0;j<501;j++)if(a[j] != 0.000 && j > 0 || j < 1 )printf("%7.3f%1s",a[j],"\n")}' > important_frames.txt

for frame in $(cat important_frames.txt)
do
	printf "0\n" | gmx_mpi trjconv -s pull.tpr -f pull.xtc -o conf${frame%%.*}.gro -b $frame -e $frame
done
