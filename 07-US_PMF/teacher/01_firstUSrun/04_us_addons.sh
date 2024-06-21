#!/bin/bash
#SBATCH --job-name=07_USPMF_addon
#SBATCH --nodes=1
#SBATCH --gres=gpu:0
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --output=slurmsysprep.out
#SBATCH --error=errorsysprep
#SBATCH --export=ALL

ml gromacs/5.1.4

extraframe=( 149 184 202 210 219 229 244 267 297 318 333 372 408 427 455 482 ) # already added: 50

if [[ -z "$extraframe" ]]; then
   echo "no frames specified to run an umbrella simulation for." ; exit 1
else
	for i in "${extraframe[@]}"
	do
		echo "running for frame $i"
		printf "0\n" | gmx_mpi trjconv -s pull.tpr -f pull.xtc -o conf$i.gro -b $i -e $i

        	gmx_mpi grompp -f npt_umbrella.mdp -c conf$i.gro -p topol.top -r conf$i.gro -n index.ndx -o npt$i.tpr
        	gmx_mpi mdrun -ntomp 4 -v -deffnm npt$i

        	gmx_mpi grompp -f md_umbrella.mdp -c npt$i.gro -t npt$i.cpt -p topol.top -r npt$i.gro -n index.ndx -o umbrella$i.tpr
        	gmx_mpi mdrun -ntomp 4 -v -deffnm umbrella$i

		echo "umbrella$i.tpr" >> tpr-files.dat
		echo "umbrella${i}_pullf.xvg" >> pullf-files.dat
        	echo "done"
	done
fi

gmx_mpi wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
