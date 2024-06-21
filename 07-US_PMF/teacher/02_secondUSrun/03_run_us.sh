#!/bin/bash
#SBATCH --job-name=07_USPMF_runUS
#SBATCH --nodes=1
#SBATCH --gres=gpu:0
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --output=slurmsysprep.out
#SBATCH --error=errorsysprep
#SBATCH --export=ALL

ml gromacs/5.1.4

for frame in $(cat important_frames.txt)
do
	i=${frame%%.*}
	echo "running for frame $i..."
	gmx_mpi grompp -f npt_umbrella.mdp -c conf$i.gro -p topol.top -r conf$i.gro -n index.ndx -o npt$i.tpr
	gmx_mpi mdrun -ntomp 4 -v -deffnm npt$i

	gmx_mpi grompp -f md_umbrella.mdp -c npt$i.gro -t npt$i.cpt -p topol.top -r npt$i.gro -n index.ndx -o umbrella$i.tpr
	gmx_mpi mdrun -ntomp 4 -v -deffnm umbrella$i

	echo "umbrella$i.tpr" >> tpr-files.dat
	echo "umbrella${i}_pullf.xvg" >> pullf-files.dat
	echo "done"
done

gmx_mpi wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
