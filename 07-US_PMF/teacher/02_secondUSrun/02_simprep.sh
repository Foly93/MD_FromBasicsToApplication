#!/bin/bash
#SBATCH --job-name=07_USPMF_simprep \
#SBATCH --nodes=1 \
#SBATCH --gres=gpu:0 \
#SBATCH --time=1-00:00:00 \
#SBATCH --ntasks=16 \
#SBATCH --cpus-per-task=1 \
#SBATCH --output=slurmsysprep.out \
#SBATCH --error=errorsysprep \
#SBATCH --export=ALL \

ml gromacs/5.1.4

# Run the pulling simulation
gmx_mpi grompp -f md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr
gmx_mpi mdrun -ntomp 4 -deffnm pull -pf pullf.xvg -px pullx.xvg

gmx_mpi pairdist -f pull.xtc -s pull.tpr -n index.ndx -sel "com of group Chain_A" -ref "com of group Chain_B" -o dist.xvg

# delete all lines starting with @
sed '/^\@/d' dist.xvg |\
# delete all lines starting with #
sed '/^\#/d' |\
# if another spacing is wanted the i+0.2 can be changed into i+spacing
awk '$2 > i+0.2{a[k++]=$1;i=i+0.2}END{for(j=0;j<501;j++)if(a[j] != 0.000 && j > 0 || j < 1 )printf("%7.3f%1s",a[j],"\n")}' > important_frames.txt	# write out every frame that has COM dist 0.2 higher than the prior one

# generate the coordinate files from the important frames
for frame in $(cat important_frames.txt)
do
	printf "0\n" | gmx_mpi trjconv -s pull.tpr -f pull.xtc -o conf${frame%%.*}.gro -b $frame -e $frame
done
