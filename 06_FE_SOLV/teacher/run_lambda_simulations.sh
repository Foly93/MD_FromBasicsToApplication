#!/bin/bash

ml gromacs/5.1.4

l=$1

OUTDIR="_output/lambda_$l"
gro_file="methane_water.gro"

mkdir -p $OUTDIR

# Simulation Runs
# loop over all different simulation steps
for sim_phase in em nvt npt prod
do
	# Concatenate the according banner to give information about the current run (just for fun)
	cat $sim_phase.banner

	# start the run of the corresponding lambda and simulation step
	# for '--gres:gpu:0' in 'submit_script.sh', set '-ntomp 6' after 'mdrun'
	gmx_mpi grompp \
		-f ${sim_phase}_$l.mdp \
		-c $gro_file \
		-p topol.top \
		-o $OUTDIR/$sim_phase.tpr \
		-po $OUTDIR/mdout$sim_phase.mdp \
		&>$OUTDIR/grompp.$sim_phase \
		|| { echo "something went wrong, check $OUTDIR/grompp.$sim_phase"; exit; }
	gmx_mpi mdrun \
		-deffnm $OUTDIR/$sim_phase &>$OUTDIR/mdrun.$sim_phase \
		|| { echo "something went wrong, check $OUTDIR/mdrun.$sim_phase"; exit; }

	gro_file="$OUTDIR/$sim_phase.gro"
done
