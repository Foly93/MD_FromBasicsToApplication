#!/bin/bash

# loop over all values of lambda
# $() is called a command substitution, where the output of that command
# can be transferred into a new command or, in this case into 
# a variable used in a for loop for more information google 'command substitution bash'
# and for more information about seq, google it or type 'man seq' in your terminal
for l in $(seq <<<USE THE CORRECT ARGUMENTS FOR 'seq'>>>)
do
	echo "Starting minimization for lambda = $l ..." 

	# Define Variables
	# the following variable is very important for the md runs. To set it correctly, of course one
	# would have to understand what it is used for. in vim you can search a pattern by typing e.g. '/pattern' + <enter>
	# you can find the next/previous instance by hitting <n>/<N> keys.
	# What are they used for? Can you set them freely or are they constrained to certain values?
	OUTDIR="lambda_$l"
	gro_file=<<<INSERT THE CORRECT VALUE FOR THIS VARIABLE>>>
	
	mkdir -p $OUTDIR

	# loop over all different simulation steps which are identified by the MDP files
	for sim_phase in <<<INSERT THE CORRECT ITERABLES>>>
	do
		# Concatenate the according banner to give information about the current run (just for fun)
		cat $sim_phase.banner 
		
		# start the run of the corresponding lambda and simulation step
		gmx grompp \
			-f ${sim_phase}_$l.mdp \
			-c $gro_file \
			-p topol.top \
			-o $OUTDIR/$sim_phase.tpr \
			-po $OUTDIR/mdout$sim_phase.mdp \
			&>$OUTDIR/grompp.$sim_phase \
			|| { echo "something went wrong, check $OUTDIR/grompp.$sim_phase"; exit; }
		gmx mdrun \
			-deffnm <<<INSERT THE CORRECT VALUE FOR THIS FLAG>>> &>$OUTDIR/mdrun.$sim_phase \
			|| { echo "something went wrong, check $OUTDIR/mdrun.$sim_phase"; exit; }

		gro_file="$OUTDIR/$sim_phase.gro"
	<<<YOU ARE NOT done YET. SOMETHING IS MISSING HERE...>>>
done
