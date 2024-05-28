#!/bin/bash

# loop over all values of lambda
# $() is called a command substitution, where the output of that command
# can be transferred into a new command or, in this case into 
# a variable used in a for loop for more information google 'command substitution bash'
# and for more information about seq, google it or type 'man seq' in your terminal
for l in $(seq <<<USE THE CORRECT ARGUMENTS FOR 'seq'>>>)
do
	mkdir -p lambda_$l
	echo "Starting minimization for lambda = $l ..." 

	# Define Variables
	# the following two variables are very important for the md runs. To set them correctly, of course one
	# would have to understand what they are used for. in vim you can search the file by typing e.g. '/wd_prev' + <enter>
	# you can find the next/previous instance by hitting <n>/<N> keys.
	# What are they used for? Can you set them freely or are they constrained to certain values?
	wd_prev=<<<INSERT THE CORRECT VALUE FOR THIS VARIABLE>>>
	sim_prev=<<<INSERT THE CORRECT VALUE FOR THIS VARIABLE>>>
	lambdadir=$FREE_ENERGY/lambda_$l	# current lambda value

	# loop over all different simulation steps
	for sim_phase in <<<INSERT THE CORRECT ITERABLES>>>
	do
		# Concatenate the according banner to give information about the current run (just for fun)
		cat $sim_phase.banner 
		
		# start the run of the corresponding lambda and simulation step
		gmx grompp \
			-f $MDP/${sim_phase}_${l}.mdp \
			-c $wd_prev/$sim_phase_prev.gro \
			-p $FREE_ENERGY/topol.top \
			-o $cwd/$sim_phase.tpr \
			-po $cwd/mdout.mdp &>$cwd/grompp.$sim_phase \
			|| { echo "something went wrong, check $cwd/grompp.$sim_phase"; exit; }
		gmx mdrun \
			-deffnm <<<SET THE OPTION FOR THIS FLAG YOURSELF>>> &>$cwd/mdrun.$sim_phase \
			|| { echo "something went wrong, check $cwd/mdrun.$sim_phase"; exit; }

		gro_file="$OUTDIR/$sim_phase.gro"
	<<<YOU ARE NOT done YET. SOMETHING IS MISSING HERE...>>>
done
