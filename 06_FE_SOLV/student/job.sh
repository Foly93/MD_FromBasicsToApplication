#!/bin/bash

# Set some environment variables 
FREE_ENERGY=$(pwd)
MDP=$FREE_ENERGY/MDP

# make the user aware of his/her choices
echo "Free energy home directory set to $FREE_ENERGY"
echo ".mdp files are stored in $MDP"

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
	for sim in <<<INSERT THE CORRECT ITERABLES>>>
	do
		# preliminaries before the simulations are executed
		cwd="$lambdadir/$sim"		# set the current working directory
		cat $sim.banner			# Concatenate the according banner to give information about the current run (just for fun) 
		mkdir -p $cwd			# create the current working directory
		
		# start the run of the corresponding lambda and simulation step
		gmx grompp -f $MDP/${sim}_${l}.mdp -c $wd_prev/$sim_prev.gro -p $FREE_ENERGY/topol.top -o $cwd/$sim.tpr -po $cwd/mdout.mdp &>$cwd/grompp.$sim || { echo "something went wrong, check $cwd/grompp.$sim"; exit; }
		gmx mdrun -deffnm <<<SET THE OPTION FOR THIS FLAG YOURSELF>>> &>$cwd/mdrun.$sim || { echo "something went wrong, check $cwd/mdrun.$sim"; exit; }

		# save the name of the current simulation step to access it in the next iteration 
		sim_prev=$sim			# name of next simulations input file
		wd_prev="$lambdadir/$sim"	# name of next simulations input directory
	<<<YOU ARE NOT done YET. SOMETHING IS MISSING HERE...>>>
done
