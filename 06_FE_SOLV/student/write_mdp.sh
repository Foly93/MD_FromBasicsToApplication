#!/bin/bash

mkdir -p MDP
echo "Should we just test the simulation workflow? [y/n]"
read TEST

# for loops a very useful tools in bash and also in almost 
# every other programming language as well
# here, a variable is set iteratively to the values listed after 'in'
# These values, are to be set by you. They must match the names already in use
# in this directory as can be seen below: ${mdp}.mdp corresponds to the mdp
# files in your exercise directory. If you are unsure, ask the tutor or the internet.
for mdp in <<<INSERT THE ITERABLES HERE>>>
do
	for i in <<<INSERT THE ITERABLES HERE>>>
	do
		# What means the '-E' option of 'sed'? Add some comments
		sed -E "s/(init_lambda_state\s+=\s+)0/\1$i/" ${mdp}.mdp > ./MDP/${mdp}_${i}.mdp

		if [[ $TEST == y ]];
		then
			# what does this command do? Add some comments
			sed -i 's/nsteps.*/nsteps = 500/g' ./MDP/${mdp}_${i}.mdp
		fi
	done
done
