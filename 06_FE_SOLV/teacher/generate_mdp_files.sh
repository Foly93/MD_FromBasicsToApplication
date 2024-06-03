#!/bin/bash

echo "Is this just a test run? (recommended for first try)? [y/n]"
read TEST

for mdp in em nvt npt prod
do
	for i in $(seq 0 20)
	do
		# '-E' enables regular expression in the sed pattern between '/../../' and groups with '()'
		sed -E "s/(init.lambda.state\s+=\s+)0/\1$i/" ${mdp}.mdp > ${mdp}_${i}.mdp

		if [[ $TEST == y ]];
		then
			# changes the step number to 500 and therefore the runs are done very quickly
			sed -i 's/nsteps.*/nsteps = 500/g' ${mdp}_${i}.mdp
		fi
	done
done
