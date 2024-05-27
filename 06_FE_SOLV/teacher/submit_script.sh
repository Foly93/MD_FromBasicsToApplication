#!/bin/bash

mkdir -p _trash
mkdir -p _output

for l in $(seq 0 20)
do
	mkdir -p _output/lambda_$l
	echo "Starting minimization for lambda = $l ..." 

	# loop over all different simulation steps
	sbatch \
	--job-name=06_FESOLV_lambda${l} \
	--nodes=1 \
	--gres=gpu:1 \
	--time=1-00:00:00 \
	--ntasks=4 \
	--cpus-per-task=1 \
	--output=_trash/lambda_${l}_slurm.out \
	--error=_trash/error_lambda_${l} \
	./run_lambda_simulations.sh $l
done
