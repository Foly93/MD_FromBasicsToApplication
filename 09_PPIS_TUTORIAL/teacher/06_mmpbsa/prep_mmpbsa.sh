#!/bin/bash

L_ids=(0156 1258)
R_id='3hmx'
reps=3

A_res='1-306'
B_res='307-503'
L_res='504'

ABL=( 'A' 'B' 'L' 'AL' 'BL')
ABLres=( $A_res $B_res $L_res $A_res,$L_res $B_res,$L_res )

for L_id in ${L_ids[@]}
do
	R_L=${R_id}_$L_id
	SIM=../05_mdsimulation/_sim/$R_L

	mkdir -p -v $R_L
	cp -v mmpbsa.sh $R_L
	chmod 777 $R_L/mmpbsa.sh
	

	for i in $(seq 0 4)
	do
		echo "parmstrip !(:${ABLres[$i]})" > cpptraj
		echo "parmwrite out $R_L/${ABL[$i]}.prmtop" >> cpptraj
		cpptraj -i cpptraj -p $SIM/$R_L.prmtop
	done
	
	for rep in $(seq 1 $reps)
	do
		for i in $(seq 3 4)
		do
			echo "strip !(:${ABLres[$i]})" > cpptraj
			echo "trajout  $R_L/${ABL[$i]}_$rep.nc" >> cpptraj
			cpptraj -i cpptraj -p $SIM/$R_L.prmtop -y $SIM/prod_center_$rep.nc
		done
	done
done
