#!/bin/bash

# declare identifier of ligand and receptor and the pdb file created during the receptor preparation
R_id='3hmx'
L_ids=(0156 1258)
R_pdb='../../01_receptor_prep/hmx.B99990001.pdb'

# declare variables of the respective directories
INP='../_inp'
SIM='../_sim'
MOL2='../../04_virtual_screening/MOL2'

# copy the receptor pdb into the input directory for leaping it later for parameter preparation
cp -v $R_pdb $INP/$R_id.pdb

# loop over the ligand ids specified above
for L_id in "${L_ids[@]}"
do
	# declare a combined variable of receptor and ligand for easier readability
	R_L=${R_id}_${L_id}
	# create the simulation directories for each ligand automatically (more consistent)
	mkdir -v -p $SIM/$R_L
	# copy the executables to the respective sim directory and substitute the correct IDs
	cp -v prep_and_run.sh $SIM/$R_L/
	cp -v convert_gmx2amb.py $SIM/$R_L/
	cp -v amb2gro_top_gro_cmap.py $SIM/$R_L/
	sed -i "s/LIGANDID/$L_id/g" $SIM/$R_L/prep_and_run.sh
	sed -i "s/RECEPTORID/$R_id/g" $SIM/$R_L/prep_and_run.sh
	chmod 777 $SIM/$R_L/prep_and_run.sh
	# copy the mol2 file of the each ligand here for manual modification (overwriting is deactivated via -n)
	cp -v -n $MOL2/Ligand_$L_id.mol2 $INP/$L_id.mol2
done
