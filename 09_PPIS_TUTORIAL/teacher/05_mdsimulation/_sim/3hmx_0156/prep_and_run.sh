#!/bin/bash

# declare identifier of ligand+receptor (done automatically by prep_setup.sh) and the path of input dir
INP='../../_inp'
L_id='0156'
R_id='3hmx'
R_L=${R_id}_${L_id}

cp_inp_files()
{
	# copy all necessary input file into the simulation dir, with verbosity turned on
	cp -v $INP/*.mdp ./
	cp -v $INP/$L_id.mol2 ./
	cp -v $INP/$R_id.pdb ./$R_L.pdb
	cp -v $CONDA_PREFIX/dat/leap/parm/gaff2.dat ./

	# change and copy the leap files into the simulation dir and rename them
	sed "s/LID/$L_id/g" $INP/ligand_leap.in > ./${L_id}_leap
	sed "s/LID/$L_id/g" $INP/complex_leap > ./${R_L}_leap
	sed -i "s/RL_ID/$R_L/g" ./${R_L}_leap
}


gen_params()
{
	# antechamber and parmchk2 are used to create prameters for the ligands
	antechamber -i $L_id.mol2 -at GAFF2 -c bcc -fi mol2 -fo mol2 -o ${L_id}_gaff2.mol2
	parmchk2 -i ${L_id}_gaff2.mol2 -f mol2 -o ${L_id}_gaff2.frcmod

	# leap the ligand alone, to generate .lib .pdb and .prmtop file
	tleap -f ${L_id}_leap
	# exchange some artifacts from tleap with the proper name (LIG)
	sed -i 's/MOL/LIG/g' $L_id.pdb
	sed -i 's/UNK/LIG/g' $L_id.pdb 
	sed -i 's/UNL/LIG/g' $L_id.pdb 
	sed -i 's/ACE/LIG/g' $L_id.pdb 
	sed -i 's/.00/LIG/g' $L_id.pdb 

	# add the ligand coordinates to the complex pdb (receptor and ligand) and correct some keywords
	cat $L_id.pdb >> $R_L.pdb
	sed -i '/END/d' $R_L.pdb
	sed -i '/CONECT/d' $R_L.pdb
	
	# leap the complex pdb and create amber run input files
	tleap -f ${R_L}_leap
}

gen_gmxtop()
{
	# amber cannot be used for running the simulations since it is not open source :(
	# Remove prior topology and prepare the new one GROMACS topology with awesome external python script.
	mv -v $R_L.gro $R_L.top ../../_trash
	printf "y\n$R_L.prmtop\n$R_L.inpcrd\n${L_id}_gaff2.frcmod\ngaff2.dat\n" | python convertTOP_amb2gmx.py $R_L
	cp $R_L.top topol.top
	
	# add forcefield specifications for water, ions and general parameters, also delete defaults section
	sed -i '/\[ defaults \]/i \#include "oplsaa.ff/forcefield.itp"' topol.top
	sed -i '/\[ defaults \]/,+3d' topol.top
	sed -i '/\[ system \]/i \#include "oplsaa.ff/spce.itp"' topol.top
	sed -i '/\[ system \]/i \#include "oplsaa.ff/ions.itp"' topol.top
}

run_gmx_sims()
{
	### students will not need 'gmx_mpi' but 'gmx'
	### not really necessary since the students have gromacs installed on  their pc's
	ml gromacs/5.1.4
	# solvate the system and generate the ions
	gmx_mpi solvate -cp $R_L.gro -cs spc216.gro -o ${R_L}_solv.gro -p topol.top -box 13 13 13
	gmx_mpi grompp -f em.mdp -c ${R_L}_solv.gro -p topol.top -o ions.tpr
	printf "SOL\n" | gmx_mpi genion -s ions.tpr -o ${R_L}_ready.gro -p topol.top -pname NA -nname CL -conc 0.1 -neutral

	# loop over repetitions 1 to 3
	for rep in 1 2 3
	do
		# set j to the initial gro file of each repetition
		j="${R_L}_ready"
		# loop over the different equilibrations and production simulations
		for i in em nvt npt pre prod
		do
			# standard gromacs mdrun execution routine
			gmx_mpi grompp -f $i.mdp -c $j.gro -p topol.top -o ${i}_$rep.tpr
			gmx_mpi mdrun -ntomp 16 -deffnm ${i}_$rep -v
			# re-set j to the current file prefix for the input file of the next iteration
			j=${i}_$rep
		done
		# center the ligand in the trajectory system 
		printf "13\n0\n" | gmx_mpi trjconv -s prod_$rep.tpr -f prod_$rep.xtc -o prod_center_$rep.xtc -pbc mol -center

		echo "convert Trajectory from GROMACS XTC to AMBER NetCDF4"
		sed "s/REP/$rep/g" convert_gmx2amb.py > convert${rep}_gmx2amb.py
		python convert${rep}_gmx2amb.py
	done
}

# functions can be executed independently; switch their execution on/off by inserting/deleting '#'
cp_inp_files
gen_params
#gen_gmxtop
#run_gmx_sims
