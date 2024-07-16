#!/bin/bash
receptor='./hmx_B99990001.pdbqt'
TEST='../02*/TEST'
PDBQT='./PDBQT'
PDB='./PDB'

mkdir -p -v $PDBQT $PDB

# Loop through all ligands
for file in $TEST/*
do
	# Fetch ligand basename
	name="${file##*/}"
	lig="${name%.pdbqt}"
	# Docking with Autodock vina
	vina --config config \
		--receptor $receptor \
		--ligand $TEST/$lig.pdbqt \
		--out $PDBQT/out_$lig.pdbqt

	# Transforming and splitting the PDBQT output to PDB format
	obabel -ipdbqt $PDBQT/out_$lig.pdbqt \
		-opdb -m -O $PDB/out_${lig}_.pdb
done
# Convert pdbqt to pdb via http://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html



