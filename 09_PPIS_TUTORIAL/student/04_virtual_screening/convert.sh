#!/bin/bash

pdbqt='./PDBQT'
pdb='./PDB'
mol2='./MOL2'
smi='./SMI'

mkdir -p $pdb $mol2 $smi

for file in $pdbqt/*
do
	# get basename via parameter expansion
	name="${file##*/}"
	lig="${name%.pdbqt}"
	# converting to pdb
	echo "$pdbqt/$lig"
	obabel -ipdbqt $pdbqt/$lig.pdbqt -opdb -O $pdb/$lig.pdb
	# converting to mol2
	obabel -ipdbqt $pdbqt/$lig.pdbqt -omol2 -O $mol2/$lig.mol2 -p 7.4
	# converting to smi
	obabel -ipdbqt $pdbqt/$lig.pdbqt -osmi -O $smi/$lig.smi
done
