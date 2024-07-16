# variable definitions for the input data
fld_path=autogrid/hmx_B99990001.maps.fld
lig_path=../02_ligand_prep/PDBQT

# variable definitions for the output directories
doc_path='docking'

# create directories if not already present
mkdir -p $doc_path

# create the file by writing the location of the maps.fld file to it
echo  $fld_path > docking.dat

# cycle through all PDBQT ligands in $lig_path
for lig in $(ls $lig_path)
do
	# get the lig_id via bash parameter expansion
	lig_id=${lig##*_}
	lig_id=${lig_id%%.*}

	# write the input and output paths into docking.dat
	echo "$lig_path/$lig" >> docking.dat
	echo "$doc_path/Ligand_$lig_id" >> docking.dat
	
	# breaks the for loop directly after the first iteration
	# this is for testing; comment out if large scale screening is desired
	break
done

# run AutoDock-GPU; The path to the executable might deviate from user to user
../AutoDock-GPU/bin/autodock_gpu_64wi --filelist docking.dat
