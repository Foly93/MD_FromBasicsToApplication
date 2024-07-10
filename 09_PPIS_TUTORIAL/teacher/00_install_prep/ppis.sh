#!/bin/bash

$MAMBA_EXE create -n ppis python=3.8
$MAMBA_EXE install -n ppis -c conda-forge -c schrodinger pymol-bundle -y
$MAMBA_EXE install -n ppis -c conda-forge ambertools=22.0 -y
$MAMBA_EXE install -n ppis -c bioconda autodock-vina=1.1.2 -y
$MAMBA_EXE install -n ppis -c bioconda autogrid=4.2.6 -y
$MAMBA_EXE install -n ppis -c conda-forge fpocket=4.0.2 -y
$MAMBA_EXE install -n ppis -c conda-forge ipykernel=6.29.0 -y
$MAMBA_EXE install -n ppis -c conda-forge ipympl=0.9.3 -y
$MAMBA_EXE install -n ppis -c conda-forge jupyter=1.0.0 -y
$MAMBA_EXE install -n ppis -c conda-forge mdanalysis=2.4.3 -y
$MAMBA_EXE install -n ppis -c conda-forge mdtraj=1.9.9 -y
$MAMBA_EXE install -n ppis -c salilab modeller=10.5 -y
$MAMBA_EXE install -n ppis -c conda-forge mplcursors=0.5.3 -y
$MAMBA_EXE install -n ppis -c conda-forge natsort=8.4.0 -y
$MAMBA_EXE install -n ppis -c conda-forge nb_conda_kernels=2.3.1 -y
$MAMBA_EXE install -n ppis -c conda-forge nglview=3.1.1 -y
$MAMBA_EXE install -n ppis -c conda-forge openbabel=3.1.1 -y
$MAMBA_EXE install -n ppis -c conda-forge rdkit=2022.09.1 -y
$MAMBA_EXE install -n ppis -c hcc autodock-gpu -y
