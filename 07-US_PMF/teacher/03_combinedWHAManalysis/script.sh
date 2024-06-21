#!/bin/bash

# One run was conducted with GROMACS 5.1.4 and the other with GROMACS 2018.3; Newer versions can read older versions but not vice versa!
ml gromacs/2018.8

gmx_mpi wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
