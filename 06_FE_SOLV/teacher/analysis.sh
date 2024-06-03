#!/bin/bash

module load gromacs/5.1.4

gmx_mpi bar -f _output/lambda_*/prod.xvg -o -oi
