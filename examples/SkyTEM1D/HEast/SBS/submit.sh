#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l walltime=01:00:00
#PBS -l wd
#PBS -N juliaAEM
#PBS -e grid.err

/home/anray/Software/transD_GP/examples/SkyTEM1D/HEast/SBS
module load openmpi
/home/547/ar0754/bin/julia parallelschema_2.jl >& outfile
