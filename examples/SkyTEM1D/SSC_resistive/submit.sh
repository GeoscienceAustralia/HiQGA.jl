#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l walltime=01:40:00
#PBS -l wd
#PBS -N juliaAEM
#PBS -e grid.err

cd /scratch/z67/ar/transD_GP/examples/SkyTEM1D/SSC_resistive
module load openmpi
mpirun -np 4 /home/547/ar0754/bin/julia fullscript.jl >& outfile
