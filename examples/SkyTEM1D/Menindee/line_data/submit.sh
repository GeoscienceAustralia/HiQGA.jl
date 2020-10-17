#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=240
#PBS -l mem=960GB
#PBS -l walltime=7:00:00
#PBS -l wd
#PBS -N juliaAEM
#PBS -e grid.err

cd /scratch/z67/ar/transD_GP/examples/SkyTEM1D/Menindee/run2
module load openmpi
mpirun -np 240 /home/547/ar0754/bin/julia doall.jl >& outfile
