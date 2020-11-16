#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=240
#PBS -l mem=960GB
#PBS -l walltime=02:00:00
#PBS -l wd
#PBS -N juliaAEM
#PBS -e grid.err
#PBS -l storage=gdata/z67+gdata/r78

module load openmpi
mpirun -np 240 julia doall.jl >& outfile
