#!/bin/bash

#PBS -P z67
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=256GB
#PBS -l walltime=00:15:00
#PBS -l wd
#PBS -N testJulia
#PBS -o grid.out
#PBS -e grid.err
#PBS -j oe

ulimit -s unlimited
ulimit -c unlimited
module load gcc/5.2.0 openmpi/3.0.1 
mpirun  -np 1 /g/data1a/zk34/julia-1.1.1/julia ./run_image_regression.jl > outfile.run

