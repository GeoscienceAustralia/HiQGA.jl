#!/bin/bash
#PBS -P cr78
#PBS -q normal
#PBS -l ncpus=240
#PBS -l mem=960GB
#PBS -l walltime=1:00:00
#PBS -l wd
#PBS -N juliaSkyTEM_grad
#PBS -e grid.err
#PBS -l storage=gdata/z67+scratch/z67
module load intel-mpi/2019.8.254
~/bin/julia --project -e 'using Pkg; Pkg.precompile()'
mpiexec -np 240 ~/bin/julia 03_parallel_run.jl >& outfile

