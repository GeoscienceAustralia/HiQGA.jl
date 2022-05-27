#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=480
#PBS -l mem=1920GB
#PBS -l walltime=1:00:00
#PBS -l wd
#PBS -N juliaSkyTEM_grad
#PBS -e grid.err
#PBS -l storage=gdata/z67+scratch/z67
module load intel-mpi/2019.8.254
~/bin/julia --project -e 'using Pkg; Pkg.precompile()'
mpiexec -np 480 ~/bin/julia 03_parallel_run.jl >& outfile
