#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=720
#PBS -l mem=2880GB
#PBS -l walltime=1:30:00
#PBS -l wd
#PBS -N juliaSkyTEM_grad
#PBS -e grid.err
#PBS -l storage=gdata/z67+scratch/z67+gdata/cr78
module load intel-mpi/2019.8.254
~/bin/julia --project -e 'using Pkg; Pkg.precompile()'
mpiexec -np 720 ~/bin/julia 03_parallel_run.jl >& outfile
