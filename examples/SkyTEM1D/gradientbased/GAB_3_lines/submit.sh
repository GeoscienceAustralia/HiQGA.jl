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
/home/547/ar0754/bin/julia-1.6.0 --project -e 'using Pkg; Pkg.precompile()'
/g/data/z67/ar0754/ompi-4.1.0/bin/mpiexec -np 480 /home/547/ar0754/bin/julia-1.6.0 03_parallel_run.jl >& outfile
