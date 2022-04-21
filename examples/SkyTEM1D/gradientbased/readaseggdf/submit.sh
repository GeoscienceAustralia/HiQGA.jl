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
# if there are race conditions run next line with active internet connection
# /g/data4/z67/ar0754/julia-1.6.0/bin/julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
/g/data/z67/ar0754/ompi-4.1.0/bin/mpiexec -np 240 /g/data4/z67/ar0754/julia-1.6.0/bin/julia 03_parallel_run.jl >& outfile

