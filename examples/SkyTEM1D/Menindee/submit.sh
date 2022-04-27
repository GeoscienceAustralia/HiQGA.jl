#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=960
#PBS -l mem=3840GB
#PBS -l walltime=12:00:00
#PBS -l wd
#PBS -N juliaSkyTEM
#PBS -e grid.err
#PBS -l storage=gdata/z67
# if there are race conditions run next line with active internet connection
#/home/547/ar0754/bin/julia-1.6.0 --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
/g/data/z67/ar0754/ompi-4.1.0/bin/mpiexec -np 960 /home/547/ar0754/bin/julia-1.6.0 03_parallel_run.jl >& outfile
