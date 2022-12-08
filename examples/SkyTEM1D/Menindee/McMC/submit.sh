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
module load intel-mpi/2019.8.254
~/bin/julia --project -e 'using Pkg; Pkg.precompile()'
mpiexec -np 960 ~/bin/julia 03_parallel_run.jl >& outfile
