#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=12
#PBS -l mem=64GB
#PBS -l walltime=01:30:00
#PBS -l wd
#PBS -N juliaAEM
#PBS -e grid.err
#PBS -l storage=gdata/z67+scratch/z67+gdata/kb5
module load intel-mpi/2019.8.254
~/bin/julia --project -e 'using Pkg; Pkg.precompile()'
mpirun -np 12 ~/bin/julia 03_parallel_run_nci.jl >& outfile
