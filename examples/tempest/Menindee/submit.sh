#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=960
#PBS -l mem=3840GB
#PBS -l walltime=16:00:00
#PBS -l wd
#PBS -N juliaAEM
#PBS -e grid.err
#PBS -l storage=gdata/z67+scratch/z67
module load intel-mpi/2019.8.254
julia --project -e 'using Pkg; Pkg.precompile()'
mpirun -np 960 /home/547/ar0754/bin/julia-1.6.0 03_parallel_run.jl >& outfile
