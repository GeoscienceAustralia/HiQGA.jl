#!/bin/bash
#PBS -P cr78
#PBS -q normal
#PBS -l ncpus=912
#PBS -l mem=3648GB
#PBS -l walltime=00:50:00
#PBS -l wd
#PBS -N juliaTEMPEST_grad
#PBS -e grid.err
#PBS -l storage=gdata/z67+scratch/z67+gdata/cr78
module load intel-mpi/2019.8.254
~/bin/julia --project -e 'using Pkg; Pkg.precompile()'
mpiexec -np 912 ~/bin/julia 03_parallel_run.jl >& outfile
