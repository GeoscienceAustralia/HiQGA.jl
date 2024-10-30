#!/bin/bash
#PBS -P ns59
#PBS -q normalsr
#PBS -l ncpus=208
#PBS -l mem=1000GB
#PBS -l walltime=00:40:00
#PBS -l wd
#PBS -N juliaSkyTEM_grad
#PBS -e grid.err
#PBS -l storage=scratch/ns59
module load intel-mpi/2021.10.0
mpiexec -np 208 $HOME/bin/julia 03_parallel_run.jl >& outfile
