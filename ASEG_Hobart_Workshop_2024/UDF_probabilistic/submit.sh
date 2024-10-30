#!/bin/bash
#PBS -P ns59
#PBS -q normalsr
#PBS -l ncpus=8320
#PBS -l mem=40000GB
#PBS -l walltime=04:00:00
#PBS -l wd
#PBS -N juliaSkyTEM_prob
#PBS -e grid.err
#PBS -l storage=scratch/ns59
module load intel-mpi/2021.10.0
mpiexec -np 8304 $HOME/bin/julia 03_parallel_run.jl >& outfile
