#!/bin/bash
#PBS -P ns59
#PBS -q normalsr
#PBS -l ncpus=312
#PBS -l mem=1500GB
#PBS -l walltime=0:20:00
#PBS -l wd
#PBS -N SkyTEM_plot
#PBS -e 0000_grid.err
#PBS -l storage=gdata/z67+scratch/z67+gdata/cr78+gdata/qi71
module load intel-mpi/2021.10.0
mpiexec -np 312 $HOME/bin/julia plot_parallel.jl >& outfile_plot
