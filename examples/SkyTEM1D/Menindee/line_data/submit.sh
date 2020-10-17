#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=1200
#PBS -l mem=4800GB
#PBS -l walltime=12:00:00
#PBS -l wd
#PBS -N juliaAEM
#PBS -e grid.err

cd /scratch/z67/ar/transD_GP/examples/SkyTEM1D/Menindee/line_data
module load openmpi
mpirun -np 1200 /home/547/ar0754/bin/julia doall.jl >& outfile
