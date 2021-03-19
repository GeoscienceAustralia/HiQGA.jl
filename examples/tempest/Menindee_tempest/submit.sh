#!/bin/bash
#PBS -P z67
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l walltime=00:20:00
#PBS -l wd
#PBS -N juliaAEM
#PBS -e grid.err
#PBS -l storage=gdata/z67
cd /scratch/z67/ar/transD_GP/examples/tempest/Menindee_tempest/
module load openmpi/4.0.1
mpirun -np 48 /home/547/ar0754/bin/julia doall.jl >& outfile
