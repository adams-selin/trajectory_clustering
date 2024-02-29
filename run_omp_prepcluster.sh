#!/bin/bash

#
#  Example submission script for prep_cluster (using openMP) on NCAR's derecho
#
#  Notes:
#   Assumed prep_cluster.f90 is compiled with the -qopenmp flag
#


#-------------------------------------------
# PBS stuff below this line

# job name:
#PBS -N prep_cluster_50z


# project code:
# (you must specify the project code to be charged when you submit a job)
#PBS -A UATE0004


# number of omp threads
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=8


# maximum wall-clock time (hh:mm:ss)
#PBS -l walltime=2:00:00 

# queue:
#PBS -q main

#-------------------------------------------
# command-line stuff below this line
# (probably should not change)

# temp directory:
#export TMPDIR=/glade/derecho/scratch/$USER/20120529/rst_000025
#mkdir -p $TMPDIR

mpiexec -d 8 ./prep_cluster_omp.exe '/glade/derecho/scratch/radams/20120529_zshape/subtraj' '20120529_zshape_ge50lt150_Wrel' >& prep_cluster50z.out


