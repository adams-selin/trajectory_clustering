#!/bin/bash

#
#  Example submission script for prep_cluster (without using openMP) on NCAR's derecho
# . Good if you have the distance.bin file already created, for example.
#
#  Notes:
#   Assumed prep_cluster.f90 is NOT compiled with the -qopenmp flag
#


#-------------------------------------------
# PBS stuff below this line

# job name:
#PBS -N prep_cluster25z


# project code:
# (you must specify the project code to be charged when you submit a job)
#PBS -A UATE0004


# number of omp threads
#PBS -l select=1:ncpus=2:mpiprocs=1


# maximum wall-clock time (hh:mm:ss)
#PBS -l walltime=3:00:00 

# queue:
#PBS -q main

#-------------------------------------------
# command-line stuff below this line
# (probably should not change)

./prep_cluster.exe '/glade/derecho/scratch/radams/20120529_zshape/subtraj' '20120529_zshape_ge25lt45_Wrel' >& prep_cluster25z.out


