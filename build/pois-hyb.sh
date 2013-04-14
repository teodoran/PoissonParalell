#!/bin/bash

#PBS -A acc-freecycle

#PBS -lnodes=3:ppn=2:default

#PBS -lwalltime=00:30:00

#PBS -lpmem=2000MB

#PBS -j oe

#PBS -m abe

#PBS -M georgekw@stud.ntnu.no

#PBS -q optimist

# Set thread affinity
KMP_AFFINITY="granularity=fine,compact"
 
cd $PBS_O_WORKDIR
OMP_NUM_THREADS=3 mpirun -npernode 4 poisson-hybrid 16384