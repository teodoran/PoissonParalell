#!/bin/bash                                                                                                                                                                                                          

#PBS -N gt-poissonHyb                                                                                                                                                                                                

# Allocate three nodes with 36 processors                                                                                                                                       
#PBS -lnodes=12:ppn=3:default                                                                                                                                                                                        

# Expect to run up to 5 minutes                                                                                                                                                                                      
#PBS -lwalltime=0:30:00                                                                                                                                                                                              

# Memory per process                                                                                                                                                                                                 
#PBS -lpmem=2000MB                                                                                                                                                                                                   

# Run on the freecycle account                                                                                                                                                                                       
#PBS -A freecycle                                                                                                                                                                                                    

# Run in the optimist queue by default                                                                                                                                                                               
#PBS -q optimist                                                                                                                                                                                                     

# Join stdout and stderr output to one file                                                                                                                                                                          
#PBS -j oe                                                                                                                                                                                                           

#PBS -m abe                                                                                                                                                                                                          
#PBS -M georgekw@stud.ntnu.no                                                                                                                                                                                        

# Change directory to dir with the job script                                                                                                                                                                        
cd ${PBS_O_WORKDIR}

# Load needed modules                                                                                                                                                                                                
module load intelcomp
module load openmpi/1.4.3-intel

# Set thread affinity                                                                                                                                                                                                
KMP_AFFINITY="granularity=fine,compact"


# Run with 12 MPI processes, each with 3 threads                                                                                                                                                                     
OMP_NUM_THREADS=3 mpirun -npernode 4 poisson-hybrid 16384
