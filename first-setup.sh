#!/bin/bash

# We want a C compiler
sudo apt-get install gcc

# We want a Fortran compiler
sudo apt-get install gfortran

# We want BLAS
sudo apt-get install libblas3gf libblas-doc libblas-dev

# We want LAPAC
sudo apt-get install liblapack3gf liblapack-doc liblapack-dev

# We also need mpi. Add some yourself.

# We want to build the project
./build.sh