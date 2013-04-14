#!/bin/bash

# We want a C compiler
sudo apt-get install gcc

# We want a Fortran compiler
sudo apt-get install gfortran

# We want BLAS
sudo apt-get install libblas3gf
sudo apt-get install libblas-doc
sudo apt-get install libblas-dev

# We want LAPACK
sudo apt-get install liblapack3gf
sudo apt-get install liblapack-doc
sudo apt-get install liblapack-dev

# We want to setup CMake
rm -r build
mkdir build
cd build

cmake ..
make
