#!/bin/bash
mkdir bin
cd bin

echo Loading intel compilers...
module load intel/compilers/11.1.059

echo Compiling fst...
ifort -c ../fst.f -o fst.o -O3

echo Compiling and linking poisson-serial...
icc -c ../poisson-serial.c -o poisson-serial.o -O3 -std=c99
icc poisson-serial.o fst.o -o poisson-serial

echo Compiling and linking poisson-mpi...
mpicc -c ../poisson-mpi.c -o poisson-mpi.o -O3 -std=c99 -lm
mpicc poisson-mpi.o fst.o -o poisson-mpi -lm

echo Compiling and linking poisson-hybrid...
mpicc -c ../poisson-hybrid.c -o poisson-hybrid.o -openmp -O3 -std=c99 -lm
mpicc poisson-hybrid.o fst.o -o poisson-hybrid -openmp -lm

cd ..
echo Done compiling!
