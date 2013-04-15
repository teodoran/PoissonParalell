#!/bin/bash
mkdir bin
cd bin

echo Compiling fst...
gfortran -c ../fst.f -o fst.o

echo Compiling and linking poisson-serial...
gcc -std=c99 -c ../poisson-serial.c -o poisson-serial.o
gfortran poisson-serial.o fst.o -o poisson-serial

echo Compiling and linking poisson-mpi...
mpicc -std=c99 -lm -c ../poisson-mpi.c -o poisson-mpi.o
mpicc -lm poisson-mpi.o fst.o -o poisson-mpi

echo Compiling and linking poisson-hybrid...
mpicc -fopenmp -std=c99 -lm -c ../poisson-hybrid.c -o poisson-hybrid.o
mpicc -fopenmp -lm poisson-hybrid.o fst.o -o poisson-hybrid

cd ..
echo Done compiling!