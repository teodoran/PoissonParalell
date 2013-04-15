#!/bin/bash
./build.sh
cd bin

echo
echo Collecting system info
nproc
grep -i memtotal /proc/meminfo
grep -i memfree /proc/meminfo

for N in 32 64 128 256 512 #1024 2048 4096
do
	echo
	echo Calculating with problemsize $N
	echo Serial:
	./poisson-serial $N
	echo MPI:
	mpirun -np 4 poisson-mpi $N
	echo Hybrid:
	mpirun -np 4 poisson-hybrid $N
done

cd ..
echo
echo Done!

