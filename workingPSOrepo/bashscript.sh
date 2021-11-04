#!/bin/bash
module load gcc/4.9.2
module load open_mpi/1.6.5
module load gsl/1.16
module load interconnect/ethernet
module load mvapich2
#bsub mpiCC  -L/cluster/apps/gsl/1.16/x86_64/gcc_4.8.2/include -std=c++14 -Wall -pedantic -O3 -march=native -o ppso_alternative copympipso.cpp -lgsl -lgslcblas -lm
export OMP_NUM_THREADS=4
#bsub -R "span[ptile=4]" -n 4 mpirun -np 4 ./ppso_alternative 10000 0 1 100
# cmake .
# make
