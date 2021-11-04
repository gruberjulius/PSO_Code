#!/bin/bash

module load gcc/4.9.2
module load open_mpi/1.6.5
module load gsl/1.6
module load interconnect/ethernet
module load mvapich2/2.1
module load hwloc/1.11.0



i=10000 # iterations
b=1 # basic file
p=100 # particles
for t in {1..33..1}
do
    export OMP_NUM_THREADS=$t
    for f in {0..5..1} # funciton id
    do
        for x in {1..50} ;
        do
		bsub -R "span[ptile=$t]" -n $t mpirun -np $t ./hello $i $f $b $p
        done
    done
done

exit 0
