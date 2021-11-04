#!/bin/bash


module load gcc/4.9.2
module load open_mpi/1.6.5


bsub cmake .
bsub make

i=10000 # iterations
for f in {0..5..1}
do
    for x in {1..50} ;
    do
	bsub ./ps_functions $i $f
    done
done
for t in {2..33..1}#2 3 4 5 7 8 9 15 16 17 31 32 33 #{2..16..1}  #number of threads
do
    export OMP_NUM_THREADS=$t
    for f in in {0..5..1} # funciton id
    do
        for x in {1..50} ;
        do
            bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./psps_functions $i $f           
        done
    done
done

exit 0