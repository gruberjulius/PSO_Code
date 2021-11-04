#!/bin/bash

module load gcc/4.9.2
module load open_mpi/1.6.5


i=10000 # iterations

for t in 2 3 4 8 16 #{2..16..1}  #number of threads
do
    export OMP_NUM_THREADS=$t
    for f in in {0..5..1} # funciton id
    do
        for x in {1..50} ;
        do
            bsub -R "span[ptile=$t]" -n $t mpirun -np $t ./paps_functions_send_move_smaller_packets $i $f
            bsub -R "span[ptile=$t]" -n $t mpirun -np $t ./paps_functions_send_move_part_unroll $i $f
            bsub -R "span[ptile=$t]" -n $t mpirun -np $t ./paps_functions_send_move $i $f
        done
    done
done

exit 0
