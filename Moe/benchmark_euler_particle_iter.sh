#!/bin/bash

module load gcc/4.9.2
module load open_mpi/1.6.5


bsub cmake .
bsub make

b=0 # basic file
for i in 1000 # 2000 4000
do
  for p in 500 1000 2000
  do
    for f in {0..5..1}
    do
        for x in {1..50}
        do
      bsub ./ps_functions $i $f $b $p
        done
    done
    for t in {1..33..1}
    do
        export OMP_NUM_THREADS=$t
        for f in {0..5..1} # funciton id
        do
            for x in {1..50}
            do
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./paps_functions_send_move $i $f $b $p
            	bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./paps_functions_precompute_vel $i $f $b $p
            	bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./paps_functions_send_move_smaller_packets $i $f $b $p
            	bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./paps_functions_precompute_vel_nonblocking $i $f $b $p
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./paps_functions_send_move_part_unroll $i $f $b $p
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./psps_functions_kristof $i $f $b $p
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./psps_functions_k_test $i $f $b $p
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./psps_functions_k_Reduce_Bcast $i $f $b $p
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./psps_functions_k_AllReduce $i $f $b $p
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./psps_functions_k_AR_smaller $i $f $b $p
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./psps_functions $i $f $b $p

            done
        done
    done
  done
done
exit 0
