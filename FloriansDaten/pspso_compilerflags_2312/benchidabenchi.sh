#!/bin/bash
#
#  std::ofstream times("times_paps.dat", std::ios::app);
#   times << iterations << ';' <<diff.count() << "\n";
#    times.close();
#module load gcc/4.9.2
#module load open_mpi/1.6.5
#bsub cmake ..
#bsub make
#i=10000 # iterations
for i in 10000
do
    for f in {0..5..1}
    do
        for x in {1..50}
        do
        bsub ./ps_functions $i $f
        done
    done
    for t in {1..65..1}
    do
        export OMP_NUM_THREADS=$t
        for f in in {0..5..1} # funciton id
        do
            for x in {1..50}
            do
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./paps_functions_send_move $i $f
                bsub -R "span[ptile=$t]" -n $t mpirun --mca btl_vader_single_copy_mechanism none -np $t ./psps_functions_kristof $i $f
            done
        done
    done
done
exit 0

