
module load gcc/4.9.2
module load open_mpi/1.6.5



i=10000 # iterations
for t in 2 4 8 16 #{2..16..1}  #number of threads
do
    export OMP_NUM_THREADS=$t
    for x in {1..50}
    do
          bsub -R "span[ptile=$t]" -n $t mpirun -np $t ./paps_functions_precompute_vel $i 0
          bsub -R "span[ptile=$t]" -n $t mpirun -np $t ./paps_functions_precompute_vel_nonblocking $i 0
    done
done

exit 0
