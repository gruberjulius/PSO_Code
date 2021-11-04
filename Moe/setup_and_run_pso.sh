#!/bin/bash
#HERE still needed from pso codes:
#  std::ofstream times("times_papso.dat", std::ios::app); 
#   times << diff.count() << " " << [bestfitnessvalue] HERE << "\n";
#    times.close();

#times_papso.dat
#times_pspso.dat
#times_pso.dat

module load gcc/4.9.2
module load open_mpi/1.6.5
#first building
rm -r build
mkdir build
cd build
module load gcc/4.9.2
module load open_mpi/1.6.5
bsub "cmake .. ; make"
mkdir results

#now computing

i=10000 # iterations
#numcalls=10 #HERE replace 10 by numcalls #number of calls to get good mean
for N in {0..5..1}  #number of functions
do
  for F in {1..1..1} #number of function calls for good mean
    do
      bsub ./ps_functions $i $N
    done
  for T in 2 3 #4 5 7 8 9 15 16 17 31 32 33 #{2..16..1}  #number of threads
  do
    export OMP_NUM_THREAD=$T
    for F in {1..1..1} #number of function calls for good mean
    do
      bsub -R "span[ptile=$T]" -n $T mpirun ./paps_functions $i $N
      bsub -R "span[ptile=$T]" -n $T mpirun ./psps_functions $i $N
      bsub -R "span[ptile=$T]" -n $T mpirun ./paps_functions_precompute_vel $i $N
      bsub -R "span[ptile=$T]" -n $T mpirun ./paps_functions_send_move $i $N
    done
  done
  cd ..
done
