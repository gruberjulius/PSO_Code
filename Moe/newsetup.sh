#!/bin/bash
#HERE still needed from pso codes:
#  std::ofstream times("times_papso.dat", std::ios::app); 
#   times << diff.count() << " " << [bestfitnessvalue] HERE << "\n";
#    times.close();

#times_papso.dat
#times_pspso.dat
#times_pso.dat


#first building
mkdir build
cd build
cmake ..
make

#now computing
RESULT_PAPSO_F=results_papso.dat
RESULT_PSPSO_F=results_pspso.dat
RESULT_PSO_F=results_pso.dat
i=1000 # iterations
numcalls=10 #HERE replace 10 by numcalls #number of calls to get good mean
for N in 1 #{1..5..1}  #number of functions
do
  FILE=func${N}
  rm -r $FILE
  mkdir $FILE
  cd $FILE
  rm $RESULT_PAPSO_F
  rm $RESULT_PSPSO_F
  rm $RESULT_PSO_F
  touch $RESULT_PAPSO_F
  touch $RESULT_PSPSO_F
  touch $RESULT_PSO_F
  for F in {1..10..1} #number of function calls for good mean
    do
      ./../ps_functions $i   #HERE implement a way to say which functions we want to evaluate
      echo "${N} 1 " >> RESULT_PSO_F
      times_pso.dat >> RESULT_PSO_F
      rm times_pso.dat
    done
  for T in {2..16..1}  #number of threads
  do
    for F in {1..10..1} #number of function calls for good mean
    do
      mpirun -np $T ../paps_functions $i   #HERE implement a way to say which functions we want to evaluate
      mpirun -np $T ../psps_functions $i   #HERE implement a way to say which functions we want to evaluate
      echo "${N} ${T} " >> RESULT_PAPSO_F
      times_papso.dat >> RESULT_PAPSO_F
      echo "${N} ${T} " >> RESULT_PSPSO_F
      times_pspso.dat >> RESULT_PSPSO_F
      rm times_papso.dat times_pspso.dat
    done
  done
  cd ..
done

