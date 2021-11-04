#!/bin/bash
#ab jz nur fp hier
#
#  std::ofstream times("times_paps.dat", std::ios::app); 
#   times << iterations << ';' <<diff.count() << "\n";
#    times.close();
#
#
#
#
#if [[ -f results/times_ps.dat ]]
#then
#    rm results/times_ps.dat
#fi
#if [[ -f results/times_psps.dat ]]
#then
#    rm results/times_psps.dat
#fi
#if [[ -f results/times_paps.dat ]]
#then
#    rm results/times_paps.dat
#fi
#
#touch results/times_ps.dat
#touch results/times_psps.dat
#touch results/times_paps.dat

"""t=2
FOLDER=results_${t}_threads
echo "$FOLDER"
rm -r $FOLDER
mkdir $FOLDER
cd $FOLDER"""

for t in 2 3 4
do
  for i in  1000 10000 # {1000..20000..1000} #1000 10000 100000 1000000
  do
    for j in {1..100}
    do
      #echo "$i " 
      #./../ps_functions $i
      mpirun -np $t ./paps_functions $i #./../paps_functions $i
      mpirun -np $t ./paps_functions_precompute_vel $i #./../paps_functions_precompute_vel $i
    done
  done
done

#python3 ../time_iter_plotting.py 

exit 0