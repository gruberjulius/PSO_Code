#!/bin/bash
#ab jz nur fp hier
#
#  std::ofstream times("times_paps.dat", std::ios::app); 
#   times << iterations << ';' <<diff.count() << "\n";
#    times.close();



FOLDER=results_benchmark
echo "$FOLDER"
rm -r $FOLDER
mkdir $FOLDER
cd $FOLDER

for t in 2 3 7 8 15 16
do
    for i in 1000 10000 100000
    do
        #if [[ -f results_benchmark/benchmark_times_paps_$i_$t.dat ]]
        #then
        #    rm results_benchmark/benchmark_times_paps_$i_$t.dat
        #fi

        #cat > results_benchmark/benchmark_times_paps_$i_$t.dat
        for x in [0..100..1]
        do
        mpirun -np $t ../paps_functions $i
        done
    done
done

exit 0