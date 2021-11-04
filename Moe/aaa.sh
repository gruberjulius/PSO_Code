i=1000
for t in 2 3 4  #2 3 4 5 7 8 9 15 16 17 31 32 33 #{2..16..1}  #number of threads
do
    for f in in 3 #{0..5..1} # funciton id
    do
        for x in {1..10}
        do
            #mpirun -np $t paps_functions_send_move_termcrit $i $f
            #mpirun -np $t paps_functions_send_move $i $f
        done
    done
done
