#only bash to call in the end

#first let PGA work
#need setup_and_run.sh  |  wscript_good  |  ea_edited.hpp  |  func1.cpp (*number of functions we test) in directory PGA
cd PGA
bash setup_and_run.sh
#timing for each function in resultPGA${N}.dat file with [timing bestfit]

echo "Finished with PGA"
cd ../MOE

#need iterator.sh | time_iter_plotting | MOE stuff in MOE
bash setup_and_run_pso.sh

# now HERE plotting and other use of data



#More or less how it should look like (NOT FINISHED)
# global_folder
# |--PGA
#    |--sferes2
#       |**results.dat #file contains [IndexOfTheFunction NumberOfThreads Time BestFitnessValueAttained] * #calls of the same executable
# |--MOE
#    |--build
#    	|**results_papso.dat #file contains [IndexOfTheFunction NumberOfThreads Time BestFitnessValueAttained] * #calls of the same executable
#    	|**results_pspso.dat #idem
#    	|**results_pso.dat #idem
