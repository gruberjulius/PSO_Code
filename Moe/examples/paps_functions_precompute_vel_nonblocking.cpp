#include <moe/moe.hpp>
#include <iostream>
#include <chrono>
#include <cmath>
#include <mpi.h>

#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include<unistd.h>

#include "functions_to_optimize.hpp"

#ifdef _MSC_VER
    #include <string>
#endif

using variable_datatype = double;

#define f(x, y) himmelblau(x, y)
#define OUTPUT_FILE "results/times_papso_pre_nonblc" //folder and start of filename
                                          //end always #iter_#threads_#seq/master



int main(int argc, char * argv[])
{

    long iterations = strtol(argv[1], NULL, 10);
    int s_case = atoi(argv[2]); //which function to use
    bool basic_file = atoi(argv[3]);//1;
        //if ture results are written to times_papso.dat
        //if not results written to times_papso_#iter_#threads.dat
    int particles = atoi(argv[4]);

    //std::cin>>iterations;

    //Initializing MPI
    MPI_Init(&argc,&argv);
    int rank, numproces;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);
/*
    //Run the sequential SPSO on the master processor first
    if(rank == 0){
        //init
        //moe::ParticleSwarm<variable_datatype> moether(50, 0.5f, 0.8f, 1.2f, 2, {-10, 10});
        moe::ParticleSwarm<variable_datatype> moether( moe::PSParameters<variable_datatype>()
                                        .withMoesPerGen(100)
                                        .withDimensions(2)
                                        .withRange({-10, 10})
                                        );
        moether.setFitnessFunction( [s_case](auto moe) -> double {
            variable_datatype   x = moe.genotype[0],
                    y = moe.genotype[1];
            //variable_datatype   result = f(x, y);
        variable_datatype result = 0;
        //cvariable_datatype result = 3*x + 2*y;

        switch(s_case) {
            case 0:
                result = Griewank_f(x,y);
                break;
            case 1:
                result = Ackley_f(x,y);
                break;
            case 2:
                result = holdertable(x,y);
                break;
            case 3:
                result = Rastrigin_f(x,y);
                break;
            case 4:
                result = Schwefel2_f(x,y);
                break;
            case 5:
                result = VERY_Complicated(x,y);
                break;
            default:
                std::cout << "Only 6 functions used (0,...,5)" << std::endl;
                break;
        }
            return result;
        });
        //timing
        auto start = std::chrono::high_resolution_clock::now();
        moether.run(iterations); //the run function can be found in include/moe/base/algorithms/particleSwarmOptimization
        auto end = std::chrono::high_resolution_clock::now();

        //output
        std::chrono::duration<variable_datatype> diff = end - start;
        auto best_moe = moether.getBestMoe();
        std::cout   << "Sequential genotype: "     << best_moe.genotype[0] << "\t" << best_moe.genotype[1] << "\n"
            << "Sequential fitness: "      << best_moe.fitness << "\n"
            << "Sequential time spent: "   << diff.count() << " seconds" << std::endl;

    std::string w_file_name = OUTPUT_FILE;
    if(!basic_file){
        std::cout<<"Results written to long filename!" << std::endl;
        w_file_name += "_";
        w_file_name += std::to_string(iterations);
        w_file_name += "_";
        w_file_name += std::to_string(numproces);
        w_file_name += "_s";

        w_file_name += ".dat";

        std::ofstream times(w_file_name, std::ios::app);
        times << s_case << " " << numproces << " " << diff.count() << " " << best_moe.fitness << "\n";
        times.close();
    }
    }
*/
    std::vector<variable_datatype> range;

    switch(s_case) { // switch case to set correct range of function
        case 0: 
            range = {-600, 600};
            break;
        case 1: 
            range = {-32, 32};
            break;
        case 2: 
            range = {-10, 10};
            break;
        case 3: 
            range = {-5.12, 5.12};
            break;
        case 4: 
            range = {-100, 100};
            break;
        case 5:
            range = {-100, 100};
            break;
        default:
            std::cout << "Only 5 functions used (0,...,4)" << std::endl;
            break;
    }
    moe::ParticleSwarmPA_prec_vel_nonblc<variable_datatype> PAPSOmoether( moe::PSParameters<variable_datatype>()
                                            .withMoesPerGen(particles)
                                            .withDimensions(2)
                                            .withRange(range)
                                            );
    PAPSOmoether.setFitnessFunction( [s_case](auto moe) -> double
        {
        variable_datatype   x = moe.genotype[0],
                y = moe.genotype[1];
        //variable_datatype   result = f(x, y);
        variable_datatype result = 0;
        //cvariable_datatype result = 3*x + 2*y;

        switch(s_case) {
            case 0:
                result = Griewank_f(x,y);
                break;
            case 1:
                result = Ackley_f(x,y);
                break;
            case 2:
                result = holdertable(x,y);
                break;
            case 3:
                result = Rastrigin_f(x,y);
                break;
            case 4:
                result = Schwefel2_f(x,y);
                break;
            case 5:
                result = VERY_Complicated(x,y);
                break;
            default:
                std::cout << "Only 6 functions used (0,...,5)" << std::endl;
                break;
        }
        return result;
        });

    //Now our PAPSO
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0){ //Master processor
        auto start = std::chrono::high_resolution_clock::now();

        PAPSOmoether.run( iterations);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<variable_datatype> diff = end - start;

        auto best_moe = PAPSOmoether.getBestMoe();

        std::cout   << "PAPSO genotype: "     << best_moe.genotype[0] << "\t" << best_moe.genotype[1] << "\n"
                << "fitness: "      << best_moe.fitness << "\n"
                << "Master time spent: "   << diff.count() << " seconds" << std::endl;

    //write to file

    std::string w_file_name = OUTPUT_FILE;
    if(!basic_file){
        std::cout<<"Results written to long filename!" << std::endl;
        w_file_name += "_";
        w_file_name += std::to_string(iterations);
        w_file_name += "_";
        w_file_name += std::to_string(particles);
        //w_file_name += "_m";
    }
    w_file_name += ".dat";

    std::ofstream times(w_file_name, std::ios::app);
    times << s_case << " " << numproces << " " << diff.count() << " " << best_moe.fitness << "\n";
    times.close();

    }
    else {  //Slave processors
        auto start = std::chrono::high_resolution_clock::now();

        PAPSOmoether.run( iterations);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<variable_datatype> diff = end - start;

        //std::cout   << "Slave " << rank << ", time spent: "   << diff.count() << " seconds" << std::endl;
    }

    //End of parallel part and program
    //std::cout<<rank<<std::endl;
    int error = MPI_Finalize();
    std::cout<<"Post finalize with error:"<< error << ", with rank " << rank <<std::endl;
    return 0; //would add this for style reasons
}
