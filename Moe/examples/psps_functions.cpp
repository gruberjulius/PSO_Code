//this is the example function for the PARALLEL SYNCHRONUS FUNCTION 

#include <moe/moe.hpp>
#include <iostream>
#include <chrono>
#include <cmath>

#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include "functions_to_optimize.hpp"


#ifdef _MSC_VER
    #include <string>
#endif

#define OUTPUT_FILE "results/times_pspso" //folder and start of filename
                                          //end always #iter_#threads_#seq/master

#define variable_datatype float

//#define f(x, y) himmelblau(x, y)

int main(int argc, char * argv[])
{

    long iterations = strtol(argv[1], NULL, 10);
    int s_case = atoi(argv[2]); //which function to use
    bool basic_file = atoi(argv[3]);//1;
        //if ture results are written to times_papso.dat
        //if not results written to times_papso_#iter_#threads.dat
    int particles = atoi(argv[4]);

    //Initializing MPI
    MPI_Init(&argc,&argv);
    int rank, numproces;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);

    
    std::cout << "Pre func decl\n";
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
            //std::cout << "Only 5 functions used (0,...,4)" << std::endl;
            break;
    }

    std::cout << range[0] << " " << range[1] << std::endl;
    moe::ParticleSwarmPS<variable_datatype> PSPSOmoether( moe::PSParameters<variable_datatype>() 
                                                    .withMoesPerGen(particles)
                                                    .withDimensions(2)
                                                    .withRange(range)
                                                    );

    PSPSOmoether.setFitnessFunction( [s_case](auto moe) -> double 
    {
        variable_datatype   x = moe.genotype[0],
                            y = moe.genotype[1];
        double result = 0;
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
                //std::cout << "Only 5 functions used (0,...,4)" << std::endl;
                break;
        }

        return result;
    });

/* ///OLD INITIALIZATION START #########################################
    moe::ParticleSwarmPS<variable_datatype> PSPSOmoether( moe::PSParameters<variable_datatype>() 
                                                    .withMoesPerGen(100)
                                                    .withDimensions(2)
                                                    .withRange({-10, 10})
                                                    );
    PSPSOmoether.setFitnessFunction( [s_case](auto moe) -> double 
    {
        variable_datatype   x = moe.genotype[0],
                            y = moe.genotype[1];
        variable_datatype result = 0;
        //variable_datatype result = 3*x + 2*y;

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
                std::cout << "Only 5 functions used (0,...,4)" << std::endl;
                break;
        }

        return result;
    });
    */ //OLD INITIALIZATION END ###################################

    //Now our PSPSO
    MPI_Barrier(MPI_COMM_WORLD);
    
    auto start = std::chrono::high_resolution_clock::now();

    //int calc_per_rank;
    //calc_per_rank = m_population.size() / size;

    //static variable_datatype max_of_iter;

    //!!!!!!!!!!!create vector to pass into run and but local bests in there!!!!!!!!!!!!
    PSPSOmoether.run( iterations);

    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> diff = end - start;

    if(rank == 0) {
        
        auto best_moe_PSPSO = PSPSOmoether.getBestMoe();

        std::cout   << "PSPSO: \n";
        std::cout   << "PSPSO genotype: "       << best_moe_PSPSO.genotype[0] << "\t" << best_moe_PSPSO.genotype[1] << "\n"
                    << "fitness: "              << best_moe_PSPSO.fitness << "\n"
                    << "time spent: "    << diff.count() << " seconds\n\n" << std::endl; 

    std::string w_file_name = OUTPUT_FILE;
    if(!basic_file){
        std::cout << "Results written to long filename!" << std::endl;
        w_file_name += "_";
        w_file_name += std::to_string(iterations);
        w_file_name += "_";
        w_file_name += std::to_string(particles);
        //w_file_name += "_m";
    }
    w_file_name += ".dat";

    std::ofstream times(w_file_name, std::ios::app); 
    times << s_case << " " << numproces << " " << diff.count() << " " << best_moe_PSPSO.fitness << "\n";
    times.close();
    }
    
    //write times to file
    /*std::ofstream times("results/times_psps.txt", std::ios::app);
    times << iterations << ';' <<diff.count() << "\n";
    times.close(); */
/*kc commented it out and added into rank == 0
    std::ofstream times("times_papso.dat", std::ios::app); 
    times << diff.count() << " " << PSPSOmoether.getBestMoe().fitness << "\n";
    times.close();
*/
    //std::cout << "Rank " << rank << " got to this point." << std::endl;       
    
    MPI_Finalize();

    return 0;
}
