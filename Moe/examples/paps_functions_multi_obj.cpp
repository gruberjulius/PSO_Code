//STOPPED IMPLEMENTING IT

#include <moe/moe.hpp>
#include <iostream>
#include <chrono>
#include <cmath>
#include <mpi.h>

#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <vector>

#ifdef _MSC_VER
    #include <string>
#endif

#define f(x, y) Ackley_f(x, y)
#define g(x, y) Griewank_f(x, y)
#define pi ((variable_datatype)std::acos(-1))

using variable_datatype = double;

variable_datatype booth(variable_datatype x, variable_datatype y)
{
    variable_datatype  op1 = x + 2*y - 7,
            op2 = 2*x + y - 5;

    return op1*op1 + op2*op2; // 
}

variable_datatype himmelblau(variable_datatype x, variable_datatype y)
{
    variable_datatype  op1 = x*x + y - 11,
            op2 = x + y*y - 7;
    
    return op1*op1 + op2*op2;
}

variable_datatype easom(variable_datatype x, variable_datatype y)
{
    variable_datatype  opcos = -std::cos(x) * std::cos(y),
            xpart = (x-pi)*(x-pi),
            ypart = (y-pi)*(y-pi),
            opexp = std::exp( -( xpart + ypart ) );

    return opcos*opexp;
}

variable_datatype holdertable(variable_datatype x, variable_datatype y)
{
    variable_datatype  optrigo = std::sin(x) * std::cos(y),
            opexp   = std::exp( std::abs( 1 - std::sqrt( x*x + y*y )/pi ) );
    
    return -std::abs( optrigo*opexp );
}

variable_datatype ackley(variable_datatype x, variable_datatype y)
{
    variable_datatype  op1 = -20 * std::exp( 0.2 * std::sqrt( 0.5 * (x*x + y*y) ) ),
            op2 = std::exp( 0.5 * ( std::cos( 2*pi*x ) + std::cos( 2*pi*y ) ) );
    
    return op1 - op2 + std::exp(1)+20;
}

//////////MY_FUNCTIONS/////////////////////////////////////////

double Sphere_f(double x, double y){
    double ret = x*x + y*y;
    return ret;
}

double Schwefel_f(double x, double y){
    double ret = 0;
    double temp = x+y;
    ret = x*x + temp*temp;
    return ret;
}

double Rosenbrock_f(double x, double y){
    double ret = 100. * (y + x*x) + (x - 1)*(x - 1);
    return ret;
}

double Schwefel2_f(double x, double y){
    double ret = x * sin(std::sqrt(std::abs(x)) +  y * sin(std::sqrt(std::abs(y))));
    return -ret;
}

double Rastrigin_f(double x, double y){
    double ret = x*x - 10 * cos(2 * pi * x) + y*y - 10 * cos(2 * pi * y);
    return 10 + ret;
}

double Ackley_f(double x, double y){
    double tmp1 = x*x + y*y;
    double tmp2 = cos(2 * pi * x) + cos(2 * pi * y);
    return 20. - 20. * std::exp(-0.2 * std::sqrt(0.5 * tmp1)) - std::exp(0.5 * tmp2);
}

double Griewank_f(double x, double y){
    double tmp1 = x*x + y*y;
    double tmp2 = cos(x / 1) * cos(y / std::sqrt(2));
    return tmp1 / 4000. - tmp2 + 1;
}

////////////////MY_FUNCTIONS_END///////////////////////////////////


int main(int argc, char * argv[])
{
    long iterations = strtol(argv[1], NULL, 10); //user defined number of iterations
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
        //moe::ParticleSwarmPA_multi_obj<variable_datatype> moether(50, 0.5f, 0.8f, 1.2f, 2, {-10, 10});
        moe::ParticleSwarmPA_multi_obj<variable_datatype> moether( moe::PSParameters<variable_datatype>() 
                                        .withMoesPerGen(100)
                                        .withDimensions(2)
                                        .withRange({-10, 10})
                                        );
        moether.setFitnessFunction( [](auto moe) -> std::vector<variable_datatype>{
            variable_datatype   x = moe.genotype[0],
                    y = moe.genotype[1];
            //variable_datatype   result = f(x, y);
            //variable_datatype result = -1. * (3*x + 2*y);
            std::vector<variable_datatype> result;
            result[0] = f(x, y); //Ackley
            result[1] = g(x, y); //Griewank
            return result;
        });
        //timing
        auto start = std::chrono::high_resolution_clock::now();
        moether.run(iterations); //the run function can be found in include/moe/base/algorithms/ParticleSwarmPA_multi_objOptimization
        auto end = std::chrono::high_resolution_clock::now();

        //output
        std::chrono::duration<variable_datatype> diff = end - start;
        auto best_moe = moether.getBestMoe();
        std::cout   << "Sequential genotype: "     << best_moe.genotype[0] << "\t" << best_moe.genotype[1] << "\n"
            << "Sequential fitness: "      << best_moe.fitness << "\n"
            << "Sequential time spent: "   << diff.count() << " seconds" << std::endl;
    }
*/    
    moe::ParticleSwarmPA_multi_obj<variable_datatype> PAPSOmoether( moe::PSParameters<variable_datatype>() 
                                            .withMoesPerGen(100)
                                            .withDimensions(2)
                                            .withRange({-10, 10})
                                            );
    PAPSOmoether.setFitnessFunction( [](auto moe) -> double
        {
        variable_datatype   x = moe.genotype[0],
                y = moe.genotype[1];
        //variable_datatype   result = f(x, y);
        ///variable_datatype result = -1. * (3*x + 2*y);
        variable_datatype result = Ackley_f(x, y);
        return result;
        });
    
    //Now our PAPSO
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0){ //Master processor
        auto start = std::chrono::high_resolution_clock::now();

        PAPSOmoether.run( iterations, argc, argv );

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<variable_datatype> diff = end - start;

        auto best_moe = PAPSOmoether.getBestMoe();

        std::cout   << "PAPSO genotype: "     << best_moe.genotype[0] << "\t" << best_moe.genotype[1] << "\n"
                << "fitness: "      << best_moe.fitness << "\n"
                << "Master time spent: "   << diff.count() << " seconds" << std::endl;

    //write to file
    
    std::string w_file_name = "results_benchmark/benchmark_times_paps_";
    w_file_name += std::to_string(iterations);
    w_file_name += "_";
    w_file_name += std::to_string(particles);
    w_file_name += ".dat";

    std::ofstream times(w_file_name, std::ios::app); 
    times << s_case << " " << numproces << " " << diff.count() << " " << best_moe.fitness << "\n";
    times.close();

    }
    else {  //Slave processors
        auto start = std::chrono::high_resolution_clock::now();

        PAPSOmoether.run( iterations, argc, argv );

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
