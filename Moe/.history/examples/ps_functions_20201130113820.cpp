#include <moe/moe.hpp>
#include <iostream>
#include <chrono>
#include <cmath>
#include <mpi.h>

#ifdef _MSC_VER
    #include <string>
#endif

#define f(x, y) himmelblau(x, y)
#define pi ((ex_float)std::acos(-1))

using ex_float = double;

ex_float booth(ex_float x, ex_float y)
{
    ex_float  op1 = x + 2*y - 7,
            op2 = 2*x + y - 5;

    return op1*op1 + op2*op2; // 
}

ex_float himmelblau(ex_float x, ex_float y)
{
    ex_float  op1 = x*x + y - 11,
            op2 = x + y*y - 7;
    
    return op1*op1 + op2*op2;
}

ex_float easom(ex_float x, ex_float y)
{
    ex_float  opcos = -std::cos(x) * std::cos(y),
            xpart = (x-pi)*(x-pi),
            ypart = (y-pi)*(y-pi),
            opexp = std::exp( -( xpart + ypart ) );

    return opcos*opexp;
}

ex_float holdertable(ex_float x, ex_float y)
{
    ex_float  optrigo = std::sin(x) * std::cos(y),
            opexp   = std::exp( std::abs( 1 - std::sqrt( x*x + y*y )/pi ) );
    
    return -std::abs( optrigo*opexp );
}

ex_float ackley(ex_float x, ex_float y)
{
    ex_float  op1 = -20 * std::exp( 0.2 * std::sqrt( 0.5 * (x*x + y*y) ) ),
            op2 = std::exp( 0.5 * ( std::cos( 2*pi*x ) + std::cos( 2*pi*y ) ) );
    
    return op1 - op2 + std::exp(1)+20;
}

int main(int argc, char * argv[])
{
    unsigned int iterations = 1000; //user defined number of iterations

    //Initializing MPI
    MPI_Init(&argc,&argv);
    int rank, numproces;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);

    //Run the sequential SPSO on the master processor first
    if(rank == 0){
        //init
        moe::ParticleSwarm<ex_float> moether(50, 0.5f, 0.8f, 1.2f, 2, {-10, 10});
        moether.setFitnessFunction( [](auto moe) -> double{
            ex_float   x = moe.genotype[0],
                    y = moe.genotype[1];
            //ex_float   result = f(x, y);
            ex_float result = 3*x + 2*y;
            return -result;
        });
        //timing
        auto start = std::chrono::high_resolution_clock::now();
        moether.run(iterations); //the run function can be found in include/moe/base/algorithms/particleSwarmOptimization
        auto end = std::chrono::high_resolution_clock::now();

        //output
        std::chrono::duration<ex_float> diff = end - start;
        auto best_moe = moether.getBestMoe();
        std::cout   << "Sequential genotype: "     << best_moe.genotype[0] << "\t" << best_moe.genotype[1] << "\n"
            << "Sequential fitness: "      << best_moe.fitness << "\n"
            << "Sequential time spent: "   << diff.count() << " seconds" << std::endl;
    }
    
    moe::ParticleSwarmPA<ex_float> PAPSOmoether( moe::PSParameters<ex_float>() 
                                            .withMoesPerGen(100)
                                            .withDimensions(2)
                                            .withRange({-10, 10})
                                            );
    PAPSOmoether.setFitnessFunction( [](auto moe) -> double{
        ex_float   x = moe.genotype[0],
                y = moe.genotype[1];
        //ex_float   result = f(x, y);
        ex_float result = 3*x + 2*y;
        return -result;
    });
//Now our PAPSO
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){ //Master processor
        auto start = std::chrono::high_resolution_clock::now();

        PAPSOmoether.run( iterations, argc, argv );

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<ex_float> diff = end - start;

        auto best_moe = PAPSOmoether.getBestMoe();

        std::cout   << "PAPSO genotype: "     << best_moe.genotype[0] << "\t" << best_moe.genotype[1] << "\n"
                << "fitness: "      << best_moe.fitness << "\n"
                << "Master time spent: "   << diff.count() << " seconds" << std::endl;
    }
    else {  //Slave processors
        auto start = std::chrono::high_resolution_clock::now();

        PAPSOmoether.run( iterations, argc, argv );

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<ex_float> diff = end - start;
 
        std::cout   << "Slave " << rank << ", time spent: "   << diff.count() << " seconds" << std::endl;
    }

    //End of parallel part and program
    MPI_Finalize();
    return 0; //would add this for style reasons
}
