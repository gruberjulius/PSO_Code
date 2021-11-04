#include"mpi.h"
//was muss übergeben werden
//input: the different values that should be evaluated in parallel
//post: a vector all fitness values for every process 

//argument: functions von den neu evaluated functions
void evaluatefunction(m_population);
    std::vector<MPI_DOUBLE> functionresults(m_population.size());

    /*
    for(unsigned int j = 0; j < m_population.size(); j++){
        m_population[j].fitness = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] ); 
    }
    */
    //probably add the mpi stuff here
    MPI_Init(&argc, &argv);

    for(unsigned int j = 0; j < m_population.size(); j++){
        functionresults[j] = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] ); 
    }

    MPI_Finalize();
    return functionresults;
}

//velocity und position update implementieren, 
//jeweiligen, local best position updaten
//
