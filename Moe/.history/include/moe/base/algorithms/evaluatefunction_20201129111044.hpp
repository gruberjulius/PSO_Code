#include"mpi.h"
//input: the different values that should be evaluated in parallel
//post: a vector all fitness values for every process 

//argument: functions von den neu evaluated functions
void evaluatefunction(std::vector<double> function_results){
    
    for(unsigned int j = 0; j < function_results.size(); j++){
        for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )
        {
            float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator );
            float   r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

            //velocity update (not to be parallelized)
            m_velocities[j][k] *= m_weight;
            m_velocities[j][k] += m_coef1*r1*( m_best_genotypes[j].genotype[k] - m_population[j].genotype[k] );
            m_velocities[j][k] += m_coef2*r2*( Algorithm<GenotypeType>::m_bestMoe.genotype[k] - m_population[j].genotype[k] );
        
            //position update (not to be parallelized)
            m_population[j].genotype[k] += m_velocities[j][k];
        
            //checks if genotype actual dimension is still in provided search space
            if( m_population[j].genotype[k] != std::max( std::min( m_population[j].genotype[k], NumericAlgorithm<GenotypeType>::m_range[1]), NumericAlgorithm<GenotypeType>::m_range[0] ) )
            {
                m_population[j].genotype = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                m_velocities[j] = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                break;
            }

            //This is the local best update (not to be parallelized)
            if( function_results[j] > m_best_genotypes[j].fitness ) 
            {
                m_population[j].fitness = function_results[j];
                m_best_genotypes[j] = m_population[j];
                //This is the global best update (not to be parallelized)
                if( m_best_genotypes[j].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness ) 
                    Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[j];
            }
        }
    }
    /*

    std::vector<MPI_DOUBLE> functionresults(m_population.size());


    for(unsigned int j = 0; j < m_population.size(); j++){
        m_population[j].fitness = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] ); 
    }
    */
    //probably add the mpi stuff here
    /*
    MPI_Init(&argc, &argv);

    for(unsigned int j = 0; j < m_population.size(); j++){
        functionresults[j] = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] ); 
    }

    MPI_Finalize();
    
    return functionresults;
    */
}

//velocity und position update implementieren, 
//jeweiligen, local best position updaten
//
