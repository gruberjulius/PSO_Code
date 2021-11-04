#pragma once

#include "NumericAlgorithmImpl.hpp"
#include "ParticleSwarm.hpp"
#include "../parameters/PSParameters.hpp"
#include <mpi.h>
#include <vector>
//decommented the file otherwise a redecleration would arise
namespace moe
{

template <typename GenotypeType>
class ParticleSwarmPS : public NumericAlgorithm<GenotypeType>
{
    public:
        ParticleSwarmPS( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions = 1, std::vector<GenotypeType> _range = { std::numeric_limits<GenotypeType>::lowest() , std::numeric_limits<GenotypeType>::max() });
        ParticleSwarmPS( const PSParameters<GenotypeType>& _parameters );

        void run( unsigned int _iterations ) override;

    protected:
        void init( unsigned int _iterations ) override;

    private:
        unsigned int    m_iterations;
        
        float           m_weight,
                        m_coef1,
                        m_coef2;

        std::vector< Moe<GenotypeType> >            m_population;
        std::vector< Moe<GenotypeType> >            m_best_genotypes;
        std::vector< std::vector<GenotypeType> >    m_velocities;
        std::uniform_real_distribution<float>       m_dist_coef;
};

template <typename GenotypeType>
ParticleSwarmPS<GenotypeType>::ParticleSwarmPS( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions, std::vector<GenotypeType> _range )
:NumericAlgorithm<GenotypeType>( _moesPerGen, _dimensions, _range ),
m_weight        ( _weight       ),
m_coef1         ( _coef1        ),
m_coef2         ( _coef2        ),
m_population    ( _moesPerGen   ),
m_best_genotypes(_moesPerGen    ),
m_velocities    (_moesPerGen    ),
m_dist_coef     ( 0.0f, 1.0f    )
{
}

template <typename GenotypeType>
ParticleSwarmPS<GenotypeType>::ParticleSwarmPS( const PSParameters<GenotypeType>& _parameters )
:ParticleSwarmPS<GenotypeType>( _parameters.moesPerGen, _parameters.inertia, _parameters.coef1, _parameters.coef2, _parameters.dimensions, _parameters.range )
{
    
}

template <typename GenotypeType>
void ParticleSwarmPS<GenotypeType>::init( unsigned int _iterations )
{
    m_iterations = _iterations;

    double max = 0.0;
    unsigned int    index = 0,
                    count = 0;
    
    for( auto& moe : m_population )
    {
        moe.genotype            = NumericAlgorithm<GenotypeType>::getRandomGenotype();
        moe.fitness             = Algorithm<GenotypeType>::m_fitnessFunction( moe );
        
        m_velocities[count]       = NumericAlgorithm<GenotypeType>::getRandomGenotype();
        
        m_best_genotypes[count].genotype  = moe.genotype;
        m_best_genotypes[count].fitness   = moe.fitness;
        if( max < moe.fitness )
        {
            max = moe.fitness;
            index = count;
        }
        count++;
    }

    Algorithm<GenotypeType>::m_bestMoe = m_population[index];
}

//PSrun for Parallel Synchronus run
template <typename GenotypeType>
void ParticleSwarmPS<GenotypeType>::run( unsigned int _iterations )
{   
    
    this->init( _iterations );

    for( unsigned int i = 0; i < m_iterations; i++ )
    {
        /*for(unsigned int j = 0; j <  m_population.size(); j++)
        {
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
        }*/
        
        int rank, size;
        int calc_per_rank;
        double start, finish, loc_elapsed, elapsed, max_of_iter;
        std::vector<double> fitness_values(m_population.size());
         
        //MPI_init(&argc,&argv)
        //MPI_Init(int *erfc, char ***argv); //got error with other  version, now it is fixed
        MPI_Init(&argc,&argv);

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        calc_per_rank = m_population.size() / size;

        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        
        if(rank != size - 1) {
            for(int j = rank * calc_per_rank; j < (rank + 1) * calc_per_rank; ++j) 
            {
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
                }
                    
                m_population[j].fitness = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] ); 

                if( m_population[j].fitness > m_best_genotypes[j].fitness ) 
                {
                    m_best_genotypes[j] = m_population[j];

                    MPI_Reduce(&m_best_genotypes[j], &max_of_iter, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                    //This is the global best update (not to be parallelized)
                    /*if(rank == 0) 
                    {
                        if( m_best_genotypes[j].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness ) 
                            Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[j];
                    } */
                        
                }
            } 
        }
        else {
            for(int j = rank * calc_per_rank; j < m_population.size(); ++j) {
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
                }
                    
                m_population[j].fitness = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] ); 

                if( m_population[j].fitness > m_best_genotypes[j].fitness ) 
                {
                    m_best_genotypes[j] = m_population[j];

                    MPI_Reduce(&m_best_genotypes[j], &max_of_iter, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                    //This is the global best update (not to be parallelized)
                    /*if(rank == 0) 
                    {
                        if( m_best_genotypes[j].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness ) 
                            Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[j];
                    } */
                        
                }
            }
        }



        finish = MPI_Wtime();
        loc_elapsed = finish - start;

        MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); //!!!!!!does reduce get rid of the local values???!!!!!!

        MPI_Finalize(); //declaring a res function like here https://riptutorial.com/de/mpi/example/8638/init---finalisieren

        if( max_of_iter > Algorithm<GenotypeType>::m_bestMoe.fitness ){ 
                            Algorithm<GenotypeType>::m_bestMoe.fitness = max_of_iter; //changed from m_bestMoe to m_bestMoe.fitness
        }
        //std::cout << "3rds slot of function_values verctor: " << fitness_values[2] << std::endl;

        //evaluatefunction(fitness_values);

        /*for( unsigned int j = 0; j < m_population.size(); j++ )
        {
            for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )
            {
                //set the random coefficients for the volocity update (r1 for cognitive and r2 for social)
                float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator ),
                        r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

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
            }

            //This is the f(x), the part we need to parallelize
            m_population[j].fitness = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] ); 
            
            //This is the local best update (not to be parallelized)
            if( m_population[j].fitness > m_best_genotypes[j].fitness ) 
            {
                m_best_genotypes[j] = m_population[j];
                //This is the global best update (not to be parallelized)
                if( m_best_genotypes[j].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness ) 
                    Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[j];
            }
        }*/
    }
    

}

} //end namespace Moe


//take out line 118 (fitness evaluation) and parallelize this process. 
//Insert values into a vector of fitness values.
//create a function with rest of population for loop, that takes this vector as an input and updates postion, velocity, local best and global best.

//evaluate(std::vector<double> fitness_values);