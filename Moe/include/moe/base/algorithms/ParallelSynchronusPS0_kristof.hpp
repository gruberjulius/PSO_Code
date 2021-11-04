#pragma once

#include "NumericAlgorithmImpl.hpp"
#include "ParticleSwarm.hpp"
#include "../parameters/PSParameters.hpp"
#include <mpi.h>
#include <vector>
#include <iostream>

//EDIT START
#include <chrono>
//EDIT END
//decommented the file otherwise a redecleration would arise
namespace moe
{

template <typename GenotypeType>
class ParticleSwarmPS_k : public NumericAlgorithm<GenotypeType>
{
    public:
        ParticleSwarmPS_k( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions = 1, std::vector<GenotypeType> _range = { std::numeric_limits<GenotypeType>::lowest() , std::numeric_limits<GenotypeType>::max() });
        ParticleSwarmPS_k( const PSParameters<GenotypeType>& _parameters );

        void run (unsigned int _iterations) override;

    protected:
        void init( unsigned int _iterations ) override;
        //EDIT START
        void init( unsigned int _iterations, int _particle_per_rank);
        std::vector<GenotypeType> getRandomGenotype();
        //EDIT END

    private:
        unsigned int    m_iterations;
        
        float           m_weight,
                        m_coef1,
                        m_coef2;

        std::vector< Moe<GenotypeType> >            m_population;
        std::vector< Moe<GenotypeType> >            m_best_genotypes;
        std::vector< std::vector<GenotypeType> >    m_velocities;
        std::uniform_real_distribution<float>       m_dist_coef;
        //EDIT START
        Moe<GenotypeType> m_best_for_rank;
        std::uniform_real_distribution<> m_dist_genotype;
        std::default_random_engine  m_generator;
        unsigned int m_dimensions;
        //EDIT END
};

template <typename GenotypeType>
ParticleSwarmPS_k<GenotypeType>::ParticleSwarmPS_k( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions, std::vector<GenotypeType> _range )
:NumericAlgorithm<GenotypeType>( _moesPerGen, _dimensions, _range ),
m_weight        ( _weight       ),
m_coef1         ( _coef1        ),
m_coef2         ( _coef2        ),
m_population    ( _moesPerGen   ),
m_best_genotypes(_moesPerGen    ),
m_velocities    (_moesPerGen    ),
m_dist_coef     ( 0.0f, 1.0f    ),
//EDIT START
m_dist_genotype( _range[0], _range[1] ),
m_generator( std::chrono::high_resolution_clock::now().time_since_epoch().count() ),
m_dimensions(_dimensions)
//EDIT END
{
}

template <typename GenotypeType>
ParticleSwarmPS_k<GenotypeType>::ParticleSwarmPS_k( const PSParameters<GenotypeType>& _parameters )
:ParticleSwarmPS_k<GenotypeType>( _parameters.moesPerGen, _parameters.inertia, _parameters.coef1, _parameters.coef2, _parameters.dimensions, _parameters.range )
{
    
}

template <typename GenotypeType>
void ParticleSwarmPS_k<GenotypeType>::init( unsigned int _iterations )
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

template <typename GenotypeType>
std::vector<GenotypeType> ParticleSwarmPS_k<GenotypeType>::getRandomGenotype()
{
    std::vector<GenotypeType> genotype;
    genotype.reserve( m_dimensions );

    while( genotype.size() < m_dimensions )
        genotype.push_back( m_dist_genotype( m_generator ) );

    return genotype;
}

template <typename GenotypeType>
void ParticleSwarmPS_k<GenotypeType>::init( unsigned int _iterations, int _particle_per_rank)
{

    //EDIT
    m_population.resize(_particle_per_rank);
    m_best_genotypes.resize(_particle_per_rank);
    m_velocities.resize(_particle_per_rank);

    //EDIT
    m_iterations = _iterations;

    double min = 1000.0;
    unsigned int    index = 0,
                    count = 0;
    
    for( auto& moe : m_population )
    {
        //EDIT START
        //std::cout << "moe loop\n";
        
        moe.genotype            = getRandomGenotype();

        moe.fitness             = Algorithm<GenotypeType>::m_fitnessFunction( moe );
        //std::cout << "moe init\n";
        m_velocities[count] =  getRandomGenotype();
        //m_velocities[count].push_back(m_dist_genotype( m_generator) );  //          
        //m_velocities[count].push_back(m_dist_genotype( m_generator) );    
        //EDIT END
        
        m_best_genotypes[count].genotype  = moe.genotype;
        m_best_genotypes[count].fitness   = moe.fitness;
        if( moe.fitness < min)
        {
            min = moe.fitness;
            index = count;
        }
        count++;
    }
    //Algorithm<GenotypeType>::m_bestMoe = m_population[index];

    //EDIT
    //temporairly store "global best for current rank"
    //to be gathered in run
    m_best_for_rank = m_population[index];
    //EDIT
}

//PSrun for Parallel Synchronus run
template <typename GenotypeType>
void ParticleSwarmPS_k<GenotypeType>::run( unsigned int _iterations)
{ 

    int rank, size;
    double min_of_iter;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int calc_per_rank;
    calc_per_rank = (int) m_population.size() / size;

    if(rank == 0) calc_per_rank = m_population.size() - (calc_per_rank * (size-1));

    //std::cout<<"Calc on rank: " << rank << " is #: " << calc_per_rank<<"\n";

    this->init( _iterations, calc_per_rank);
    if(rank == 0){
        std::vector<GenotypeType> gen = {m_best_for_rank.genotype[0], m_best_for_rank.genotype[1]};
        Moe<GenotypeType> current_particle = {gen, m_best_for_rank.fitness};
        Algorithm<GenotypeType>::m_bestMoe = current_particle;
    }
    /*for(int i = 0; i < calc_per_rank; ++i){
        std::cout << "Rank: " << rank << ", x: " << m_best_genotypes[i].genotype[0] << ", y: " << m_best_genotypes[i].genotype[1]<< ", fit: " << m_best_genotypes[i].fitness  << std::endl;
    }
    std::cout<<"Init on rank: " << rank << " is done...\n";*/

    double global_best_send_g[3] = {m_best_for_rank.genotype[0], m_best_for_rank.genotype[1], m_best_for_rank.fitness};
    double global_best_recv_g[size*3];
    double global_best_send_b[3] = {m_best_for_rank.genotype[0], m_best_for_rank.genotype[1], m_best_for_rank.fitness}; //b for broadcast, g for gather

    /*std::cout << "Checkpoint 2: managed to initialize the population!" << std::endl;
    std::cout << "size: " << size << std::endl;
    std::cout << "calc_per_rank: " << calc_per_rank << std::endl;
    std::cout << "Population size: " << m_population.size() << std::endl;
    std::cout << "Number of iterations: " << m_iterations << std::endl;*/

    //MPI_Barrier(MPI_COMM_WORLD);
    
    //std::cout << "Checkpoint 3: passed the barrier!" << std::endl;

    for( unsigned int i = 0; i < m_iterations; ++i )
    {
        //MPI_Barrier(MPI_COMM_WORLD);
        //std::cout << "Iteration: " << i << " started on rank: " << rank << std::endl;
        MPI_Gather(global_best_send_g, 3, MPI_DOUBLE, global_best_recv_g, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Gather(global_best_send_g, 3, MPI_DOUBLE, global_best_recv_g, size*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //I SPENT A GOOD FLIPPING 2 HOURS ON FINDING THIS GODDAMN "size*"
        //I'M PISSED
        //REALLY PISSED

        //std::cout<<"Gather on rank: " << rank << " is done...\n";

        /*if(rank ==0){
            for(int i = 0; i<size; ++i){
                std::cout<<global_best_recv_g[(i+2)*3]<<std::endl;
            }
        }*/
        //std::cout << "Checkpoint -1: Gather done!" << std::endl;

        if(rank == 0) 
        {
            unsigned k = 2;
            //std::cout<< "PLEASSSSSSSSSSSSEEEEEEEEEEEEE\n";
            //std::cout<< "global_best_recv_g[2] = " << global_best_recv_g[2]<<"\n";
            min_of_iter = global_best_recv_g[2];
            //std::cout<<"first recieved in gather: " << min_of_iter << "\n";
            for(unsigned j = 0; j < size; ++j) {
                if(global_best_recv_g[j*3+2] < min_of_iter) {
                    min_of_iter = global_best_recv_g[j*3+2];
                    k = j*3+2;
                    //std::cout << "k: " << k << "\n";
                    //std::cout << "min_of_iter: " << min_of_iter << "         m_best_genotypes[j].fitness: " << m_best_genotypes[j].fitness << "           Algorithm<GenotypeType>::m_bestMoe.fitness: " << Algorithm<GenotypeType>::m_bestMoe.fitness << std::endl;
                }
            }
            //std::cout << "m_bestMoe: " << Algorithm<GenotypeType>::m_bestMoe.fitness << ", min_of_iter: " << min_of_iter << std::endl;
            //std::cout << "size of recv: " << global_best_recv_g[k] << std::endl;
            if( min_of_iter < Algorithm<GenotypeType>::m_bestMoe.fitness )
                { 
                    //std::cout << "In if, check global best, k:" << k << "\n";
                    
                    std::vector<GenotypeType> gen = {global_best_recv_g[k-2], global_best_recv_g[k - 1]};
                    Moe<GenotypeType> current_particle = {gen, global_best_recv_g[k]};
                    //std::cout << "one more step\n";
                    //Moe<GenotypeType> tmp_moe;
                    //tmp_moe.genotype[0] = global_best_recv_g[k - 2];
                    //tmp_moe.genotype[1] = global_best_recv_g[k - 1];
                    //tmp_moe.fitness = global_best_recv_g[k];
                    //std::cout << "In if, check global best, finished tmp\n";
                    //std::cout << "min_of_iter: " << min_of_iter << "         m_best_genotypes[j].fitness: " << m_best_genotypes[k].fitness << "           Algorithm<GenotypeType>::m_bestMoe.fitness: " << Algorithm<GenotypeType>::m_bestMoe.fitness << std::endl;
                    //Algorithm<GenotypeType>::m_bestMoe.fitness = min_of_iter; //changed from m_bestMoe to m_bestMoe.fitness
                    
                    //Algorithm<GenotypeType>::m_bestMoe.fitness = global_best_recv_g[k];//tmp_moe;
                    //Algorithm<GenotypeType>::m_bestMoe.genotype[0] = global_best_recv_g[k - 2];
                    //Algorithm<GenotypeType>::m_bestMoe.genotype[1] = global_best_recv_g[k - 1];
                    Algorithm<GenotypeType>::m_bestMoe = current_particle;
                }

            //std::cout << "m_bestMoe: " << Algorithm<GenotypeType>::m_bestMoe.fitness << std::endl;
            //std::cout<< "0 start\n";

            global_best_send_b[0] = Algorithm<GenotypeType>::m_bestMoe.genotype[0];
            //std::cout<< "0 done\n";
            global_best_send_b[1] = Algorithm<GenotypeType>::m_bestMoe.genotype[1];
            global_best_send_b[2] = Algorithm<GenotypeType>::m_bestMoe.fitness;
            //std::cout <<"Is this your problem bro?\n";
        }

        //std::cout << "Checkpoint -1/2: pre Bcast!" << std::endl;

        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(global_best_send_b, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //std::cout << "Checkpoint 0: Bcast done!" << std::endl;

        m_best_for_rank.genotype[0] = global_best_send_b[0];
        m_best_for_rank.genotype[1] = global_best_send_b[1];
        m_best_for_rank.fitness = global_best_send_b[2];
        //std::cout << "Rank: " << rank << ", iter: " << i << ", best val: " << m_best_for_rank.fitness << std::endl;
        //std::cout << "Rank: " << rank << " best moe fit: " << m_best_for_rank.fitness << std::endl;
        
        //std::cout << "Checkpoint 1/2: m_best_for_rank, after Bcast finito!" << std::endl;

        //std::cout<<m_velocities.size();
        for(unsigned int j = 0; j < calc_per_rank; j++) 
        {
            for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )
            {
                //std::cout << "Checkpoint 1: entered double for loop!" << std::endl;
                float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator );
                float   r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

                //velocity update (not to be parallelized)
                m_velocities[j][k] *= m_weight;
                m_velocities[j][k] += m_coef1*r1*( m_best_genotypes[j].genotype[k] - m_population[j].genotype[k] );
                m_velocities[j][k] += m_coef2*r2*( m_best_for_rank.genotype[k] - m_population[j].genotype[k] );//instead of m_best_for_rank.genotype
            
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
            //std::cout << "Checkpoint 2: exited double for loop!" << std::endl;
                
            m_population[j].fitness = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] ); 

            //std::cout << "Checkpoint 3: evaluated funtion!" << std::endl;
            //std::cout << "Rank " << rank << " EVALUATED THE FITNESS FUNCTIONS" << std::endl;

            // Check if new value is local (for that particle) or global best
            if( m_population[j].fitness < m_best_genotypes[j].fitness ) 
            {
                m_best_genotypes[j] = m_population[j];

                if( m_best_genotypes[j].fitness < m_best_for_rank.fitness ) 
                    m_best_for_rank = m_best_genotypes[j];                        
            }
            //std::cout << "Checkpoint 4: finished global update!" << std::endl;

            global_best_send_g[0] = m_best_for_rank.genotype[0];
            global_best_send_g[1] = m_best_for_rank.genotype[1];
            global_best_send_g[2] = m_best_for_rank.fitness;
        } 
            
    }
    

}

} //end namespace Moe

