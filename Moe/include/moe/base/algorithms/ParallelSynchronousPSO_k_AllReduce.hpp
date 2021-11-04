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
class ParticleSwarmPS_k_aR : public NumericAlgorithm<GenotypeType>
{
    public:
        ParticleSwarmPS_k_aR( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions = 1, std::vector<GenotypeType> _range = { std::numeric_limits<GenotypeType>::lowest() , std::numeric_limits<GenotypeType>::max() });
        ParticleSwarmPS_k_aR( const PSParameters<GenotypeType>& _parameters );

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
        //Moe<GenotypeType> m_best_for_rank;
        std::uniform_real_distribution<> m_dist_genotype;
        std::default_random_engine  m_generator;
        unsigned int m_dimensions;
        //EDIT END
};

template <typename GenotypeType>
ParticleSwarmPS_k_aR<GenotypeType>::ParticleSwarmPS_k_aR( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions, std::vector<GenotypeType> _range )
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
ParticleSwarmPS_k_aR<GenotypeType>::ParticleSwarmPS_k_aR( const PSParameters<GenotypeType>& _parameters )
:ParticleSwarmPS_k_aR<GenotypeType>( _parameters.moesPerGen, _parameters.inertia, _parameters.coef1, _parameters.coef2, _parameters.dimensions, _parameters.range )
{
    
}

template <typename GenotypeType>
void ParticleSwarmPS_k_aR<GenotypeType>::init( unsigned int _iterations )
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
std::vector<GenotypeType> ParticleSwarmPS_k_aR<GenotypeType>::getRandomGenotype()
{
    std::vector<GenotypeType> genotype;
    genotype.reserve( m_dimensions );

    while( genotype.size() < m_dimensions )
        genotype.push_back( m_dist_genotype( m_generator ) );

    return genotype;
}

template <typename GenotypeType>
void ParticleSwarmPS_k_aR<GenotypeType>::init( unsigned int _iterations, int _particle_per_rank)
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
        if( min > moe.fitness )
        {
            min = moe.fitness;
            index = count;
        }
        count++;
    }
    Algorithm<GenotypeType>::m_bestMoe = m_population[index];

    //EDIT
    //temporairly store "global best for current rank"
    //to be gathered in run
    //m_best_for_rank = m_population[index];
    //EDIT
}
/*
struct particlereduce{
        double genotype1;
        double genotype2;
        double fitness;
        };
//void min_fit(particlereduce*, particlereduce*, int*, MPI_Datatype*);
void min_fit(particlereduce* invec, particlereduce* inoutvec, int* len, MPI_Datatype *datatype){
    for(int i = 0; i < *len; ++i){
        if(invec[i].fitness < inoutvec[i].fitness){
            inoutvec[i] = invec[i];
        }
    }
}*/
//PSrun for Parallel Synchronus run
template <typename GenotypeType>
void ParticleSwarmPS_k_aR<GenotypeType>::run( unsigned int _iterations)
{ 

    int rank, size;
    double max_of_iter;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    

    int blocklen[2]={2,1};
    MPI_Aint displ[2]={0,16};
    MPI_Datatype types[2]={MPI_DOUBLE,MPI_DOUBLE};
    MPI_Datatype MinMoeType;
    MPI_Type_create_struct(
        2,
        blocklen,
        displ,
        types,
        &MinMoeType
    );
    MPI_Type_commit(&MinMoeType);


    MPI_Op min_moe_op;
    MPI_Op_create( (MPI_User_function *)min_fit,0,&min_moe_op);

    
    int calc_per_rank;
    calc_per_rank = (int) m_population.size() / size;

    if(rank == 0) calc_per_rank = m_population.size() - (calc_per_rank * (size-1));

    //std::cout<<"Calc on rank: " << rank << " is #: " << calc_per_rank<<"\n";

    this->init( _iterations, calc_per_rank);
    if(rank == 0){
        std::vector<GenotypeType> gen = {Algorithm<GenotypeType>::m_bestMoe.genotype[0], Algorithm<GenotypeType>::m_bestMoe.genotype[1]};
        Moe<GenotypeType> current_particle = {gen, Algorithm<GenotypeType>::m_bestMoe.fitness};
        Algorithm<GenotypeType>::m_bestMoe = current_particle;
    }
    /*for(int i = 0; i < calc_per_rank; ++i){
        std::cout << "Rank: " << rank << ", x: " << m_best_genotypes[i].genotype[0] << ", y: " << m_best_genotypes[i].genotype[1]<< ", fit: " << m_best_genotypes[i].fitness  << std::endl;
    }
    std::cout<<"Init on rank: " << rank << " is done...\n";*/

    particlereduce global_best_send_g = {Algorithm<GenotypeType>::m_bestMoe.genotype[0], Algorithm<GenotypeType>::m_bestMoe.genotype[1], Algorithm<GenotypeType>::m_bestMoe.fitness};
    particlereduce global_best_recv_g[size];

    /*std::cout << "Checkpoint 2: managed to initialize the population!" << std::endl;
    std::cout << "size: " << size << std::endl;
    std::cout << "calc_per_rank: " << calc_per_rank << std::endl;
    std::cout << "Population size: " << m_population.size() << std::endl;
    std::cout << "Number of iterations: " << m_iterations << std::endl;*/

    //MPI_Barrier(MPI_COMM_WORLD);
    
    //std::cout << "Checkpoint 3: passed the barrier!" << std::endl;

    for( unsigned int i = 0; i < m_iterations; ++i )
    {
        
        MPI_Allreduce(&global_best_send_g, &global_best_recv_g, 1, MinMoeType, min_moe_op, MPI_COMM_WORLD);


        //std::cout << "Checkpoint 0: Bcast done!" << std::endl;

        Algorithm<GenotypeType>::m_bestMoe.genotype[0] = global_best_send_g.genotype1;
        Algorithm<GenotypeType>::m_bestMoe.genotype[1] = global_best_send_g.genotype2;
        Algorithm<GenotypeType>::m_bestMoe.fitness = global_best_send_g.fitness;
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
                m_velocities[j][k] += m_coef2*r2*( Algorithm<GenotypeType>::m_bestMoe.genotype[k] - m_population[j].genotype[k] );//instead of m_best_for_rank.genotype
            
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

                if( m_best_genotypes[j].fitness < Algorithm<GenotypeType>::m_bestMoe.fitness ) 
                    Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[j];                        
            }
            //std::cout << "Checkpoint 4: finished global update!" << std::endl;

            global_best_send_g.genotype1 = Algorithm<GenotypeType>::m_bestMoe.genotype[0];
            global_best_send_g.genotype2 = Algorithm<GenotypeType>::m_bestMoe.genotype[1];
            global_best_send_g.fitness = Algorithm<GenotypeType>::m_bestMoe.fitness;
        } 
            
    }
    
    //Algorithm<GenotypeType>::m_bestMoe = current_particle;

    MPI_Op_free(&min_moe_op);


}

} //end namespace Moe

