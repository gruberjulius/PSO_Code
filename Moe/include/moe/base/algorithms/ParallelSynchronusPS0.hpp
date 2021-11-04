#pragma once

#include "NumericAlgorithmImpl.hpp"
#include "ParticleSwarm.hpp"
#include "../parameters/PSParameters.hpp"
#include <mpi.h>
#include <vector>
#include <iostream>
#include <random>
//decommented the file otherwise a redecleration would arise
namespace moe
{

template <typename GenotypeType>
class ParticleSwarmPS : public NumericAlgorithm<GenotypeType>
{
    public:
        ParticleSwarmPS( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions = 1, std::vector<GenotypeType> _range = { std::numeric_limits<GenotypeType>::lowest() , std::numeric_limits<GenotypeType>::max() });
        ParticleSwarmPS( const PSParameters<GenotypeType>& _parameters );

        void run (unsigned int _iterations) override;

    protected:
        void init( unsigned int _iterations ) override;
        //void init( unsigned int _iterations, unsigned int _partialSize );

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

/*template <typename GenotypeType>
void ParticleSwarmPS<GenotypeType>::init( unsigned int _iterations, unsigned int _partialSize )
{
    m_iterations = _iterations;

    double max = 0.0;
    unsigned int    index = 0;
                    //count = 0;


    
    for(int i = 0; i < _partialSize; ++i)
    {
        m_population[i].genotype            = NumericAlgorithm<GenotypeType>::getRandomGenotype();
        m_population[i].fitness             = Algorithm<GenotypeType>::m_fitnessFunction( m_population[i] );

        //std::cout << "PeePee and this:    " << _partialSize << "and " << m_population[i].genotype[1] << std::endl;
        
        m_velocities[i]       = NumericAlgorithm<GenotypeType>::getRandomGenotype();
        
        m_best_genotypes[i].genotype  = m_population[i].genotype;
        m_best_genotypes[i].fitness   = m_population[i].fitness;
        if( max < m_population[i].fitness )
        {
            max = m_population[i].fitness;
            index = i;
        }
        //count++;
    }
    for(int i = _partialSize; i < m_population.size(); ++i)
    {
        m_population[i].genotype            = NumericAlgorithm<GenotypeType>::getRandomGenotype();
        m_population[i].fitness             = Algorithm<GenotypeType>::m_fitnessFunction( m_population[i] );

        //std::cout << "PeePee and this:    " << _partialSize << "and " << m_population[i].genotype[1] << std::endl;
        
        m_velocities[i]       = NumericAlgorithm<GenotypeType>::getRandomGenotype();
        
        m_best_genotypes[i].genotype  = m_population[i].genotype;
        m_best_genotypes[i].fitness   = m_population[i].fitness;
        if( max < m_population[i].fitness )
        {
            max = m_population[i].fitness;
            index = i;
        }
    }

    Algorithm<GenotypeType>::m_bestMoe = m_population[index];
}
*/
//PSrun for Parallel Synchronus run
template <typename GenotypeType>
void ParticleSwarmPS<GenotypeType>::run( unsigned int _iterations)
{ 

    int rank, size;
    double max_of_iter;
    int pop_size, rest;
    int calc_per_rank;
    int iter;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank == 0) {
        this->init( _iterations );
        pop_size = m_population.size();
        calc_per_rank = (int) pop_size / size;
        //std::cout << "iterations: " << _iterations << std::endl;
        //rest = pop_size - (size - 1) * calc_per_rank;
        iter = _iterations;
    }
    MPI_Bcast(&pop_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_per_rank, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);

    //std::cout << "Rank: " << rank << "        pop_size: " << pop_size << std::endl;
    //std::cout << "Rank: " << rank << "        calc_per_rank: " << calc_per_rank << std::endl;
    
    //MPI stuff for reduce
    struct bestmoe_s{
        float genotype1;
        float genotype2;
        double fitness;
    };
    int blocklen_2[2]={2,1};
    MPI_Aint displ_2[2]={0,8};
    MPI_Datatype types_2[2]={MPI_FLOAT,MPI_DOUBLE};  //[ffff ffff dddddddd]
    MPI_Datatype MoeType;
    MPI_Type_create_struct(
        2,
        blocklen_2,
        displ_2,
        types_2,
        &MoeType
    );
    MPI_Type_commit(&MoeType);
    //END MPI stuff for gather

    struct particlescatter{
      float genotype1;
      float genotype2;
      float local_best_x;
      float local_best_y;
      float global_best_x;
      float global_best_y;
      float velocity_x;
      float velocity_y;
      float weight;
      float coeff1;
      float coeff2;
      double fitness;
      double local_best_fitness;
      double global_best_fitness;
      int idx;
    };


    int blocklen[3]={11,3,1};
    MPI_Aint displ[3]={0,44,68};
    MPI_Datatype types[3]={MPI_FLOAT,MPI_DOUBLE,MPI_INT};
    MPI_Datatype ScatterType;
    MPI_Type_create_struct(
      3,
      blocklen,
      displ,
      types,
      &ScatterType
  );
    MPI_Type_commit(&ScatterType);

    particlescatter particles_to_scatter[pop_size];
    particlescatter particle_to_recieve[calc_per_rank];

    if(rank == 0) 
    {
        
        /*this->init( _iterations );
        pop_size = m_population.size();
        calc_per_rank = (int) pop_size / size;
        rest = pop_size - (size - 1) * calc_per_rank;*/

        for(int j = 0; j < pop_size; ++j){
            particles_to_scatter[j].genotype1 = m_population[j].genotype[0];
            particles_to_scatter[j].genotype2 = m_population[j].genotype[1];
            particles_to_scatter[j].local_best_x = m_best_genotypes[j].genotype[0];
            particles_to_scatter[j].local_best_y = m_best_genotypes[j].genotype[1];
            particles_to_scatter[j].local_best_fitness = m_best_genotypes[j].fitness;
            particles_to_scatter[j].global_best_x = Algorithm<GenotypeType>::m_bestMoe.genotype[0];
            particles_to_scatter[j].global_best_y = Algorithm<GenotypeType>::m_bestMoe.genotype[1];
            particles_to_scatter[j].global_best_fitness = Algorithm<GenotypeType>::m_bestMoe.fitness;
            particles_to_scatter[j].velocity_x = m_velocities[j][0];
            particles_to_scatter[j].velocity_y = m_velocities[j][1];
            particles_to_scatter[j].fitness = m_population[j].fitness;
            particles_to_scatter[j].weight = m_weight;
            particles_to_scatter[j].coeff1 = m_coef1;
            particles_to_scatter[j].coeff2 = m_coef2;
            particles_to_scatter[j].idx = j;
        }
    }


    //int part_on_slave = (int) m_population.size() / (size-1); //number of particles to be evalutated on each slave
    int part_on_master = pop_size-(size-1)*calc_per_rank;
    if(rank == 0) calc_per_rank = part_on_master;

    //std::cout << "Pre-Scatter" << std::endl;

    MPI_Scatter(&particles_to_scatter, calc_per_rank, ScatterType, &particle_to_recieve, calc_per_rank, ScatterType, 0, MPI_COMM_WORLD);

    //std::cout << "Post-Scatter" << std::endl;

    //std::vector< Moe<GenotypeType> > m_population_v2[calc_per_rank];
    //std::vector< Moe<GenotypeType> > m_best_genotypes_v2[calc_per_rank];

    
    //int calc_per_rank;
    //calc_per_rank = (int) m_population.size() / size;
    
    //MPI_Scatter(particles_to_scatter, 4, MPI_DOUBLE, particle_to_recieve, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //float global_best = -1000000000;

    //this->init( _iterations, calc_per_rank);

    for( unsigned int i = 0; i < iter; ++i )
        {
            //std::cout << "Iteration " << i << " started now!" << std::endl;

            bestmoe_s best_of_iter[calc_per_rank];
            //best_of_iter.reserve(calc_per_rank);
            //std::cout << "iterations started" << std::endl;
            //for every process: update, eval, local_bestmoe
            for( unsigned int j = 0; j < calc_per_rank; ++j )
            {
                //for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )
                //{

                    float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator ),
                            r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );
                    
                    //velocity update x direction
                    particle_to_recieve[j].velocity_x *= particle_to_recieve[j].weight;
                    particle_to_recieve[j].velocity_x += particle_to_recieve[j].coeff1*r1*(particle_to_recieve[j].local_best_x - particle_to_recieve[j].genotype1);
                    particle_to_recieve[j].velocity_x += particle_to_recieve[j].coeff2*r2*(particle_to_recieve[j].global_best_x - particle_to_recieve[j].genotype1);
                    particle_to_recieve[j].genotype1 += particle_to_recieve[j].velocity_x;

                    //velocity update y direction
                    particle_to_recieve[j].velocity_y *= particle_to_recieve[j].weight;
                    particle_to_recieve[j].velocity_y += particle_to_recieve[j].coeff1*r1*(particle_to_recieve[j].local_best_y - particle_to_recieve[j].genotype2);
                    particle_to_recieve[j].velocity_y += particle_to_recieve[j].coeff2*r2*(particle_to_recieve[j].global_best_y - particle_to_recieve[j].genotype2);
                    particle_to_recieve[j].genotype2 += particle_to_recieve[j].velocity_y;

                    //std::cout << "Pos and Vel updated" << std::endl;
                    //std::cout << "x: " << particle_to_recieve[j].genotype1 << std::endl;

                    /*m_velocities[j][k] *= m_weight;
                    m_velocities[j][k] += m_coef1*r1*( m_best_genotypes[j].genotype[k] - m_population[j].genotype[k] );
                    m_velocities[j][k] += m_coef2*r2*( Algorithm<GenotypeType>::m_bestMoe.genotype[k] - m_population[j].genotype[k] );
                    m_population[j].genotype[k] += m_velocities[j][k];*/
                    //checks if genotype actual dimension is still in provided search space

                    if( particle_to_recieve[j].genotype1 != std::max( std::min( particle_to_recieve[j].genotype1, NumericAlgorithm<GenotypeType>::m_range[1]), NumericAlgorithm<GenotypeType>::m_range[0] ) ) 
                    {
                    //if(particle_to_recieve[j].genotype1 < NumericAlgorithm<GenotypeType>::m_range[0] || particle_to_recieve[j].genotype1 > NumericAlgorithm<GenotypeType>::m_range[1])
                    //{
                        std::vector<GenotypeType> temp_pos, temp_vel;
                        temp_pos.reserve(2);
                        temp_vel.reserve(2);
                        temp_pos = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                        particle_to_recieve[j].genotype1 = temp_pos[0];
                        particle_to_recieve[j].genotype2 = temp_pos[1];
                        temp_vel = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                        particle_to_recieve[j].velocity_x = temp_vel[0];
                        particle_to_recieve[j].velocity_y = temp_vel[1];
                        //break;
                    }
                    else if( particle_to_recieve[j].genotype2 != std::max( std::min( particle_to_recieve[j].genotype2, NumericAlgorithm<GenotypeType>::m_range[1]), NumericAlgorithm<GenotypeType>::m_range[0] ) ) 
                    {
                    //else if(particle_to_recieve[j].genotype2 < NumericAlgorithm<GenotypeType>::m_range[0] || particle_to_recieve[j].genotype2 > NumericAlgorithm<GenotypeType>::m_range[1])
                    //{
                        std::vector<GenotypeType> temp_pos, temp_vel;
                        temp_pos.reserve(2);
                        temp_vel.reserve(2);
                        temp_pos = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                        particle_to_recieve[j].genotype1 = temp_pos[0];
                        particle_to_recieve[j].genotype2 = temp_pos[1];
                        temp_vel = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                        particle_to_recieve[j].velocity_x = temp_vel[0];
                        particle_to_recieve[j].velocity_y = temp_vel[1];
                        //break;
                    }

                    //std::cout << "x: " << particle_to_recieve[j].genotype1 << std::endl;

                    //std::cout << "if out of bounds" << std::endl;

                    /*if( m_population[j].genotype[k] != std::max( std::min( m_population[j].genotype[k], NumericAlgorithm<GenotypeType>::m_range[1]), NumericAlgorithm<GenotypeType>::m_range[0] ) )
                    {
                        m_population[j].genotype = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                        m_velocities[j] = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                        break;
                    }*/
                //}

                Moe<GenotypeType> temp;
                //temp.genotype = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                //std::cout << "temp.genotype[0]: " << temp.genotype[0] << std::endl;
                std::vector<GenotypeType> genotype_mine;
                genotype_mine.reserve( 2 );
                genotype_mine.push_back(particle_to_recieve[j].genotype1);
                genotype_mine.push_back(particle_to_recieve[j].genotype2);
                //genotype_mine[0] = particle_to_recieve[j].genotype1;
                //genotype_mine[1] = particle_to_recieve[j].genotype2;
                //temp.genotype[0] = particle_to_recieve[j].genotype1;
                //temp.genotype[1] = particle_to_recieve[j].genotype2;
                temp.genotype = genotype_mine;
                temp.fitness = particle_to_recieve[j].fitness;
                particle_to_recieve[j].fitness = Algorithm<GenotypeType>::m_fitnessFunction( temp );
                //std::cout << "fitness: " << particle_to_recieve[j].fitness << std::endl;

                //std::cout << "updated fitness" << std::endl;

                //m_population[j].fitness = Algorithm<GenotypeType>::m_fitnessFunction( m_population[j] );
                if( particle_to_recieve[j].fitness < particle_to_recieve[j].local_best_fitness)
                {
                    //int GLOBAL_BEST_MOE = 0;
                    particle_to_recieve[j].local_best_x = particle_to_recieve[j].genotype1;
                    particle_to_recieve[j].local_best_y = particle_to_recieve[j].genotype2;
                    particle_to_recieve[j].local_best_fitness = particle_to_recieve[j].fitness;
                    //if( m_best_genotypes[j].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness )
                        //Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[j];
                    /*if(particle_to_recieve[j].fitness < particle_to_recieve[j].global_best_fitness)
                    {
                        particle_to_recieve[j].global_best_x = particle_to_recieve[j].genotype1;
                        particle_to_recieve[j].global_best_y = particle_to_recieve[j].genotype2;
                        particle_to_recieve[j].global_best_fitness = particle_to_recieve[j].fitness;

                        //NEEDDD TOOO IMPLEMENT SEND GLOBAL BEST TO RANK 0 AND BACK
                    }*/
                }

                //std::cout << "updated local and global best" << std::endl;

                bestmoe_s g_best;
                g_best.genotype1 = particle_to_recieve[j].genotype1;
                g_best.genotype2 = particle_to_recieve[j].genotype2;
                g_best.fitness = particle_to_recieve[j].fitness;

                best_of_iter[j] = g_best;

                /*bestmoe_s best_moe_recieved[pop_size]; //IE ANZAHL PROZESSOREN

                //std::cout << "Pre-Gather" << std::endl;

                //const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,MPI_Datatype recvtype, int root, MPI_Comm comm
                MPI_Gather(&g_best, 1, MoeType, best_moe_recieved, 1, MoeType, 0, MPI_COMM_WORLD);

                //std::cout << "Post-Gather" << std::endl;
                
                bestmoe_s send_back;
                send_back.genotype1 = 0;
                send_back.genotype2 = 0;
                send_back.fitness = 0;
                if(rank == 0) 
                {
                    for(auto part: best_moe_recieved)
                    {
                        if(part.fitness > Algorithm<GenotypeType>::m_bestMoe.fitness)
                        {
                            Algorithm<GenotypeType>::m_bestMoe.genotype[0] = part.genotype1;
                            Algorithm<GenotypeType>::m_bestMoe.genotype[1] = part.genotype2;
                            Algorithm<GenotypeType>::m_bestMoe.fitness = part.fitness;

                            send_back.genotype1 = part.genotype1;
                            send_back.genotype2 = part.genotype2;
                            send_back.fitness = part.fitness;
                        }
                    }
                    
                }*/

                //std::cout << "Pre-Bcast" << std::endl;
                //MPI_Bcast(&send_back, 1, MoeType, 0, MPI_COMM_WORLD);
                //std::cout << "Post-Bcast" << std::endl;
                //particle_to_recieve[j].global_best_x = send_back.genotype1;
                //particle_to_recieve[j].global_best_y = send_back.genotype2;
                //particle_to_recieve[j].global_best_fitness = send_back.fitness;
                /*if( m_population[j].fitness > m_best_genotypes[j].fitness ) 
                {
                    m_best_genotypes[j] = m_population[j];
                }*/
                /*if( m_population[j].fitness > particle_to_recieve[j].local_best_fitness)
                {
                    //int GLOBAL_BEST_MOE = 0;
                    m_best_genotypes[j] = m_population[j];
                    //if( m_best_genotypes[j].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness )
                        //Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[j];
                }*/
            }
            //Now we want to broadcast the best moe
            bestmoe_s best_moe_recieved[pop_size]; //IE ANZAHL PROZESSOREN
            //best_moe_recieved.reserve(pop_size);

            //std::cout << "Pre-Gather" << std::endl;

            //const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,MPI_Datatype recvtype, int root, MPI_Comm comm
            MPI_Gather(&best_of_iter, calc_per_rank, MoeType, best_moe_recieved, calc_per_rank, MoeType, 0, MPI_COMM_WORLD);

            //std::cout << "Post-Gather" << std::endl;
                
            bestmoe_s send_back;
            //send_back.genotype1 = 0;
            //send_back.genotype2 = 0;
            //send_back.fitness = 0;
            if(rank == 0) 
            {
                for(auto part: best_moe_recieved)
                {
                    send_back.genotype1 = Algorithm<GenotypeType>::m_bestMoe.genotype[0];
                    send_back.genotype2 = Algorithm<GenotypeType>::m_bestMoe.genotype[1];
                    send_back.fitness = Algorithm<GenotypeType>::m_bestMoe.fitness;
                    if(part.fitness < Algorithm<GenotypeType>::m_bestMoe.fitness)
                    {
                        Algorithm<GenotypeType>::m_bestMoe.genotype[0] = part.genotype1;
                        Algorithm<GenotypeType>::m_bestMoe.genotype[1] = part.genotype2;
                        Algorithm<GenotypeType>::m_bestMoe.fitness = part.fitness;

                        send_back.genotype1 = part.genotype1;
                        send_back.genotype2 = part.genotype2;
                        send_back.fitness = part.fitness;
                    }
                }
                    
            }
            MPI_Bcast(&send_back, 1, MoeType, 0, MPI_COMM_WORLD);

            for( unsigned int j = 0; j < calc_per_rank; ++j )
            {
                particle_to_recieve[j].global_best_x = send_back.genotype1;
                particle_to_recieve[j].global_best_y = send_back.genotype2;
                particle_to_recieve[j].global_best_fitness = send_back.fitness;
            }
        }

    //std::cout << "rank: " << rank << " got here!" << std::endl;
















    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*for( unsigned int i = 0; i < m_iterations; ++i )
    {   
        //if(rank < (size - 1)) 
        //{
            for(unsigned int j = (unsigned) rank * calc_per_rank; j < (rank + 1) * calc_per_rank; j++) 
            {
                for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )
                {
                    float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator );
                    float   r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

                    std::cout << "r1: " << r1 << "            r2: " << r2 << std::endl;

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
                }
            } 
        //}
        //else {
            for(unsigned int j = (unsigned) rank * calc_per_rank; j < m_population.size(); j++) {
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
                }
            }
        //}

        if(rank == 0) 
        {
            unsigned k = 0;
            max_of_iter = m_best_genotypes[0].fitness;
            for(unsigned j = 1; j < m_population.size(); ++j) {
                if(m_best_genotypes[j].fitness >= max_of_iter) {
                    max_of_iter = m_best_genotypes[j].fitness;
                    k = j;
                }
            }
            if( max_of_iter > Algorithm<GenotypeType>::m_bestMoe.fitness )
                { 
                    std::cout << "max_of_iter: " << max_of_iter << "         m_best_genotypes[j].fitness: " << m_best_genotypes[k].fitness << "           Algorithm<GenotypeType>::m_bestMoe.fitness: " << Algorithm<GenotypeType>::m_bestMoe.fitness << std::endl;
                    Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[k];
                }
        }
    }
}*/

} 


//take out line 118 (fitness evaluation) and parallelize this process. 
//Insert values into a vector of fitness values.
//create a function with rest of population for loop, that takes this vector as an input and updates postion, velocity, local best and global best.

//evaluate(std::vector<double> fitness_values);

} //end namespace Moe