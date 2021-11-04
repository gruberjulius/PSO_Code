/*
Trying out smaller packets with MPI_DATAYPE struct
FLorian
*/



#pragma once

#include "NumericAlgorithmImpl.hpp"
#include "../parameters/PSParameters.hpp"

//EDIT START
#include <queue>
#include <mpi.h>
#include <iostream>
//EDIT END

namespace moe
{

template <typename GenotypeType>
class ParticleSwarmPA_send_move_smaller_packets : public NumericAlgorithm<GenotypeType>
{
    public:
        ParticleSwarmPA_send_move_smaller_packets( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions = 1, std::vector<GenotypeType> _range = { std::numeric_limits<GenotypeType>::lowest() , std::numeric_limits<GenotypeType>::max() });
        ParticleSwarmPA_send_move_smaller_packets( const PSParameters<GenotypeType>& _parameters );

        void run( unsigned int _iterations, int& argc, char * argv[] ); //was override
        void run (unsigned int _iterations) override;

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
ParticleSwarmPA_send_move_smaller_packets<GenotypeType>::ParticleSwarmPA_send_move_smaller_packets( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions, std::vector<GenotypeType> _range )
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
ParticleSwarmPA_send_move_smaller_packets<GenotypeType>::ParticleSwarmPA_send_move_smaller_packets( const PSParameters<GenotypeType>& _parameters )
:ParticleSwarmPA_send_move_smaller_packets<GenotypeType>( _parameters.moesPerGen, _parameters.inertia, _parameters.coef1, _parameters.coef2, _parameters.dimensions, _parameters.range )
{

}

template <typename GenotypeType>
void ParticleSwarmPA_send_move_smaller_packets<GenotypeType>::init( unsigned int _iterations )
{
    m_iterations = _iterations;

    double min = 1000.0;
    unsigned int    index = 0,
                    count = 0;

    for( auto& moe : m_population )
    {
        moe.genotype            = NumericAlgorithm<GenotypeType>::getRandomGenotype();
        moe.fitness             = Algorithm<GenotypeType>::m_fitnessFunction( moe );

        m_velocities[count]       = NumericAlgorithm<GenotypeType>::getRandomGenotype();

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
}

template <typename GenotypeType>
void ParticleSwarmPA_send_move_smaller_packets<GenotypeType>::run( unsigned int _iterations){std::cout<<"The wrong \" run \" function has been called in ParticleSwarmPA_send_move_smaller_packets.cpp"<<std::endl;}
template <typename GenotypeType>
void ParticleSwarmPA_send_move_smaller_packets<GenotypeType>::run( unsigned int _iterations, int& argc, char * argv[])
{
    int rank, numproces;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);

//FLORIAN EDIT START
    struct particlescatter{
      float genotype1;
      float genotype2;
      int idx;
    };

    int blocklen[2]={2,1};
    MPI_Aint displ[2]={0,8};
    MPI_Datatype types[2]={MPI_FLOAT,MPI_INT};
    MPI_Datatype ScatterType;
    MPI_Type_create_struct(
      2,
      blocklen,
      displ,
      types,
      &ScatterType
  );
    MPI_Type_commit(&ScatterType);

    particlescatter particles_to_scatter[(numproces-1)];
    particlescatter particle_to_recieve;

    /*double *particles_to_scatter = new double [(numproces-1)*4];
    double *particle_to_recieve = new double [4];*/
    //double particles_to_scatter[(numproces-1)*3];
    // particle_to_recieve[3];
    if(rank == 0){
        this->init( _iterations );
        particles_to_scatter[0].genotype1 = 69;
        particles_to_scatter[0].genotype2 = 69;
        particles_to_scatter[0].idx = 69;
        for(int j = 1; j<(numproces); ++j){
            particles_to_scatter[j].genotype1 = m_population[j-1].genotype[0];
            particles_to_scatter[j].genotype2 = m_population[j-1].genotype[1];
            particles_to_scatter[j].idx = j-1;
        }
//        std::cout << "Pre scatter\n";

        /*for ( int i = 0; i < 4*(numproces); ++i){
            std::cout << particles_to_scatter[i] << " ";
        }
        std::cout<<"this was master ^ " << std::endl;//*/
    }
    MPI_Scatter(particles_to_scatter, 1, ScatterType, &particle_to_recieve, 1, ScatterType, 0, MPI_COMM_WORLD);
//FLORIAN EDIT END
    //delete[] particles_to_scatter;
    if(rank == 0){
        //initialization
//        std::cout << "Post scatter\n";

        int num_particles = m_population.size();

        std::queue<int> particlequeue;//, available_threads; //queue-s of particle indexes and available thread ranks
        //initialization of particle queue with numproces-1,...,#particles
        for(int i = numproces-1; i < num_particles; ++i){ //numproces-1 because we already scattered
            particlequeue.push(i);
        }

        //unsigned int step = 0;
        unsigned int curr_sum_iters = 0; //current # of particles that have been evaluated (used for step, to have a rough bound on limitations)
//        std::cout<<"Checkpoint 1: pre-while loop"<<std::endl;

        float our_current_particle[2];
        double fitness_f;

        //std::cout<<"# Iterations: "<<m_iterations<<std::endl;

        double modified_iterations = num_particles * m_iterations;
        while(curr_sum_iters < modified_iterations)//was step < m_iterations
        {
//                std::cout<<"Checkpoint 2: while loop step: "<< step <<std::endl;

                MPI_Status status;
                MPI_Recv(&fitness_f, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//need MPI_Datatype of Particle

                int current_idx = status.MPI_TAG; //particle index recieved in tag. Send particle idx + 1, to know when to terminate
                //std::cout << current_idx << "is da current index" << std::endl;
                m_population[current_idx].fitness = fitness_f;
                //std::cout<< "recieved fitness: " << our_current_particle[2] << ", particle number: " << current_idx << std::endl;
                //std::cout << "Current idx: " << current_idx <<std::endl;
                particlequeue.push(current_idx);
                //available_threads.push(status.MPI_SOURCE);
                //std::cout <<"Particle queue is empty " << particlequeue.empty() << std::endl;

                int new_idx = particlequeue.front(); //get new particle idx
                particlequeue.pop();//deletes accessed particle
                //std::cout<<"Particle que size" << particlequeue.size()<<std::endl;

                Moe<GenotypeType> current_particle = m_population[new_idx];
                our_current_particle[0] = current_particle.genotype[0];
                our_current_particle[1] = current_particle.genotype[1];

                //std::cout<<"Master: Pre send" << std::endl;
                double tag = new_idx + 1.;
                MPI_Send(our_current_particle, 2, MPI_FLOAT, status.MPI_SOURCE, tag, MPI_COMM_WORLD);
                //std::cout<<"Master: Post send" << std::endl;

                //std::cout<<"Checkpoint 3: this sh!t is poppin"<<std::endl;
                int j = current_idx;
                for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )// 2D -> k = 0, 1; UNROLL ?
                {
                    float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator ),
                            r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

                    m_velocities[j][k] *= m_weight;
                    m_velocities[j][k] += m_coef1*r1*( m_best_genotypes[j].genotype[k] - m_population[j].genotype[k] );
                    m_velocities[j][k] += m_coef2*r2*( Algorithm<GenotypeType>::m_bestMoe.genotype[k] - m_population[j].genotype[k] );

                    m_population[j].genotype[k] += m_velocities[j][k];

                    //checks if current position is still in provided search space
                    if( m_population[j].genotype[k] != std::max( std::min( m_population[j].genotype[k], NumericAlgorithm<GenotypeType>::m_range[1]), NumericAlgorithm<GenotypeType>::m_range[0] ) )
                    {
                        m_population[j].genotype = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                        m_velocities[j] = NumericAlgorithm<GenotypeType>::getRandomGenotype();
                        break;
                    }
                }
                //std::cout<<"Checkpoint 4: post-4 loop"<<std::endl;

                if( m_population[current_idx].fitness < m_best_genotypes[current_idx].fitness )
                {
                    m_best_genotypes[current_idx] = m_population[current_idx];
                    if( m_best_genotypes[current_idx].fitness < Algorithm<GenotypeType>::m_bestMoe.fitness )
                        Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[current_idx];
                }
                //std::cout<<"Checkpoint 8: full while loop iteration"<<std::endl;

                curr_sum_iters += 1;
                //step = curr_sum_iters / num_particles;
//                std::cout<<"Step: "<<step<< ", curr_sum_iter: " << curr_sum_iters << ", num_particles: " << num_particles << ", # iter: " << m_iterations<<std::endl;
        }
        //std::cout<<"Checkpoint 9: while loop exited"<<std::endl;
        //loop to kill all slave processes
        for(int id = 1; id < numproces; ++id){
            MPI_Status status;
            double fitness_f;
            MPI_Recv(&fitness_f, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int particle_number = status.MPI_TAG;
            m_population[particle_number].fitness = fitness_f;
            if( m_population[particle_number].fitness < m_best_genotypes[particle_number].fitness )//OPTIMISE WITH OR WITHOUT CHECK CONVERGENCE
                {
                    m_best_genotypes[particle_number] = m_population[particle_number];
                    if( m_best_genotypes[particle_number].fitness < Algorithm<GenotypeType>::m_bestMoe.fitness )
                        Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[particle_number];
                }
            //std::cout<<"terminator send to: "<<status.MPI_SOURCE<<std::endl;
            MPI_Send(our_current_particle, 2, MPI_FLOAT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        }

        //std::cout<<"Checkpoint 10: terminator for-loop exited"<<std::endl;

    }


    if(rank != 0){ //slave process

        //double *current_position = new double[3];
        float current_position[2];
        //double time;
        //double t1 = MPI_Wtime();

        /*for ( int i = 0; i < 4; ++i){
            std::cout << particle_to_recieve[i] << " ";
        }
        std::cout<<std::endl;//*/
//FLORIAN EDIT START
        int curr_particle_number;

        std::vector<GenotypeType> gen = {particle_to_recieve.genotype1, particle_to_recieve.genotype2};

        Moe<GenotypeType> current_particle = {gen, 42.}; //current fitness doesn't matter for function
        double fitness_f = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );

        //std::cout << "Rank: " << rank << " pre-first send\n";
        MPI_Send(&fitness_f, 1, MPI_DOUBLE, 0, particle_to_recieve.idx, MPI_COMM_WORLD);
//FLORIAN EDIT END
        //std::cout << "Rank: " << rank << " first send done\n";
        while(1){

            MPI_Status status;
            MPI_Recv(current_position, 2, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//need MPI_Datatype of Particle
            //std::cout<<"Recieved tag: " << status.MPI_TAG << std::endl;
            if(status.MPI_TAG == 0){
                //std::cout<<"BREAKED rank: " <<rank<<std::endl;
                break;
            }

            int particle_number = status.MPI_TAG - 1;//decr to make tag into actual particle index

            //std::cout<< "current fitness: " << current_position[2] << ", rank: " << rank << std::endl;
            gen = {current_position[0], current_position[1]};
            current_particle = {gen, 42.};
            fitness_f = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );

            //std::cout << "Rank: " << rank << ", pre send" << std::endl;
            MPI_Send(&fitness_f, 1, MPI_DOUBLE, 0, particle_number, MPI_COMM_WORLD);
            //std::cout<< "new fitness: " << current_position[2] << ", particle number: " << particle_number << std::endl;
            //double t2 = MPI_Wtime();
            //time = t2 - t1;

        }
        //std::cout << "Rank: " << rank << " ran for: " << time << "seconds\n";
        //delete[] current_position;

    }

   //delete[] particle_to_recieve;


}

}
