/*
Florian
transforming the papso_precompute_vel into non-blocking  sends to slaves
-> no Iscatter
Find changes:
//FLORIAN EDIT START

//FLORIAN EDIT END

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
class ParticleSwarmPA_prec_vel_nonblc : public NumericAlgorithm<GenotypeType>
{
    public:
        ParticleSwarmPA_prec_vel_nonblc( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions = 1, std::vector<GenotypeType> _range = { std::numeric_limits<GenotypeType>::lowest() , std::numeric_limits<GenotypeType>::max() });
        ParticleSwarmPA_prec_vel_nonblc( const PSParameters<GenotypeType>& _parameters );

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
ParticleSwarmPA_prec_vel_nonblc<GenotypeType>::ParticleSwarmPA_prec_vel_nonblc( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions, std::vector<GenotypeType> _range )
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
ParticleSwarmPA_prec_vel_nonblc<GenotypeType>::ParticleSwarmPA_prec_vel_nonblc( const PSParameters<GenotypeType>& _parameters )
:ParticleSwarmPA_prec_vel_nonblc<GenotypeType>( _parameters.moesPerGen, _parameters.inertia, _parameters.coef1, _parameters.coef2, _parameters.dimensions, _parameters.range )
{

}

template <typename GenotypeType>
void ParticleSwarmPA_prec_vel_nonblc<GenotypeType>::init( unsigned int _iterations )
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
void ParticleSwarmPA_prec_vel_nonblc<GenotypeType>::run( unsigned int _iterations)
{
    int rank, numproces;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);

    /*double *particles_to_scatter = new double [(numproces-1)*4];
    double *particle_to_recieve = new double [4];*/
    //int scattersize = (numproces-1)*6;
    double particles_to_scatter[(numproces-1)*8];
    double particle_to_recieve[8];
    if(rank == 0){
        this->init( _iterations );
        //initialize first "entry" of particles_to_scatter as dummy
        //dummy will be recieved (scatter) by master, where we don't need data
        for(int i = 0; i < 8; ++i){
            particles_to_scatter[i] = 69;
        }
        /*//loop unrolled
        particles_to_scatter[0] = 69;
        particles_to_scatter[1] = 69;
        particles_to_scatter[2] = 69;
        particles_to_scatter[3] = 69;
        particles_to_scatter[4] = 69;
        particles_to_scatter[5] = 69;
        particles_to_scatter[6] = 69;
        particles_to_scatter[7] = 69;*/
        for(int j = 1; j<(numproces); ++j){
            particles_to_scatter[j*8] = m_population[j-1].genotype[0];//x
            particles_to_scatter[j*8 + 1] = m_population[j-1].genotype[1];//y
            particles_to_scatter[j*8 + 2] = m_best_genotypes[j-1].fitness;//best val
            particles_to_scatter[j*8 + 3] = m_velocities[j-1][0];//u
            particles_to_scatter[j*8 + 4] = m_velocities[j-1][1];//v
            particles_to_scatter[j*8 + 5] = m_best_genotypes[j-1].genotype[0];//best_x
            particles_to_scatter[j*8 + 6] = m_best_genotypes[j-1].genotype[1];//best_y
            particles_to_scatter[j*8 + 7] = j-1;//particle index
        }
//        //std::cout << "Pre scatter\n";

        /*for ( int i = 0; i < 4*(numproces); ++i){
            //std::cout << particles_to_scatter[i] << " ";
        }
        //std::cout<<"this was master ^ " << std::endl;//*/
    }
    //MPI_ScatterV MAYBE
    //FLORIAN EDIT START
    MPI_Request unused_request;
    MPI_Scatter(particles_to_scatter, 8, MPI_DOUBLE, particle_to_recieve, 8, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //OLD: MPI_Scatter(particles_to_scatter, 8, MPI_DOUBLE, particle_to_recieve, 8, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //FLORIAN EDIT END
    //delete[] particles_to_scatter;
    if(rank == 0){
        //initialization
//        //std::cout << "Post scatter\n";

        int num_particles = m_population.size();

        std::queue<int> particlequeue;//, available_threads; //queue-s of particle indexes and available thread ranks
        //initialization of particle queue with numproces-1,...,#particles
        for(int i = numproces-1; i < num_particles; ++i){ //numproces-1 because we already scattered
            particlequeue.push(i);
        }

        //unsigned int step = 0;
        unsigned int curr_sum_iters = 0; //current # of particles that have been evaluated (used for step, to have a rough bound on limitations)
//        //std::cout<<"Checkpoint 1: pre-while loop"<<std::endl;

        double our_send_particle[7];
        //[x, y, local_best_fitness, u, v, local_best_x, local_best_y]
        double our_recv_particle[6];
        //[local_best_x, local_best_y, fitness, u, v, r2]

        //std::cout<<"# Iterations: "<<m_iterations<<std::endl;

        double modified_iterations = num_particles * m_iterations;
        while(curr_sum_iters < modified_iterations)//was step < m_iterations
        {
//                //std::cout<<"Checkpoint 2: while loop step: "<< step <<std::endl;

                MPI_Status status;
                MPI_Recv(our_recv_particle, 6, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//need MPI_Datatype of Particle

                //particle index recieved in tag. Send particle idx + 1, to know when to terminate
                int current_idx = status.MPI_TAG;

                //std::cout<< "recieved fitness: " << our_recv_particle[2] << ", particle number: " << current_idx << std::endl;
                //std::cout << "Current idx: " << current_idx <<std::endl;
                particlequeue.push(current_idx);
                //available_threads.push(status.MPI_SOURCE);
                //std::cout <<"Particle queue is empty " << particlequeue.empty() << std::endl;
                int new_idx = particlequeue.front(); //get new particle idx
                particlequeue.pop();//deletes accessed particle
                //std::cout<<"Particle que size" << particlequeue.size()<<std::endl;
                //std::cout<<"Checkpoint 3: this sh!t is poppin"<<std::endl;

                Moe<GenotypeType> current_particle = m_population[new_idx];
                //fill send buffer with:
                //current particle x and y
                our_send_particle[0] = current_particle.genotype[0];
                our_send_particle[1] = current_particle.genotype[1];
                //best local fitness
                our_send_particle[2] = m_best_genotypes[current_idx].fitness;
                //current velocities
                our_send_particle[3] = m_velocities[new_idx][0];
                our_send_particle[4] = m_velocities[new_idx][1];
                //best local x and y
                our_send_particle[5] = m_best_genotypes[current_idx].genotype[0];
                our_send_particle[6] = m_best_genotypes[current_idx].genotype[1];

                //std::cout<<"Master: Pre send" << std::endl;
                double tag = new_idx + 1.;
                //FLORIAN EDIT START

                MPI_Isend(our_send_particle, 7, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD,&unused_request);
                //OLD: MPI_Send(our_send_particle, 7, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD);
                //FLORIAN EDIT START
                //std::cout<<"Master: Post send" << std::endl;

                //update fitness and velocities
                m_population[current_idx].fitness = our_recv_particle[2];
                m_velocities[current_idx][0] = our_recv_particle[3];
                m_velocities[current_idx][1] = our_recv_particle[4];
                float r2 = our_recv_particle[5];

                //checks if function evaluation from slave is new best
                if(our_recv_particle[2] < m_best_genotypes[current_idx].fitness){
                //update m_best_genotypes(local best) with current particle
                    m_best_genotypes[current_idx].genotype[0] = our_recv_particle[0];
                    m_best_genotypes[current_idx].genotype[1] = our_recv_particle[1];
                    m_best_genotypes[current_idx].fitness = our_recv_particle[2];
                    //checks if function evaluation from slave is new global best
                    if( m_best_genotypes[current_idx].fitness < Algorithm<GenotypeType>::m_bestMoe.fitness )
                        Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[current_idx];
                }

                int j = current_idx;//finish update of previous particle
                for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )// 2D -> k = 0, 1; UNROLL ?
                {

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




                //std::cout<<"Checkpoint 8: full while loop iteration"<<std::endl;

                curr_sum_iters += 1;
                //step = curr_sum_iters / num_particles;
                ////std::cout<<"Step: "<<step<< ", curr_sum_iter: " << curr_sum_iters << ", num_particles: " << num_particles << ", # iter: " << m_iterations<<std::endl;
        }
        //std::cout<<"Checkpoint 9: while loop exited"<<std::endl;

        //loop to kill all slave processes
        for(int id = 1; id < numproces; ++id){
            MPI_Status status;
            MPI_Recv(our_recv_particle, 6, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int particle_number = status.MPI_TAG;
            if(our_recv_particle[2] < m_best_genotypes[particle_number].fitness){//OPTIMIZE WITH OR WITHOUT THIS
            //update m_best_genotypes(local best) with current particle
                m_best_genotypes[particle_number].genotype[0] = our_recv_particle[0];
                m_best_genotypes[particle_number].genotype[1] = our_recv_particle[1];
                m_best_genotypes[particle_number].fitness = our_recv_particle[2];
                //checks if function evaluation from slave is new global best
                if( m_best_genotypes[particle_number].fitness < Algorithm<GenotypeType>::m_bestMoe.fitness )
                    Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[particle_number];
            }
            //std::cout<<"terminator send to: "<<status.MPI_SOURCE<<std::endl;
            //FLORIAN EDIT START
            MPI_Isend(our_recv_particle, 7, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD,&unused_request);
            //FLORIAN EDIT END
        }

        //std::cout<<"Checkpoint 10: terminator for-loop exited"<<std::endl;

    }


    if(rank != 0){ //slave process

        //double *current_position = new double[3];
        double current_position[8];
        //double time;
        //double t1 = MPI_Wtime();

        /*for ( int i = 0; i < 4; ++i){
            //std::cout << particle_to_recieve[i] << " ";
        }
        //std::cout<<std::endl;//*/
        //std::cout<<"m_weight" << m_coef1<<std::endl;

        std::copy(particle_to_recieve, particle_to_recieve + 8, current_position);
        int curr_particle_number = (int)(double)particle_to_recieve[7];

        double best_val = current_position[2];//save best value in temp

        std::vector<double> gen = {current_position[0], current_position[1]};

        Moe<double> current_particle = {gen, current_position[2]};
        current_position[2] = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );

        //temp variables for velocities and best positions
        double vel_x, vel_y, best_x, best_y;
        vel_x = current_position[3];
        vel_y = current_position[4];

        if(best_val >= current_position[2]){//no change in bestval
            best_x = current_position[5];
            best_y = current_position[6];
        }else{ //new fitness is local max
            best_x = current_position[0];
            best_y = current_position[1];
        }
        float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator ),
                r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

        vel_x *= m_weight;
        vel_x += m_coef1*r1*( best_x - current_position[3] );
        vel_y *= m_weight;
        vel_y += m_coef1*r1*( best_y - current_position[4] );

        current_position[3] = vel_x;
        current_position[4] = vel_y;
        current_position[5] = r2;

        //std::cout << "Rank: " << rank << " pre-first send\n";
        //std::cout << rank << std::endl;
        //std::cout << vel_x << " " <<vel_y << std::endl;
        MPI_Send(current_position, 6, MPI_DOUBLE, 0, curr_particle_number, MPI_COMM_WORLD);
        //std::cout << "Rank: " << rank << " first send done\n";
        while(1){

            MPI_Status status;
            MPI_Recv(current_position, 7, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//need MPI_Datatype of Particle
            //std::cout<<"Recieved tag: " << status.MPI_TAG << std::endl;
            if(status.MPI_TAG == 0){
                //std::cout<<"BREAKED rank: " <<rank<<std::endl;
                break;
            }

            curr_particle_number = status.MPI_TAG - 1;//decr to make tag into actual particle index

            //std::cout<< "current fitness: " << current_position[2] << ", rank: " << rank << std::endl;

            best_val = current_position[2];
            gen = {current_position[0], current_position[1]};
            current_particle = {gen, current_position[2]};
            current_position[2] = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );

            r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator );
            r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

            vel_x = current_position[3];
            vel_y = current_position[4];

            if(best_val >= current_position[2]){//no change in bestval
                best_x = current_position[5];
                best_y = current_position[6];
            }else{ //new fitness is local max
                best_x = current_position[0];
                best_y = current_position[1];
            }
            vel_x *= m_weight;
            vel_x += m_coef1*r1*( best_x - current_position[3] );
            vel_y *= m_weight;
            vel_y += m_coef1*r1*( best_y - current_position[4] );

            current_position[3] = vel_x;
            current_position[4] = vel_y;
            current_position[5] = r2;

            //std::cout << "Rank: " << rank << ", pre send" << std::endl;
            MPI_Send(current_position, 6, MPI_DOUBLE, 0, curr_particle_number, MPI_COMM_WORLD);
            //std::cout<< "new fitness: " << current_position[2] << ", particle number: " << curr_particle_number << std::endl;
            //double t2 = MPI_Wtime();
            //time = t2 - t1;

        }
        //std::cout << "Rank: " << rank << " ran for: " << time << "seconds\n";
        //delete[] current_position;

    }

   //delete[] particle_to_recieve;


}

}
