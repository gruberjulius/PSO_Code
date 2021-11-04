// STOPPED IMPLEMENTING

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
class ParticleSwarmPA_prec_vel_2 : public NumericAlgorithm<GenotypeType>
{
    public:
        ParticleSwarmPA_prec_vel_2( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions = 1, std::vector<GenotypeType> _range = { std::numeric_limits<GenotypeType>::lowest() , std::numeric_limits<GenotypeType>::max() });
        ParticleSwarmPA_prec_vel_2( const PSParameters<GenotypeType>& _parameters );

        void run (unsigned int _iterations, unsigned int dim);
        void run (unsigned int _iterations);

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
ParticleSwarmPA_prec_vel_2<GenotypeType>::ParticleSwarmPA_prec_vel_2( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions, std::vector<GenotypeType> _range )
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
ParticleSwarmPA_prec_vel_2<GenotypeType>::ParticleSwarmPA_prec_vel_2( const PSParameters<GenotypeType>& _parameters )
:ParticleSwarmPA_prec_vel_2<GenotypeType>( _parameters.moesPerGen, _parameters.inertia, _parameters.coef1, _parameters.coef2, _parameters.dimensions, _parameters.range )
{
    
}

template <typename GenotypeType>
void ParticleSwarmPA_prec_vel_2<GenotypeType>::init( unsigned int _iterations )
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
void ParticleSwarmPA_prec_vel_2<GenotypeType>::run( unsigned int _iterations, unsigned int dim)
{
    // dim should == NumericAlgorithm<GenotypeType>::m_dimensions
    int rank, numproces;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);

    /*double *particles_to_scatter = new double [(numproces-1)*4];
    double *particle_to_recieve = new double [4];*/
    //int scattersize = (numproces-1)*6;

    //dim is the # of dimensions
    unsigned int dim2 = 2 * dim;
    unsigned int dim3 = 3 * dim;
    unsigned int vec_size = dim3 + 2;
    double particles_to_scatter[(numproces-1)*vec_size];
    //position, velocity, best_position, fitness, particle#
    double particle_to_recieve[vec_size]; 
    //WARNING: CHANGED ORDER COMPARED TO PREVIOUS IMPL.
    if(rank == 0){
        this->init( _iterations );
        //initialize first "entry" of particles_to_scatter as dummy
        //dummy will be recieved (scatter) by master, where we don't need data
        for(unsigned int i = 0; i < vec_size; ++i){
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
        for(unsigned int j = 1; j<(numproces); ++j){
            for(unsigned int k = 0; k < dim; ++k){
                particles_to_scatter[j*vec_size + k] = m_population[j-1].genotype[0];//x
                particles_to_scatter[j*vec_size + dim + k] = m_velocities[j-1][0];//u
                particles_to_scatter[j*vec_size + dim2 + k] = m_best_genotypes[j-1].genotype[0];//best_x
            }
            particles_to_scatter[vec_size - 2] = m_best_genotypes[j-1].fitness;//best val
            particles_to_scatter[vec_size - 1] = j-1;//particle index
        }
//        std::cout << "Pre scatter\n";

        /*for ( int i = 0; i < 4*(numproces); ++i){
            std::cout << particles_to_scatter[i] << " ";
        }
        std::cout<<"this was master ^ " << std::endl;//*/
    }
    //MPI_ScatterV MAYBE
    MPI_Scatter(particles_to_scatter, vec_size, MPI_DOUBLE, particle_to_recieve, vec_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //delete[] particles_to_scatter;
        
        int send_size = vec_size - 1;
        int recv_size = dim2;

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

//WARNING: CHANGED ORDER COMPARED TO PREVIOUS IMPL.
        double our_send_particle[recv_size];
        //[x, y, u, v, local_best_x, local_best_y, local_best_fitness]
        double our_recv_particle[recv_size];
        //[local_best_x, local_best_y, u, v, fitness, r2]

        //std::cout<<"# Iterations: "<<m_iterations<<std::endl;
        
        int fit_idx = dim3 + 1;

        double modified_iterations = num_particles * m_iterations;
        while(curr_sum_iters < modified_iterations)//was step < m_iterations
        {
//                std::cout<<"Checkpoint 2: while loop step: "<< step <<std::endl;

                MPI_Status status;
                MPI_Recv(our_recv_particle, recv_size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//need MPI_Datatype of Particle

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
                for(unsigned int i = 0; i < dim; ++i){
                    our_send_particle[i] = current_particle.genotype[i];
                }
                //current velocities
                for(unsigned int i = dim; i < dim2; ++i){
                    our_send_particle[i] = m_velocities[new_idx][i];
                }
                //best local x and y
                for(unsigned int i = dim2; i < dim3; ++i){
                    our_send_particle[dim2 + i] = m_best_genotypes[current_idx].genotype[dim2 + i];
                }
                //best local fitness
                our_send_particle[fit_idx] = m_best_genotypes[current_idx].fitness;

                //std::cout<<"Master: Pre send" << std::endl;
                double tag = new_idx + 1.;
                MPI_Send(our_send_particle, send_size, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD);
                //std::cout<<"Master: Post send" << std::endl;

                //update fitness and velocities
                m_population[current_idx].fitness = our_recv_particle[fit_idx];

                m_velocities[current_idx][0] = our_recv_particle[3];
                m_velocities[current_idx][1] = our_recv_particle[4];

                //checks if function evaluation from slave is new best
                if(our_recv_particle[fit_idx] > m_best_genotypes[current_idx].fitness){
                //update m_best_genotypes(local best) with current particle
                    m_best_genotypes[current_idx].genotype[0] = our_recv_particle[0];
                    m_best_genotypes[current_idx].genotype[1] = our_recv_particle[1];
                    m_best_genotypes[current_idx].fitness = our_recv_particle[fit_idx];
                    //checks if function evaluation from slave is new global best
                    if( m_best_genotypes[current_idx].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness )
                        Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[current_idx];
                }

                int j = current_idx;//finish update of previous particle
                for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )// 2D -> k = 0, 1; UNROLL ?
                {

                    float   r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

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
                //std::cout<<"Step: "<<step<< ", curr_sum_iter: " << curr_sum_iters << ", num_particles: " << num_particles << ", # iter: " << m_iterations<<std::endl;
        }
        //std::cout<<"Checkpoint 9: while loop exited"<<std::endl;
        
        //loop to kill all slave processes
        for(int id = 1; id < numproces; ++id){
            MPI_Status status;
            MPI_Recv(our_recv_particle, recv_size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int particle_number = status.MPI_TAG;
            if(our_recv_particle[2] > m_best_genotypes[particle_number].fitness){//OPTIMIZE WITH OR WITHOUT THIS
            //update m_best_genotypes(local best) with current particle
                m_best_genotypes[particle_number].genotype[0] = our_recv_particle[0];
                m_best_genotypes[particle_number].genotype[1] = our_recv_particle[1];
                m_best_genotypes[particle_number].fitness = our_recv_particle[2];
                //checks if function evaluation from slave is new global best
                if( m_best_genotypes[particle_number].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness )
                    Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[particle_number];
            }
            //std::cout<<"terminator send to: "<<status.MPI_SOURCE<<std::endl;
            MPI_Send(our_recv_particle, send_size, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        }

        //std::cout<<"Checkpoint 10: terminator for-loop exited"<<std::endl;

    }


    if(rank != 0){ //slave process

        //double *current_position = new double[3];
        double current_position[vec_size];
        //double time;
        //double t1 = MPI_Wtime();

        /*for ( int i = 0; i < 4; ++i){
            std::cout << particle_to_recieve[i] << " ";
        }
        std::cout<<std::endl;//*/
        //std::cout<<"m_weight" << m_coef1<<std::endl;
        
        std::copy(particle_to_recieve, particle_to_recieve + vec_size, current_position);
        int curr_particle_number = (int)(double)particle_to_recieve[send_size];

        int fit_idx = dim3 + 1; //index of fitness value in current_particle vector
        double best_val = current_position[fit_idx];//save best value in temp

        
        std::vector<double> gen;
        //std::copy(particle_to_recieve, particle_to_recieve + dim, gen);
        for(unsigned int i = 0; i < dim; ++i){
            gen.push_back(particle_to_recieve[i]);
        }


        Moe<double> current_particle = {gen, 42};//particle value not important for func_eval
        current_position[fit_idx] = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );

        //temp variables for velocities and best positions
        double vel[dim];
        double best_pos[dim];
        for(unsigned int i = 0; i < dim; ++i){
            vel[i] = current_position[dim + i];
        }

        
        if(best_val >= current_position[fit_idx]){//no change in bestval
            for(unsigned int i = 0; i < dim; ++i){
                best_pos[i] = current_position[dim2 + i];
            }
        }else{ //new fitness is local max
            for(unsigned int i = 0; i < dim; ++i){
                best_pos[i] = current_position[i];
            }
        }
        float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

        for( unsigned int k = 0; k < dim; k++ )// 2D -> k = 0, 1; UNROLL ?
        {
            vel[k] *= m_weight;
            vel[k] += m_coef1*r1*( best_pos[k] - current_position[k] ); //current_position[k] is the actual position of k_th dim
        }

        std::copy(vel, vel + dim, current_position + dim);

        //std::cout << "Rank: " << rank << " pre-first send\n";
        MPI_Send(current_position, recv_size, MPI_DOUBLE, 0, curr_particle_number, MPI_COMM_WORLD);
        //std::cout << "Rank: " << rank << " first send done\n";
        while(1){

            MPI_Status status;
            MPI_Recv(current_position, send_size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//need MPI_Datatype of Particle
            //std::cout<<"Recieved tag: " << status.MPI_TAG << std::endl;
            if(status.MPI_TAG == 0){
                //std::cout<<"BREAKED rank: " <<rank<<std::endl;
                break;
            }
            
            curr_particle_number = status.MPI_TAG - 1;//decr to make tag into actual particle index

            //std::cout<< "current fitness: " << current_position[2] << ", rank: " << rank << std::endl;

            best_val = current_position[fit_idx];//save best value in temp

            //std::vector<double> gen;
            //std::copy(particle_to_recieve, particle_to_recieve + dim, gen);
            for(unsigned int i = 0; i < dim; ++i){
                gen.push_back(particle_to_recieve[i]);
            }

            current_particle = {gen, 42};//particle value not important for func_eval
            current_position[fit_idx] = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );

            //temp variables for velocities and best positions

            for(unsigned int i = 0; i < dim; ++i){
                vel[i] = current_position[dim + i];
            }

            
            if(best_val >= current_position[fit_idx]){//no change in bestval
                for(unsigned int i = 0; i < dim; ++i){
                    best_pos[i] = current_position[dim2 + i];
                }
            }else{ //new fitness is local max
                for(unsigned int i = 0; i < dim; ++i){
                    best_pos[i] = current_position[i];
                }
            }

            r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

            for( unsigned int k = 0; k < dim; k++ )// 2D -> k = 0, 1; UNROLL ?
            {
                vel[k] *= m_weight;
                vel[k] += m_coef1*r1*( best_pos[k] - current_position[k] ); //current_position[k] is the actual position of k_th dim
            }

            std::copy(vel, vel + dim, current_position + dim);
            
            //std::cout << "Rank: " << rank << ", pre send" << std::endl;
            MPI_Send(current_position, recv_size, MPI_DOUBLE, 0, curr_particle_number, MPI_COMM_WORLD);
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
