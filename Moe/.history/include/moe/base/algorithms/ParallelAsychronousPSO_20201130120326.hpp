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
class ParticleSwarmPA : public NumericAlgorithm<GenotypeType>
{
    public:
        ParticleSwarmPA( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions = 1, std::vector<GenotypeType> _range = { std::numeric_limits<GenotypeType>::lowest() , std::numeric_limits<GenotypeType>::max() });
        ParticleSwarmPA( const PSParameters<GenotypeType>& _parameters );

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
ParticleSwarmPA<GenotypeType>::ParticleSwarmPA( unsigned int _moesPerGen, float _weight, float _coef1, float _coef2, unsigned int _dimensions, std::vector<GenotypeType> _range )
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
ParticleSwarmPA<GenotypeType>::ParticleSwarmPA( const PSParameters<GenotypeType>& _parameters )
:ParticleSwarmPA<GenotypeType>( _parameters.moesPerGen, _parameters.inertia, _parameters.coef1, _parameters.coef2, _parameters.dimensions, _parameters.range )
{
    
}

template <typename GenotypeType>
void ParticleSwarmPA<GenotypeType>::init( unsigned int _iterations )
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
void ParticleSwarmPA<GenotypeType>::run( unsigned int _iterations){std::cout<<"The wrong \" run \" function has been called in ParticleSwarmPA.cpp"<<std::endl;}
template <typename GenotypeType>
void ParticleSwarmPA<GenotypeType>::run( unsigned int _iterations, int& argc, char * argv[])
{
    
//    MPI_Init(&argc, &argv); 

    int rank, numproces;
    //char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);

    //if(rank == 0){ //master process
        //MPI_Send(buf, 2, MPI_DOUBLE, next_available, 1, MPI_COMM_WORLD);//need MPI_Datatype of Particle*            
//        std::cout << "My rank is: " << rank << std::endl;
        //MPI_Recv(buf, 1, MPI_DOUBLE, next_available, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //}
    double *particles_to_scatter = new double [(numproces-1)*4];
    double *particle_to_recieve = new double [4];
    if(rank == 0){
        this->init( _iterations );
        for(int j = 0; j<numproces-1; ++j){
            particles_to_scatter[j*4] = m_population[j].genotype[0];
            particles_to_scatter[j*4 + 1] = m_population[j].genotype[1];
            particles_to_scatter[j*4 + 2] = m_population[j].fitness;
            particles_to_scatter[j*4 + 3] = j;
        }
//        std::cout << "Pre scatter\n";

    }
    MPI_Scatter(particles_to_scatter, 4, MPI_DOUBLE, particle_to_recieve, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(rank == 0){
        //initialization
//        std::cout << "Post scatter\n";

        int num_particles = m_population.size();
//        std::vector<int> particle_idx_vec;
//        particle_idx_vec.reserve(num_particles);
//        std::cout<<"Checkpoint 0 (rank 0): pre-for particle_idx_vec"<<std::endl;

/*        for(int i = 0; i < num_particles; ++i){
            particle_idx_vec[i] = i;
        }*/

        /*double all_positions[num_particles * 2];
        for(int i = 0; i < num_particles; ++i){
            all_positions[i*2] = m_population[i].genotype[0];
            all_positions[i*2+1] = m_population[i].genotype[1];
        }*///maybe we should just send particle/positions. This way we are not using the Population class implemented by the code
/*        bool running_threads[numproces-1];
        if(numproces < numparticles){
            std::fill(running_threads, running_threads + numproces-1,1);
        }else{
            std::fil_n(running_threads,running_threads + numproces-1,1);
            for(int i = numproces-1; i < numparticles; ++i){
                running_threads[i] = 0;
            }
        }*/
        std::queue<int> particlequeue;//, available_threads; //queue-s of particle indexes and available thread ranks
        //initialization of particle queue with numproces-1,...,#particles
        for(int i = numproces-1; i < num_particles; ++i){ //numproces-1 because we already scattered
            particlequeue.push(i);
        }
        /*
        for(int i = 1; i < numproces; ++i){ //initialization of available threads queue with 0,...,numproces
            available_threads.push(i);
        }*/
        unsigned int step = 0;
        int curr_sum_iters = 0; //current # of particles that have been evaluated (used for step, to have a rough bound on limitations)
//        std::cout<<"Checkpoint 1: pre-while loop"<<std::endl;

        double our_current_particle[3];
        
        while(step < m_iterations)
        {
//                std::cout<<"Checkpoint 2: while loop step: "<< step <<std::endl;

                MPI_Status status;
                MPI_Recv(&our_current_particle, 3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//need MPI_Datatype of Particle

                int current_idx = status.MPI_TAG - 1;
                m_population[current_idx].fitness = our_current_particle[2];
                
                particlequeue.push(current_idx);
                //available_threads.push(status.MPI_SOURCE);

                current_idx = particlequeue.front(); //get new particle idx
                particlequeue.pop();//deletes accessed particle
//                std::cout<<"Checkpoint 3: this sh!t is poppin"<<std::endl;
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
//                std::cout<<"Checkpoint 4: post-4 loop"<<std::endl;
                Moe<GenotypeType> current_particle = m_population[current_idx];
                our_current_particle[0] = current_particle.genotype[0];
                our_current_particle[1] = current_particle.genotype[1];
                our_current_particle[2] = current_particle.fitness;
                MPI_Send(&our_current_particle, 3, MPI_DOUBLE, status.MPI_SOURCE, current_idx + 1, MPI_COMM_WORLD);


                if( m_population[current_idx].fitness > m_best_genotypes[current_idx].fitness )
                {
                    m_best_genotypes[current_idx] = m_population[current_idx];
                    if( m_best_genotypes[current_idx].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness )
                        Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[current_idx];
                }
//                std::cout<<"Checkpoint 8: full while loop iteration"<<std::endl;

                curr_sum_iters += 1;
                step = curr_sum_iters / num_particles;
        }
//        std::cout<<"Checkpoint 9: while loop exited"<<std::endl;
        
        //loop to kill all slave processes
        for(int id = 1; id < numproces; ++id){
            MPI_Status status;
            MPI_Recv(&our_current_particle, 3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if( m_population[status.MPI_TAG-1].fitness > m_best_genotypes[status.MPI_TAG-1].fitness )//OPTIMISE WITH OR WITHOUT CHECK CONVERGENCE
                {
                    m_best_genotypes[status.MPI_TAG-1] = m_population[status.MPI_TAG-1];
                    if( m_best_genotypes[status.MPI_TAG-1].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness )
                        Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[status.MPI_TAG-1];
                }
            MPI_Send(&our_current_particle, 3, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        }

/*
        int finished = 1; //TRUE
        MPI_Bcast(&finished,1,MPI_INT,0,MPI_COMM_WORLD);*/
        
    }
/*double current_position[];
    }
    if(rank != 0){
        double current_position = malloc(2 * sizeof(double));
        double all_positions[];
    }
    MPI_Scatter(&all_positions, 2 * num_particles, MPI_DOUBLE, &current_position, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD); //may have problem, because current_position doesn't exist in rank=0 process
    if(rank == 0){
    }//*/

    if(rank != 0){ //slave process

        double current_position[3];
        double time;

        std::copy(particle_to_recieve, particle_to_recieve + 2, current_position);
        int curr_particle_number = particle_to_recieve[3];

        std::vector<double> gen = {current_position[0], current_position[1]};

        Moe<double> current_particle = {gen, current_position[2]};
        current_position[2] = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );

        std::cout << "Rank: " << rank << " pre-first send\n";
        MPI_Send(&current_position, 3, MPI_DOUBLE, 0, curr_particle_number + 1, MPI_COMM_WORLD);
        std::cout << "Rank: " << rank << " first send done\n";
        while(1){
            /*
            int finished; //TRUE
            MPI_Bcast(&finished,1,MPI_INT,0,MPI_COMM_WORLD);
            if(finished){
                break;
            }*/

            MPI_Status status;
            MPI_Recv(&current_position, 3, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//need MPI_Datatype of Particle

            if(status.MPI_TAG == 0){
                break;
            }
            //std::cout<< "x: " << current_position[0] << ", y: " << current_position[1] << std::endl;
            double t1 = MPI_Wtime();
            gen = {current_position[0], current_position[1]};
            current_particle = {gen, current_position[2]};
            current_position[2] = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );


            MPI_Send(&current_position, 3, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
            double t2 = MPI_Wtime();
            time = t2 - t1;
            
        }
        std::cout << "Rank: " << rank << " ran for: " << time << "seconds\n";

/*
        Moe<GenotypeType> current_particle;

        MPI_Recv(&current_particle, 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//need MPI_Datatype of Particle

        current_particle.fitness = Algorithm<GenotypeType>::m_fitnessFunction( current_particle );

        MPI_Send(&current_particle, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
*/
    }
///////////////////////
/*
    for( unsigned int i = 0; i < m_iterations; i++ )
    {
        for( unsigned int j = 0; j < m_population.size(); j++ )
        {
            for( unsigned int k = 0; k < NumericAlgorithm<GenotypeType>::m_dimensions; k++ )// 2D -> k = 0, 1;
            {
                float   r1 = m_dist_coef( Algorithm<GenotypeType>::m_generator ),
                        r2 = m_dist_coef( Algorithm<GenotypeType>::m_generator );

                m_velocities[j][k] *= m_weight;
                m_velocities[j][k] += m_coef1*r1*( m_best_genotypes[j].genotype[k] - m_population[j].genotype[k] );
                m_velocities[j][k] += m_coef2*r2*( Algorithm<GenotypeType>::m_bestMoe.genotype[k] - m_population[j].genotype[k] );
                
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
                if( m_best_genotypes[j].fitness > Algorithm<GenotypeType>::m_bestMoe.fitness )
                    Algorithm<GenotypeType>::m_bestMoe = m_best_genotypes[j];
            }
        }
    }
    //*/
//    MPI_Finalize();

}

}
