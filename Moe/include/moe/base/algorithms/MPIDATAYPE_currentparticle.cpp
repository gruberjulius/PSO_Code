
#pragma once

#include "NumericAlgorithmImpl.hpp"
#include "../parameters/PSParameters.hpp"
#include <cassert>
#include <queue>
#include <mpi.h>
#include <iostream>

struct particle{
  float genotype11;
  float genotype12;
  double fitness;
  float velocities1;
  float velocities2;
  float genotype21;
  float genotype22;
  int idx;
};

particle particles_to_scatter[numproces-1];
particle particle_to_recieve;

particles_to_scatter[j].genotype11 = m_population[j-1].genotype[0];//x
particles_to_scatter[j].genotype12 = m_population[j-1].genotype[1];//y
particles_to_scatter[j].fitness = m_best_genotypes[j-1].fitness;//best val
particles_to_scatter[j].velocities1 = m_velocities[j-1][0];//u
particles_to_scatter[j].velocities2 = m_velocities[j-1][1];//v
particles_to_scatter[j].genotype21 = m_best_genotypes[j-1].genotype[0];//best_x
particles_to_scatter[j].genotype22 = m_best_genotypes[j-1].genotype[1];//best_y
particles_to_scatter[j].idx = j-1;//particle index


MPI_Datatype ScatterType;
MPI_Type_create_struct(
  4,
  [2,1,4,1],
  [0,8,16,32],
  [MPI_FLOAT,MPI_DOUBLE,MPI_FLOAT,MPI_INT],
  &ScatterType
);
MPI_Type_commit(&ScatterType);

MPI_Scatter(particles_to_scatter, 1, ScatterType, particle_to_recieve, 1, ScatterType, 0, MPI_COMM_WORLD);


//TYPE 2 smaller
struct particle{
  float genotype11;
  float genotype12;
  int idx;
};

particle particles_to_scatter[numproces-1];
particle particle_to_recieve;

particles_to_scatter[j].genotype11 = m_population[j-1].genotype[0];//x
particles_to_scatter[j].genotype12 = m_population[j-1].genotype[1];//y
particles_to_scatter[j].idx = j-1;//particle index


MPI_Datatype ScatterType;
MPI_Type_create_struct(
  2,
  [2,1],
  [0,8],
  [MPI_FLOAT,MPI_INT],
  &ScatterType
);
MPI_Type_commit(&ScatterType);

MPI_Scatter(particles_to_scatter, 1, ScatterType, particle_to_recieve, 1, ScatterType, 0, MPI_COMM_WORLD);



//STUFF FOR KRISTOF IMPL. //////////////////////////////////////////////
/*#include <stddef.h>  // or <cstddef> for C++

struct particle_d
{
   double x;
   double y;
   double fit;
};

MPI_Datatype createRecType()
{
    // Set-up the arguments for the type constructor
    MPI_Datatype new_type;

    int count = 3;
    int blocklens[] = { 1,1,1 };

    MPI_Aint indices[3];
    indices[0] = (MPI_Aint)offsetof(struct particle_d, x);
    indices[1] = (MPI_Aint)offsetof(struct particle_d, y);
    indices[2] = (MPI_Aint)offsetof(struct particle_d, fit);

    MPI_Datatype old_types[] = {MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_struct(count,blocklens,indices,old_types,&new_type);
    MPI_Type_commit(&new_type);

    return new_type;
}

double mpi_particle_max(){

}

MPI_Op MPI_PARTICLE_MAX;
int MPI_PARTICLE_MAX_create MPI_Op_create((MPI_User_function *)mpi_particle_max, 1, &MPI_PARTICLE_MAX);

/*
int MPIAPI MPI_Op_create(
  _In_  MPI_User_function *function,
        int               commute,
  _Out_ MPI_Op            *op
);
*/