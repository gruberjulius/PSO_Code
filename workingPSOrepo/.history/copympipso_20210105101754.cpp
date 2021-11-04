#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include<sys/time.h>
#include"functionstooptimize.hpp"
#include<iostream>
#include<fstream>


#define OUTPUT_FILE "results/times_alternative" 

//#define PI 3.14159265358979323846, done in the other file
double nDimensions, mVelocity, nIterations, seed;
double x_min ;//= -32.768;
double x_max ;//= 32.768;
struct timeval TimeValue_Start;
struct timezone TimeZone_Start;
struct timeval TimeValue_Final;
struct timezone TimeZone_Final;
long time_start, time_end;
double time_overhead;
//double ackley(double x[], double nDimensions);
/*
double ackley(double x[], double nDimensions) {
    double c = 2*M_PI;
    double b = 0.2;
    double a = 20;
    double sum1 = 0;
    double sum2 = 0;
    int i;
    for (i=0; i<nDimensions; i++) {
        sum1 = sum1 + gsl_pow_2(x[i]);
        sum2 = sum2 + cos(c*x[i]);
    }
    double term1 = -a * exp(-b*sqrt(sum1/nDimensions));
    double term2 = -exp(sum2/nDimensions);
    return term1 + term2 + a + M_E;
}
double Griewanks_function(double x[], double nDimensions);
*/


double x;

double y;

double (*function)(double, double);


//double result(double x, double y){
    


    
//    return 0;
//} //, how did we do this?


int main(int argc, char *argv[]) {

    long iterations = strtol(argv[1], NULL, 10);
    int s_case = atoi(argv[2]); //which function to use
    bool basic_file = atoi(argv[3]);//1;
        //if ture results are written to times_papso.dat
        //if not results written to times_papso_#iter_#threads.dat
    int particles = atoi(argv[4]);

    double range[2];


    double (*function)(double, double);


    //double result(double x, double y){
    
    switch(s_case) { // switch case to set correct range of function
        case 0: 
            function =  Griewank_f;
            break;
        case 1: 
            function = Ackley_f;
            break;
        case 2: 
            function =   holdertable;
            break;
        case 3: 
            function =  Rastrigin_f;
            break;
        case 4: 
            function =  Schwefel2_f;
            break;
        case 5:
            function =  VERY_Complicated;
            break;
        default:
            //std::cout << "Only 5 functions used (0,...,4)" << std::endl;
            break;
    }
    
    //std::vector<variable_datatype> range;
    switch(s_case) { // switch case to set correct range of function
        case 0: 
            //range[2] = {-600, 600};
            range[0] = -600;
            range[1] = 600;
            //result(x,y ) = Griewank_f(x,y);
            break;
        case 1: 
            //double range[2] = {-32, 32};
            range[0] = -32;
            range[1] = 32;
            //double result(double x, double y) = Ackley_f(x,y);
            break;
        case 2: 
            //double range[2] = {-10, 10};
            range[0] = -10;
            range[1] = 10;
            //double result = holdertable(x,y);
            break;
        case 3: 
            //double range[2] = {-5.12, 5.12};
            range[0] = -5.12;
            range[1] = 5.12;
            //result = Rastrigin_f(x,y);
            break;
        case 4: 
            //double range[2] = {-100, 100};
            range[0] = -100;
            range[1] = 100;
            //result = Schwefel2_f(x,y);
            break;
        case 5:
            //double range[2] = {-100, 100};
            range[0] = -100;
            range[1] = 100;
            //result = VERY_Complicated(x,y);
            break;
        default:
            //std::cout << "Only 5 functions used (0,...,4)" << std::endl;
            break;
    }
    x_min = range[0];
    x_max = range[1];



    int i,j;
    double nParticles;
    //Argument handling START
    for(i=1; i < argc-1; i++) {
        if (strcmp(argv[i], "-D") == 0)
            nDimensions = strtol(argv[i+1],NULL,10);
        else  if (strcmp(argv[i], "-m") == 0)
            nParticles = strtol(argv[i+1],NULL,10);
        else  if (strcmp(argv[i], "-V") == 0)
            mVelocity = strtol(argv[i+1],NULL,10);
        else  if (strcmp(argv[i], "-i") == 0)
            nIterations = strtol(argv[i+1],NULL,10);
        else  if (strcmp(argv[i], "-s") == 0)
            seed = strtol(argv[i+1],NULL,10);
    }
    if (nDimensions == 0)
        nDimensions = 2;
    if (nParticles == 0)
        nParticles = 1000; //how many particles do we use?
    if (mVelocity == 0)
        mVelocity = 60;
    if (nIterations == 0)
        nIterations = 10000;
    if (seed == 0)
        seed = 1;

  gettimeofday(&TimeValue_Start, &TimeZone_Start);   
int size,myrank,distributed_particles;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&size);
MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
if(myrank==0)
{
    distributed_particles=(int)nParticles/size;
    //printf("%d distributed_particles\n",distributed_particles );
}
MPI_Bcast(&distributed_particles,1,MPI_INT,0,MPI_COMM_WORLD);
if(myrank==0)
{

distributed_particles+=(int)nParticles%size;
//printf("%d distributed_particles\n",distributed_particles );

}
    double result[(int)distributed_particles];
    int step;
    double a,b;
    double c1, c2, rho1, rho2, w, fit;
    c1 = c2 = 1.496;
    w = 0.7298;
    int recievingdata[((int)nDimensions+1)*size];
    int sendingdata[(int)nDimensions+1];
    //Random number generator initialization
    gsl_rng_env_setup();
    gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, time(0));

    double positions[(int)distributed_particles][(int)nDimensions];
    double velocities[(int)distributed_particles][(int)nDimensions];
    double pBestPositions[(int)distributed_particles][(int)nDimensions];    
    double pBestFitness[(int)distributed_particles];
    double gBestPosition[(int)nDimensions];
    double gBestFitness = DBL_MAX;

    //particle initialization
    for (i=0; i<distributed_particles; i++) {
        for (j=0; j<nDimensions; j++) {
            a = x_min + (x_max - x_min) *  gsl_rng_uniform(r);
            b = x_min + (x_max - x_min) *  gsl_rng_uniform(r);
            positions[i][j] = a;
            pBestPositions[i][j] = a;
            velocities[i][j] = (a-b) / 2.;
        }
        double x = positions[i][0];
        double y = positions[i][1];
        pBestFitness[i] = ackley(x,y);
        if (pBestFitness[i] < gBestFitness) {
            memmove((void *)gBestPosition, (void *)&positions[i], sizeof(double) * nDimensions);
            gBestFitness = pBestFitness[i];
        } 
    }

    //actual calculation
    for (step=0; step<nIterations; step++) {
        for (i=0; i<distributed_particles; i++) {
            for (j=0; j<nDimensions; j++) {
                // calculate stochastic coefficients
                rho1 = c1 * gsl_rng_uniform(r);
                rho2 = c2 * gsl_rng_uniform(r);
                // update velocity
                velocities[i][j] = w * velocities[i][j] + \
                rho1 * (pBestPositions[i][j] - positions[i][j]) +  \
                rho2 * (gBestPosition[j] - positions[i][j]);
                // update position
                positions[i][j] += velocities[i][j];

                if (positions[i][j] < x_min) {
                    positions[i][j] = x_min;
                    velocities[i][j] = 0;
                } else if (positions[i][j] > x_max) {
                    positions[i][j] = x_max;
                    velocities[i][j] = 0;
                }

            }
            double x = positions[i][0];
            double y = positions[i][1];

            // update particle fitness
            fit = ackley(x,y);
            // update personal best position?
            if (fit < pBestFitness[i]) {
                pBestFitness[i] = fit;
                // copy contents of positions[i] to pos_b[i]
                memmove((void *)&pBestPositions[i], (void *)&positions[i],
                    sizeof(double) * nDimensions);
            }
            // update gbest??
            if (fit < gBestFitness) {
                // update best fitness
                gBestFitness = fit;
                // copy particle pos to gbest vector
                memmove((void *)gBestPosition, (void *)&positions[i],
                    sizeof(double) * nDimensions);
            }
        }
        for(int k=0;k<(nDimensions);k++)
            sendingdata[k]=gBestPosition[k];
        sendingdata[(int)nDimensions]=gBestFitness;
        MPI_Gather(&sendingdata,nDimensions+1, MPI_INT,&recievingdata,nDimensions+1, MPI_INT, 0, MPI_COMM_WORLD);
        if(myrank==0)
        {
            int min=gBestFitness;
            int pos=-1;
            for(int k=0;k<size;k++)
            { //printf("%d\n",recievingdata[k*((int)nDimensions+1)+((int)nDimensions)] );
                if(min>=recievingdata[k*((int)nDimensions+1)+((int)nDimensions)])
                    {
                        min=recievingdata[k*((int)nDimensions+1)+((int)nDimensions)];
                        pos=k*((int)nDimensions+1);
                    }   
            }
            gBestFitness=min;
            for(int k=pos;k<nDimensions+pos;k++)
                gBestPosition[k-pos]=recievingdata[k];                  
        }
        MPI_Bcast(&gBestPosition,nDimensions,MPI_INT,0,MPI_COMM_WORLD);
       // MPI_Bcast(&gBestFitness,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    if(myrank==0)
    {
        printf("Result: %f\n", gBestFitness);
                gettimeofday(&TimeValue_Final, &TimeZone_Final);
        time_start = TimeValue_Start.tv_sec * 1000000 + TimeValue_Start.tv_usec;
        time_end = TimeValue_Final.tv_sec * 1000000 + TimeValue_Final.tv_usec;
        time_overhead = (time_end - time_start)/1000000.0;
        printf("\n Time in Seconds (T) : %lf\n",time_overhead);
    }   
    gsl_rng_free(r);
    MPI_Finalize();

    //begin modification
    int numproces = 1; // need to update this, we are using more than one processor
    std::string w_file_name = OUTPUT_FILE; //make this C ready
    if(!basic_file){
        std::cout<<"Results written to long filename!" << std::endl;
        w_file_name += "_";
        w_file_name += std::to_string(iterations); //make this C ready
        w_file_name += "_";
        w_file_name += std::to_string(particles); //make this C ready
        //w_file_name += "_m";
    }
    w_file_name += ".dat";

    std::ofstream times(w_file_name, std::ios::app); //make this C ready
    times << s_case << " " << numproces << " " << time_overhead << " " << gBestFitness << "\n";
    times.close();


    //end modification
}
 //mpicc mpipso.c -lm -lgsl -lgslcblas -o mpi