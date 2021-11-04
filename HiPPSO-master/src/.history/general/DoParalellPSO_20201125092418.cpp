#include "general/main.h"

#include <algorithm>
#include <dirent.h>
#include <arbitrary_precision_calculation/random_number_generator.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

#include "function/function.h"
#include "general/check_condition.h"
#include "general/configuration.h"
#include "general/general_objects.h"
#include "arbitrary_precision_calculation/operations.h"
#include "arbitrary_precision_calculation/configuration.h"
#include "general/particle.h"
#include "general/visualization.h"
#include "neighborhood/neighborhood.h"
#include "position_and_velocity_updater/position_and_velocity_updater.h"
#include "statistics/direct_statistics.h"
#include "statistics/statistics.h"

#include "DoParalellPSO.h"
#include "main.h"



namespace highprecisionpso {

Statistics* DoParalellPso(Statistics* statistics, time_t& LAST_BACKUP, time_t& LAST_RUN_CHECK, int& argc, char * argv[]) {
    //////////////////////////////////EDIT START, MPI initialization
    MPI_Init(&argc, &argv); 

    int rank, numproces;
    //char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);

    if(rank == 0){ //master process
    //////////////////////////////EDIT END
        std::vector<Particle*>* &swarm = statistics->swarm;
        long long & step = statistics->current_iteration;

        int lastNumberOfmpft = arbitraryprecisioncalculation::mpftoperations::GetNumberOfMpftValuesInUse() - arbitraryprecisioncalculation::mpftoperations::GetNumberOfMpftValuesCached();
        int startstep = step;

        LAST_BACKUP = time(NULL);
        LAST_RUN_CHECK = 0;

        sort(configuration::g_preserve_backup_times.begin(),
                configuration::g_preserve_backup_times.end());
        unsigned int nextBackupStepId = 0;
        while(nextBackupStepId < configuration::g_preserve_backup_times.size() 
                && step > configuration::g_preserve_backup_times[nextBackupStepId]){
            ++nextBackupStepId;
        }
        arbitraryprecisioncalculation::Configuration::ResetIncreasePrecisionRecommended();
        if(configuration::g_debug_swarm_activated && step == 0){
            visualization::VisualizeCurrentSwarm();
        }
        /////////////////////////////////////////////////EDIT START
        std::queue<int> particlequeue, available_threads; //queue-s of particle indexes and available thread ranks
        for(int i = 0; i < configuration::g_particles; ++i){ //initialization of particle queue with 0,...,#particles
            particlequeue.push(i);
        }
        for(int i = 1; i < numproces; ++i){ //initialization of available threads queue with 0,...,numproces
            available_threads.push(i);
        }        
        int curr_sum_iters = 0; //current # of particles that have been evaluated (used for step, to have a rough bound on limitations)
        while(step < configuration::g_max_iterations){//while iterations
            void * buf; //needed for the MPI + mpf_t
            int current_idx = particlequeue.front(); //get new particle idx
            particlequeue.pop();//deletes accessed particle

            Particle* current_particle = (*swarm)[current_idx]; //current particle declaration
            current_particle->UpdatePosition(); //update position and velocity

            while(available_threads.empty()){}//to be sure we have a thread to work on
            int next_available = available_threads.front(); //next available thread, from queue for available threads
            available_threads.pop();//deletes accessed thread

            std::vector<double> vel = arbitraryprecisioncalculation::vectoroperations::MpftToDouble(current_particle->GetPosition());
            double* position_mpi = vel.data();
            
            //transform to non blocking, buffers, #makethiscodegreatagain
            mpf_set_ui( * position_mpi); 
            buf = allocbuf_mpf(mpf_get_prec(*position_mpi),1);
            pack_mpf(*position_mpi, 1, buf);
            MPI_Send(buf, 2, MPI_DOUBLE, next_available, 1, MPI_COMM_WORLD);//need MPI_Datatype of Particle*            
            double cost_value; //make into array to have threads buffer to avoid deadlock
        //  MPI_Send(cost_value, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD); the stuff that is sended

            buf = allocbuf_mpf(mpf_get_prec(&cost_function), 1);
            //MPI_Recv(buf, 1, MPI_MPF, 0, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(buf, 1, MPI_DOUBLE, next_available, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //passed by reference as recommended
            unpack_mpf(buf, &cost_function, 1);
            //from belowMPI_Recv(current_particle, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            particlequeue.push(current_idx);//putting particle back in loop

            time_t currentTime = time(NULL);//check convergence
            if (difftime(currentTime, LAST_RUN_CHECK)
                    > configuration::g_time_between_run_checks) {
                if (!AllowedToRun()) {
                    Shutdown();
                    return statistics;
                }
                LAST_RUN_CHECK = currentTime;
            }
            curr_sum_iters++;
            step = curr_sum_iters / configuration::g_particles;

        }

/////////////////////////////////////////////////EDIT END
/**        if(step == 0){
            statistics->EvaluateStatistics();
        }
        while (step < configuration::g_max_iterations) {

            time_t currentTime = time(NULL);
            if (difftime(currentTime, LAST_RUN_CHECK)
                    > configuration::g_time_between_run_checks) {
                if (!AllowedToRun()) {
                    Shutdown();
                    return statistics;
                }
                LAST_RUN_CHECK = currentTime;
            }
            if (difftime(currentTime, LAST_BACKUP)
                    > configuration::g_time_between_backups) {
                std::string tmpFilename = configuration::g_file_prefix + ".backup";
                WriteCurrentState(tmpFilename);
            }

            if (nextBackupStepId < configuration::g_preserve_backup_times.size()
                    && configuration::g_preserve_backup_times[nextBackupStepId] == step){
                std::string tmpFilename = configuration::g_file_prefix + ".backup";
                WriteCurrentState(tmpFilename);
                std::ostringstream os;
                os << configuration::g_file_prefix << ".S" << step << ".backup";
                tmpFilename = os.str();
                WriteCurrentState(tmpFilename);

                while(nextBackupStepId < configuration::g_preserve_backup_times.size()
                        && step >= configuration::g_preserve_backup_times[nextBackupStepId]){
                    ++nextBackupStepId;
                }
            }


            if (lastNumberOfmpft != arbitraryprecisioncalculation::mpftoperations::GetNumberOfMpftValuesInUse() - arbitraryprecisioncalculation::mpftoperations::GetNumberOfMpftValuesCached()) {
                if (step > 2 + startstep) {
                    std::cerr << "step " << step << std::endl << "last mpf_t "
                        << lastNumberOfmpft << std::endl << "curr mpf_t "
                        << arbitraryprecisioncalculation::mpftoperations::GetNumberOfMpftValuesInUse() - arbitraryprecisioncalculation::mpftoperations::GetNumberOfMpftValuesCached() << std::endl;
                }
                lastNumberOfmpft = arbitraryprecisioncalculation::mpftoperations::GetNumberOfMpftValuesInUse() - arbitraryprecisioncalculation::mpftoperations::GetNumberOfMpftValuesCached();
            }
            for (int id = 0; id < configuration::g_particles; id++) {
                (*swarm)[id]->UpdatePosition();
                if(configuration::g_update_global_attractor_mode == configuration::UPDATE_GLOBAL_ATTRACTOR_MODE_EACH_PARTICLE){
                    configuration::g_neighborhood->ProceedAllUpdates();
                }
                if(arbitraryprecisioncalculation::Configuration::isIncreasePrecisionRecommended()){
                    arbitraryprecisioncalculation::Configuration::ResetIncreasePrecisionRecommended();
                    arbitraryprecisioncalculation::mpftoperations::IncreasePrecision();
                }
            }
            if(configuration::g_update_global_attractor_mode == configuration::UPDATE_GLOBAL_ATTRACTOR_MODE_EACH_ITERATION){
                configuration::g_neighborhood->ProceedAllUpdates();
            }

            step++;

            if(configuration::g_debug_swarm_activated){
                visualization::VisualizeCurrentSwarm();
            }

            if(arbitraryprecisioncalculation::Configuration::isIncreasePrecisionRecommended()){
                arbitraryprecisioncalculation::Configuration::ResetIncreasePrecisionRecommended();
                arbitraryprecisioncalculation::mpftoperations::IncreasePrecision();
            }

            statistics->EvaluateStatistics();

        }

        {
            std::string tmpFilename = configuration::g_file_prefix + ".backup";
            WriteCurrentState(tmpFilename);
        }
        FILE* logging = fopen((configuration::g_file_prefix + ".log").c_str(), "a");
        time_t rawtime = time(NULL);
        struct tm * timeinfo = localtime(&rawtime);
        std::string tmpstr(asctime(timeinfo));
        if (tmpstr[tmpstr.size() - 1] == '\n')
            tmpstr = tmpstr.substr(0, tmpstr.size() - 1);
        fprintf(logging, "finished %s with %lld steps\n", tmpstr.c_str(),
                statistics->current_iteration);
        fclose(logging);
    }
**/
////////////////////////////////
    if(rank != 0){ //slave process
        void * buf; //needed for the mpf_t + MPI

        buf = allocbuf_mpf(mpf_get_prec(&current_particle), 1);


        mpf_t* cost_value;
        Particle* current_particle;
        MPI_Recv(&current_particle, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//need MPI_Datatype of Particle*
        unpack_mpf(buf, &current_particle, 1);

        std::vector<mpf_t*> currPos = current_particle->GetPosition();
        GeometricAverageReduceOperation R;
        cost_value = R.GeometricAverageReduceOperation::Evaluate(currPos);  //should be passed by reference
        mpf_set_ui(&cost_value, 1);
        buf = allocbuf_mpf(mpf_get_prec(cost_value), 1);
        pack_mpf(cost_value, 1, buf);
        MPI_Send(&cost_value, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    MPI_Finalize();

	return statistics;
}

}
} //for closing the namespace