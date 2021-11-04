
#ifndef DO_PARALELL_PSO_H_
#define DO_PARALELL_PSO_H_

#include <mpi.h>

#include "general/general_objects.h"
#include "statistics/statistics.h"

namespace highprecisionpso {
/**
* @brief This method proceeds the activities of the PSO algorithm.
*
* Here the PSO iterations are controlled.  At each iteration the update of each
* position of the particles is started.  The actual update of the global
* attractor is initialized here.  Statistical evaluations, the swarm debugging,
* backups and the check whether the program is allowed to run are started here.
*
* @param statistics The statistics where data about the current PSO execution
* are stored.
*
* @return The final statistic.
*/
Statistics* DoParalellPso(Statistics* statistics, time_t& LAST_BACKUP, time_t& LAST_RUN_CHECK, int& argc, char * argv[]);

}

#endif
