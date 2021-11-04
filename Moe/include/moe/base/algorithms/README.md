For Greg:

RUNNING THE CODE:

    Everything to be done in dphpc/project_codes/Moe
    cmake . then make

    mpirun paps_functions with #iterations as input
    eg: mpirun -np 2 paps_functions 1000


ABOUT ParallelAsynchronousPSO.hpp:

    We implemented a paralell asynchronous PSO following the algorithm from https://rcnl.rice.edu/PDFs/ijnme2006.pdf
    
    Rank == 0 is the master process. It updates positions and velocities, checks bounds for new positions and checks wether new value is a new best value
    Rank != 0 are slave processes. They do the function evaluation.

The program is currently a work in progress

    First we scatter a particle to each slave process.
    Then in the master we start a while loop (that makes shure we roughly evaluate each particle iteration times)
        we recieve a particle with a new fitness and based on that update the position and velocity
        we check bounds
        we send a new particle to a slave
        we check if new position is the new local & global best
    Then we have a termination for loop
        usually the Tag in the MPI_SEND is the particle index + 1. In this loop we send tag 0 and that kills slave processes

    
    In Slave processes we do evaluation on the scattered particle
    Then enter an loop
        Revieve a particle
        Check if tag == 0: if yes exit loop
        evaluate particle and send
