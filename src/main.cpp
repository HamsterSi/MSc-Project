/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                 MSc in High Performance Computing, EPCC                      *
 *                      The University of Edinburgh                             *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  main.cpp
 * Brief: The entry of the program. It is responsible for initialising and finalising
 *        the MPI environment, and it also statrts the ED/MD simualtion by calling 
 *        functions for the master process and worker processes.
 */

#include <iostream>
#include "mpi.h"

#include "simulation.hpp"

using namespace std;

int main(int argc, char *argv[]){
    
    int rank, size;
    
    // Initialise the MPI environment
    MPI_Init(&argc, &argv);
    
    // Get MPI rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Requires at least 2 processes
    if (size < 2) {
        cout << "The program requires at least 2 processes." << endl;
        exit(1);
    }

    // Start simulation. If The MPI rank equals 0 then it is the master process
    // and it calls the function "master_Code()", other MPI processes are workers
    // and they are used to calculate ED & NB forces.
    if (rank == 0) master_Code();
    if (rank != 0) worker_Code();
    
    // Finalise the MPI environemnt
    MPI_Finalize();
    
    return 0;
}



















