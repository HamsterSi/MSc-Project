/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                  MSc in High Performance Computing, EPCC                     *
 *                       The University of Edinburgh                            *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  main.cpp
 * Brief: The entry of the program. It is responsible for initialising and finalising
 *        the MPI environment. It also divides the MPI processes into the master and
 *        the workers, distribute them different works (by calling different functions)
 */

#include <iostream>
#include "mpi.h"

#include "simulation.hpp"

using namespace std;

int main(int argc, char *argv[]){
    
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // The ED/MD simulation requires at least 2 processes
    if (size < 2) {
        cout << "The program requires at least 2 processes." << endl;
        exit(1);
    }

    // Start simulation. If The MPI rank equals 0 then it is the master process
    // and it calls the function "master_Code()", other MPI processes are workers
    // and they are used to calculate ED & NB forces.
    if (rank == 0) master_Code();
    if (rank != 0) worker_Code();
    
    MPI_Finalize();
    
    return 0;
}



















