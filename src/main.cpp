/* 
 * This is the entry of the program. It initialises and finalises 
 * the MPI library and the process pool.
 *
 * The returned status code from fucntion "processPoolInit()" can
 * decide the master process and worker processes, and call
 * according functions such as "master_Code" and "worker_Code" for them.
 *
 * The main funciton also set up the MPI message passing buffer size.
 */

#include <iostream>
#include "mpi.h"

#include "simulation.hpp"

using namespace std;

int main(int argc, char *argv[]){
    
    int rank, size;
    
    // Initialise the MPI environment
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (size < 2) {
        cout << "Requires at least two processes." << endl;
        exit(1);
    }

    // Start processes
    if (rank == 0) master_Code(); // The master initialises the simulation
    if (rank != 0) worker_Code(); // Force calculation workers
    
    // Finalise the MPI environemnt
    MPI_Finalize();
    
    return 0;
}



















