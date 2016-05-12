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

#include "pool/pool.h"
#include "simulation.hpp"
#include "mpilibrary.hpp"

using namespace std;

int main(int argc, char *argv[]){
    
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size, status_Code, buffer_Size = 40960;
    
    // Initialise the parallel environment
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &size);
    
    MPI_Buffer_attach(malloc(buffer_Size), buffer_Size);
    //status_Code = processPoolInit();
    
    if (size < 2) {
        cout << "Requires at least two processes." << endl;
        exit(1);
    }

    // Start the processes
    if (status_Code == 2) master_Code(); // The master initialises the simulation
    if (status_Code == 1) worker_Code(); // Force calculation workers
    
    // Finalise the environemnt
    //processPoolFinalise();
    MPI_Finalize();
    
    return 0;
}



















