
#include "simulation.hpp"

/*
 * Function: Manage the master working progress
 */
void master_Code(void) {
    
    MPI_Comm comm = MPI_COMM_WORLD;
    int i, j, rank, size;
    
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    cout << "Master rank initialising MPI program..." << endl;
    cout << "Size of MPI processes: " << size << endl;
    
    Master_Management master;
    
    clock_t begin_Time = clock();
    
    master.master_ED_Forces();
    
    master.master_NB_Forces();
    
    master.master_LV_Forces();
    
    master.master_Total_Forces();
    
    master.master_Velocities();
    
    master.master_Coordinates();
    
    clock_t end_Time = clock();
    double time_Usage = double(end_Time - begin_Time) / CLOCKS_PER_SEC;
    cout << "Time usage for the simualtion: " << time_Usage << endl;
}

/*
 * Function: Manage the worker working progress
 */
void worker_Code()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size, signal, worker_Status = 1;
    
    Worker_Management worker;
    
    while (worker_Status)
    {
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
        
        if (rank <= ((size+1)/2)) {
            cout << "Rank " << rank << "simulates ED forces" << endl;
            signal = worker.worker_ED_Forces();
            if (signal == 1) worker_Status = 1;//workerSleep();
            
        } else {
            cout << "Rank " << rank << "simulates NB forces" << endl;
            signal = worker.worker_NB_Forces();
            if (signal == 1) worker_Status = 1;//workerSleep();
            
        }
    }
}

