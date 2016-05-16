
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
    
    master.data_Sending();
    
    master.tetrad_Sending();
    
    master.force_Passing();
    
    //master.LV_Forces();
    
    master.total_Forces();
    
    master.velocities();
    
    master.coordinates();
    
    clock_t end_Time = clock();
    double time_Usage = double(end_Time - begin_Time) / CLOCKS_PER_SEC;
    cout << "Time usage for the simualtion: " << time_Usage << endl;
    
}

/*
 * Function: Manage the worker working progress
 */
void worker_Code()
{
    int rank, flag, signal = 1;
    Worker_Management worker;
    
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
    MPI_Comm_rank(comm, &rank);
    
    while (signal)
    {
        MPI_Iprobe(0, MPI_ANY_TAG, comm, &flag, &status);
        if (flag)
        {
            switch (status.MPI_TAG) {
                case TAG_DATA:
                    cout << "Rank " << rank << "receives data" << endl;
                    worker.data_Receiving();
                    break;
                    
                case TAG_TETRAD:
                    cout << "Rank " << rank << "receives tetrads" << endl;
                    worker.tetrad_Receiving();
                    break;
                    
                case TAG_ED:
                    cout << "Rank " << rank << "computes ED forces" << endl;
                    worker.ED_Calculation();
                    break;
                    
                case TAG_NB:
                    cout << "Rank " << rank << "computes NB forces" << endl;
                    worker.NB_Calculation();
                    break;
                    
                case TAG_DEATH:
                    signal = 0;
                    break;
            }
        }
    }
}



