
#include "simulation.hpp"

/*
 * Function: Manage the master working progress
 */
void master_Code(void) {
    
    Master_Management master;
    
    clock_t begin_Time = clock();
    
    master.initialise();
    
    master.parameters_Sending();
    
    master.tetrads_Sending();
    
    master.force_Calculation();
    
    master.velocity_Calculation();
    
    master.coordinate_Calculation();
    
    master.finalise();
    
    clock_t end_Time = clock();
    double time_Usage = double(end_Time - begin_Time) / CLOCKS_PER_SEC;
    cout << "Time usage for the simualtion: " << time_Usage << "\n" << endl;
    
}

/*
 * Function: Manage the worker working progress
 */
void worker_Code() {
    
    int flag, signal;
    Worker_Management worker;
    MPI_Status status;
    
    while (1)
    {
        MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        if (flag)
        {
            if (status.MPI_TAG == TAG_DATA) {
                worker.parameters_Receiving();
                
            } else if (status.MPI_TAG >= TAG_TETRAD) {
                worker.tetrads_Receiving();
                
            } else if (status.MPI_TAG == TAG_ED) {
                worker.ED_Calculation();
                
            } else if (status.MPI_TAG == TAG_NB) {
                worker.NB_Calculation();
                
            } else if (status.MPI_TAG == TAG_DEATH) {
                MPI_Recv(&signal, 1, MPI_INT, 0, TAG_DEATH, MPI_COMM_WORLD, &status);
                break;
                
            }
        }
    }
    
}



