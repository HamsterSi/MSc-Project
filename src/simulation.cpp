
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
    
    master.tetrad_Received_Signal();
    
    master.force_Passing();
    
    master.calculate_Total_Forces();
    
    master.velocities();
    
    master.coordinates();
    
    master.finalise();
    
    clock_t end_Time = clock();
    double time_Usage = double(end_Time - begin_Time) / CLOCKS_PER_SEC;
    cout << "Time usage for the simualtion: " << time_Usage << "\n" << endl;
    
}

/*
 * Function: Manage the worker working progress
 */
void worker_Code() {
    int flag, signal = 1;
    
    MPI_Status status;
    Worker_Management worker;
    
    while (signal)
    {
        MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        if (flag)
        {
            if (status.MPI_TAG == TAG_DATA) {
                worker.parameters_Receiving();
                
            } else if (status.MPI_TAG >= TAG_TETRAD && status.MPI_TAG < TAG_ED) {
                worker.tetrads_Receiving();
                
            } else if (status.MPI_TAG == TAG_ED) {
                worker.ED_Calculation();
                
            } else if (status.MPI_TAG == TAG_NB) {
                worker.NB_Calculation();
                
            } else if (status.MPI_TAG == TAG_DEATH) {
                signal = 0;
                
            }
        }
    }
    
}



