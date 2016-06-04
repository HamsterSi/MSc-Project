
#include "simulation.hpp"

/*
 * Function: Manage the master working progress
 */
void master_Code(void) {
    
    Master master;
    
    clock_t begin_Time = clock();
    
    master.initialise();
    
    master.send_Parameters();
    master.send_Tetrads();
    
    master.force_Passing();
    
    master.cal_Velocities();
    master.cal_Coordinate();
    master.data_Processing();
    
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
    Worker worker;
    MPI_Status status;
    
    while (signal == 1)
    {
        // Test if there is any message arrived
        MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        if (flag)
        {
            // Arrived message tag indicates to receive parameters
            if (status.MPI_TAG == TAG_DATA) {
                worker.recv_Parameters();
                
            // Arrived message tag indicates to receive tetrads
            } else if (status.MPI_TAG >= TAG_TETRAD) {
                worker.recv_Tetrads();
                
            // Arrived message tag indicates to do the ED force calculation
            } else if (status.MPI_TAG == TAG_ED) {
                worker.ED_Calculation();
                
            // Arrived message tag indicates to do the NB force calculation
            } else if (status.MPI_TAG == TAG_NB) {
                worker.NB_Calculation();
                
            // Arrived message tag indicates to stop work
            } else if (status.MPI_TAG == TAG_DEATH) {
                signal = worker.terminate();
                
            }
        }
    }
    
}



