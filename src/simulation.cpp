
#include "simulation.hpp"

/*
 * Function: Manage the master working progress
 */
void master_Code(void) {

    Master master;
    
    clock_t begin_Time = clock();
    
    master.initialise();
    
    
    for (; master.io.iteration < master.io.nsteps; master.io.iteration++) {

        cout << "\nIteration: " << master.io.iteration << endl << endl;
        cout << "Sending data..." << endl;
        cout << ">>> Master sending parameters..." << endl;
        master.send_Parameters();
        cout << ">>> Master sending tetrads..." << endl;
        master.send_Tetrads();
        
        cout << "Calculation..." << endl;
        cout << ">>> Calculating ED & NB forces..." << endl;
        master.force_Calculation();
        
        cout << ">>> Calculating Velocities..." << endl;
        master.cal_Velocities();
        cout << ">>> Calculating Coordinates..." << endl;
        //master.cal_Coordinate();
        cout << ">>> Velocities & Coordinates processing..." << endl;
        master.data_Processing();
        
        cout << "Writing files..." << endl << endl;
        //if (master.io.iteration % master.io.frequency == 0) {
            master.write_Files();
        //}
    }
    
    master.finalise();
    
    clock_t end_Time = clock();
    double time_Usage = double(end_Time - begin_Time) / CLOCKS_PER_SEC;
    cout << "Time usage for simualtion: " << time_Usage << endl << endl;
    
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



