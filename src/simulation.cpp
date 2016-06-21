
#include "simulation.hpp"

/*
 * Function: Manage the master working progress
 */
void master_Code(void) {

    Master master;
    
    clock_t begin_Time = clock();
    
    master.initialise();
    
    cout << "\nSending data...\n>>> Master sending parameters..." << endl;
    master.send_Parameters();
    
    cout << ">>> Master sending tetrads..." << endl;
    master.send_Tetrads();
    
    for (int istep = 0, icyc = 0; icyc < master.io.ncycs; icyc++) {

        cout << "\nIteration: " << icyc << endl << endl;
        cout << "Generate pair lists..." << endl;
        master.generate_Pair_Lists();
        
        //for (int i = 0; i < master.io.ntsync; i++) {

            cout << "Calculation...\n>>> Calculating ED & NB forces..." << endl;
            master.cal_Forces();
            
            cout << ">>> Calculating Velocities..." << endl;
            master.cal_Velocities();
            
            cout << ">>> Calculating Coordinates..." << endl;
            master.cal_Coordinate();
            
            cout << ">>> Velocities & Coordinates processing..." << endl;
            master.data_Processing();
    
        //}
        
        istep += master.io.ntsync;
        
        cout << "Writing files..." << endl << endl;
        master.write_Energy(istep);
        master.write_Forces();
        master.write_Trajectories();
        /*
        if (istep % master.io.ntwt == 0) {
            master.write_Energy(istep);
            master.write_Forces();
        }
        if (istep % master.io.ntpr == 0) {
            master.write_Trajectories();
        }*/
        
        master.send_Vels_n_Crds();
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
                
            // Receive the parameters of all tetrads
            } else if (status.MPI_TAG >= TAG_TETRAD && status.MPI_TAG <= TAG_CRDS) {
                worker.recv_Tetrads();
                
            // Receive the velocities and coordinates of tetrads
            } else if (status.MPI_TAG >= TAG_CRDS) {
                worker.recv_Vels_n_Crds();
                
            // Message tag indicates to do the ED force calculation
            } else if (status.MPI_TAG == TAG_ED) {
                worker.ED_Calculation();
                
            // Message tag indicates to do the NB force calculation
            } else if (status.MPI_TAG == TAG_NB) {
                worker.NB_Calculation();
                
            // Indicates to stop work
            } else if (status.MPI_TAG == TAG_DEATH) {
                signal = worker.terminate();
                
            }
        }
    }
    
}



