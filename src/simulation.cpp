/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                 MSc in High Performance Computing, EPCC                      *
 *                      The University of Edinburgh                             *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  simulation.cpp
 * Brief: Implementation of two functions for the master process and worker processes
 */

#include "simulation.hpp"


void master_Code(void) {

    Master master;
    
    clock_t begin_Time = clock();
    
    master.initialise();
    
    //cout << "\nSending data...\n>>> Master sending parameters..." << endl;
    master.send_Parameters();
    
    //cout << ">>> Master sending tetrads..." << endl;
    master.send_Tetrads();
    
    for (int istep = 0, icyc = 0; icyc < master.io.ncycs; icyc++) {

        cout << "\nIteration: " << icyc << endl;
        //cout << "\n\nGenerate pair lists..." << endl;
        
        // Generate pair lists
        master.generate_Pair_Lists();
        
        for (int i = 0; i < master.io.ntsync; i++) {

            cout << i << endl;
            //cout << "Calculation...\n>>> Calculating ED & NB forces..." << endl;
            master.cal_Forces();
            
            //cout << ">>> Calculating Velocities..." << endl;
            master.cal_Velocities();
            
            //cout << ">>> Calculating Coordinates..." << endl;
            master.cal_Coordinate();
    
        }
        
        istep += master.io.ntsync;
        
        //cout << ">>> Processing Velocities & Coordinates..." << endl;
        master.data_Processing();
        
        //cout << "Writing files..." << endl << endl;
        if (istep % master.io.ntwt == 0) {
            master.write_Energy(istep);
            //master.write_Forces();
            //master.write_Trajectories(istep-master.io.ntsync);
        }
        if (istep % master.io.ntpr == 0) {
            master.write_Crds();
        }
        
    }
    
    master.finalise();
    
    clock_t end_Time = clock();
    double time_Usage = double(end_Time - begin_Time) / CLOCKS_PER_SEC;
    cout << "Time usage for simualtion: " << time_Usage << endl << endl;
    
}



void worker_Code(void) {
    
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
            } else if (status.MPI_TAG >= TAG_TETRAD) {
                worker.recv_Tetrads();

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



