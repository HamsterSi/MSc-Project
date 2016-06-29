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
    
    clock_t start_Time = clock();
    
    // Master initialises the simulaton, including reading data from input files,
    // allocate memory for arrays, set the frequencies of iterations, etc.
    master.initialise();
    
    cout << "\nMaster sending data to Workers." << endl;
    cout << ">>> Sending parameters..." << endl;
    master.send_Parameters();
    
    cout << ">>> Sending tetrads..." << endl;
    master.send_Tetrads();
    cout << "Data sending completed." << endl;
    cout << "\nStart simulation...\n" << endl;
    
    for (int istep = 0, icyc = 0; icyc < 2; icyc++) {//master.io.ncycs; icyc++) {//

        cout << "\nIteration: " << icyc << endl;
        
        // Generate pair lists
        master.generate_Pair_Lists();
        
        for (int i = 0; i < master.io.ntsync; i++) {//1; i++) {//

            cout << i << endl;
        
            // Master sends coordinates and tetrad indexes to workers and worker
            // calculate ED/NB forces for tetrads
            master.cal_Forces();
            
            // Master calculates velocities for all tetrads
            master.cal_Velocities();
            
            // Master calculates coordinates for all tetrads
            master.cal_Coordinate();
    
        }
        
        istep += master.io.ntsync;
        
        // The velocities & coordinates of tetrads are gathered and processed
        master.data_Processing();
        
        // Write out file with certain frequencies
        if (istep % master.io.ntwt == 0) {
            
            // Write out the energies and tmeperature
            master.write_Energy(istep);
            
            //master.write_Forces();
            
            //master.write_Trajectories(istep-master.io.ntsync);
        }
        if (istep % master.io.ntpr == 0) {
            
            // Write out a new crd file
            master.write_Crds();
        }
        
    }
    
    // Finalise the simulation, send signal to workers to stop their work
    master.finalise();
    cout << "Simulation ended." << endl;
    
    double time_Usage = double (clock() - start_Time) / CLOCKS_PER_SEC;
    cout << "Time usage of simualtion: " << time_Usage << endl << endl;
    
}



void worker_Code(void) {
    
    int flag, signal = 1;
    Worker worker;
    MPI_Status status;
    
    worker.recv_Parameters();
    
    worker.recv_Tetrads();
    
    while (signal == 1) {
        
        MPI_Recv(&(worker.buffer[0][0]), 2 * (3 * worker.max_Atoms + 2), MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        if (status.MPI_TAG == TAG_ED) {
            worker.ED_Calculation();
            
        } else if (status.MPI_TAG == TAG_NB) {
            worker.NB_Calculation();
            
        } else if (status.MPI_TAG == TAG_SIGNAL) {
            signal = 0;
            
        }
    }

    /*
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
            } else if (status.MPI_TAG == TAG_SIGNAL) {
                signal = worker.terminate();
                
            }
        }
    }*/
    
}



