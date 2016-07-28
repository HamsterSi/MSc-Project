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
    
    cout << "\nMaster sending data to workers.\n>>> Sending parameters..." << endl;
    master.send_Parameters();
    
    cout << ">>> Sending tetrads..." << endl;
    master.send_Tetrads();
    cout << "Data sending completed.\n\nStart simulation...\n" << endl;
    
    for (int istep = 0, icyc = 0; icyc < 1; icyc++) {//master.io.ncycs; icyc++) {//

        //cout << "\nIteration: " << icyc << endl;
        
        master.generate_Pair_Lists(); // Generate pair lists
        
        for (int i = 0; i < 100; i++) {//master.io.ntsync; i++) {//

            cout << i << endl;
        
            master.cal_Forces();     // Calculate ED/NB forces
            master.cal_Velocities(); // Calculate velocities
            master.cal_Coordinate(); // Calculate coordinates
    
        }
        
        istep += master.io.ntsync;
        
        // Process the velocities & coordinates of tetrads together
        master.data_Processing();
        
        // Write out energies & tmeperature, and the trajectories of DNA
        if (istep % master.io.ntwt == 0) {
            master.write_Energy(istep);
            master.write_Trajectories(istep-master.io.ntsync);
        }
        
        // Write out a new crd file
        if (istep % master.io.ntpr == 0) {
            master.write_Crds();
        }
        
    }
    
    // Finalise the simulation, send signal to workers to stop their work
    master.finalise();

    double time_Usage = double (clock() - start_Time) / CLOCKS_PER_SEC;
    cout << "Simulation ended.\nTime usage of simualtion: " << time_Usage << endl << endl;
    
}



void worker_Code(void) {
    
    Worker worker;
    
    // Receive parameters & initialisation
    worker.recv_Parameters();
    
    // Receive all tetrads
    worker.recv_Tetrads();
    
    // Start ED/NB force calculation
    worker.force_Calculation();
    
}



