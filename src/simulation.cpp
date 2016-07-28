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
    
    master.initialise();      // Initialises simulaton (reading files, set parameters, etc.)
    
    master.send_Parameters(); // Send parameters to all workers
    
    master.send_Tetrads();    // Send tetrads to all workers
    
    for (int istep = 0, icyc = 0; icyc < 2; icyc++) {//master.io.ncycs; icyc++) {//

        //cout << "\nIteration: " << icyc << endl;
        
        master.generate_Pair_Lists(); // Generate pair lists
        
        for (int i = 0; i < 100; i++) {//master.io.ntsync; i++) {//

            cout << i << endl;
        
            master.cal_Forces();     // Calculate ED/NB forces
            
            master.cal_Velocities(); // Calculate velocities
            
            master.cal_Coordinate(); // Calculate coordinates
    
        }
        
        istep += master.io.ntsync;
        
        master.data_Processing(); // Divide velocities & coordinates by 4
        
        if (istep % master.io.ntwt == 0) { master.write_Info(istep); } // Write energy & trajectory
        
        if (istep % master.io.ntpr == 0) { master.write_Crds(); } // Update crd file
        
    }
    
    master.finalise(); // Finalise simulation, send signal to terminate all workers
    
}



void worker_Code(void) {
    
    Worker worker;
    
    worker.recv_Parameters();   // Receive parameters & initialisation
    
    worker.recv_Tetrads();      // Receive all tetrads
   
    worker.force_Calculation(); // Start ED/NB force calculation
    
}



