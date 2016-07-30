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
    
    // Initialises simulaton (reading files, set parameters, etc.)
    master.initialise();
    
    // Master sends simulation parameters & tetrads
    master.send_Parameters();
    master.send_Tetrads();
    
    for (int istep = 0; istep < master.io.nsteps; istep += master.io.ntsync) {//1000; istep += master.io.ntsync) {//

        // Generate pair lists for NB force calculation
        master.generate_Pair_Lists();
        
        // Loop to calculate forces, then updating the velocities & coordinates
        for (int i = 0; i < master.io.ntsync; i++) {//1; i++) {//
            cout << i << endl;
            
            master.calculate_Forces();
            master.update_Velocities();
            master.update_Coordinates();
    
        }
        
        // Divide velocities & coordinates by 4
        master.data_Processing();
        
        // Write out the energy & temperature, trajectory, update crd file
        if (istep % master.io.ntwt == 0) { master.write_Info(istep); }
        if (istep % master.io.ntpr == 0) { master.write_Crds(); }
        
    }
    
    // Finalise simulation, terminate all workers
    master.finalise();
    
    double time_Usage = double (clock() - start_Time) / CLOCKS_PER_SEC;
    cout << "Time usage of simualtion: " << time_Usage << endl << endl;
    
}



void worker_Code(void) {
    
    Worker worker;
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Receive parameters & initialisation
    worker.recv_Parameters();
    
    // Receive all tetrads
    worker.recv_Tetrads();
   
    // Start ED/NB force calculation
    if (rank == 1) { worker.NB_Force_Processing(); }
    else { worker.force_Calculation(); }
    
}



