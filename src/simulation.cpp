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
    
    for (int istep = 0, icyc = 0; icyc < 1; icyc++) {//master.io.ncycs; icyc++) {//

        //cout << "\nIteration: " << icyc << endl;
        
        master.generate_Pair_Lists(); // Generate pair lists
        
        for (int i = 0; i < master.io.ntsync; i++) {//1; i++) {//

            cout << i << endl;
        
            master.cal_Forces(); // Master sends tetrad indexes & coordinates to workers, worker send ED/NB forces back
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
    cout << "Simulation ended." << endl;
    
    double time_Usage = double (clock() - start_Time) / CLOCKS_PER_SEC;
    cout << "Time usage of simualtion: " << time_Usage << endl << endl;
    
}



void worker_Code(void) {
    
    int flag, signal = 1;
    Worker worker;
    MPI_Status status;
    MPI_Request send_Request, recv_Request;
    
    // Receive parameters & initialisation
    worker.recv_Parameters();
    
    // Receive all tetrads
    worker.recv_Tetrads();
    
    // Receive tetrad indexes & coordinates for Ed/NB calculation
    MPI_Irecv(&(worker.recv_Buf[0][0]), 2 * (3 * worker.max_Atoms + 2), MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_Request);
    MPI_Wait(&recv_Request, &status);
    
    if (status.MPI_TAG == TAG_ED) { worker.ED_Calculation(&send_Request); }
    if (status.MPI_TAG == TAG_NB) { worker.NB_Calculation(&send_Request); }
    
    while (signal == 1) {
        
        MPI_Irecv(&(worker.recv_Buf[0][0]), 2 * (3 * worker.max_Atoms + 2), MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_Request);
        
        MPI_Wait(&send_Request, &status);
        MPI_Wait(&recv_Request, &status);
        
        if (status.MPI_TAG == TAG_ED) { worker.ED_Calculation(&send_Request); }
        if (status.MPI_TAG == TAG_NB) { worker.NB_Calculation(&send_Request); }
        if (status.MPI_TAG == TAG_SIGNAL) { signal = 0; }
        
    }
    
}



