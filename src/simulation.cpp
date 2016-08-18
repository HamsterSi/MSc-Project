/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                  MSc in High Performance Computing, EPCC                     *
 *                       The University of Edinburgh                            *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  simulation.cpp
 * Brief: The implementation of the two functions declared
 */

#include "simulation.hpp"
 
 
void master_Code(void) {

    Master master;
    
    clock_t start_Time = clock();
    
    master.initialise();
    
    master.send_Parameters();
    master.send_Tetrads();
    
    for (int istep = 0; istep < master.io.nsteps; istep += master.io.ntsync) {

        cout << "istep: " << istep << endl;
        
        master.generate_Pair_Lists(); 
        master.generate_Indexes();
        master.send_Workload_Indexes();
        
        for (int i = 0; i < master.io.ntsync; i++) {
            
            cout << i << endl;
            
            master.calculate_Forces();
            master.update_Velocity();
            master.update_Coordinate();
    
        }
    
        master.merge_Vels_n_Crds();
        
        if (istep % master.io.ntwt == 0) { master.write_Info(istep); }
        if (istep % master.io.ntpr == 0) { master.write_Crds(); }
        
    }
    
    master.finalise();
    
    double time_Usage = double (clock() - start_Time) / CLOCKS_PER_SEC;
    cout << "Simulation ended.\nTime usage of simualtion: " << time_Usage << endl << endl;
	
}



void worker_Code(void) {
    
    Worker worker;

    worker.recv_Parameters();
    worker.recv_Tetrads();
    
    worker.recv_Messages();
    
}




