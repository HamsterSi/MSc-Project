//
//  main.cpp
//  
//
//  Created by Zhuowei Si on 05/04/2016.
//
//

#include <iostream>
#include <string>
#include <vector>
#include "mpi.h"

#include "parameters.hpp"
#include "mpilibrary.hpp"
#include "io.hpp"

using namespace std;

int main(int argc, char *argv[]){
    
    MPI_Comm comm = MPI_COMM_WORLD;
    int i, j, rank, size;
    string prm_File = "./data//GC90c12.prm";
    string crd_File = "./data//GC90_6c.crd";
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    
    if (rank == 0) {
        
        IO io;
        EDMD edmd;
        
        cout << "Master rank initialising MPI program..." << endl;
        cout << Size of MPI processes: " << size << endl;
        
        // Read prm file (Initialise tetrads)
        io.read_Prm(prm_File);
        
        // Read crd file ("true" is for circular DNA and "flase" is for linar DNA)
        io.read_Crd(crd_File, true);

        // Read in initial coordinates
        io.read_Initial_Crds();
        
        cout << "Data read completed." << endl;
        
        // Generate pair lists
        /*for (i = 0; i < io.prm.num_Tetrads; i++) {
            int num_Pairs = 0;
            for (j = 0; j < io.prm.num_Tetrads; j++) {
                
            }
        }*/
        
        /* Statrt MD simulation */
        edmd.parameters_Setting();
        edmd.generate_Stochastic_Term();
        
        clock_t begin_Time = clock();
        
        /* Distrubute ED forces calculation
         * Send tetrads parameters to processes
         * Receive new forces */
        //MPI_Send();
        //MPI_Recv();
        
        /* Distrubut VDW forces calculation
         * Send parameters
         * Receive new forces */
        //MPI_Send();
        //MPI_Recv();
        
        MPI_Barrier(comm);
        
        // Update veolcities & Berendsen temperature control
        edmd.update_Velocities();
        
        // Update coordinates
        edmd.update_Coordinates();
        
        clock_t end_Time = clock();
        double time_Usage = double(end_Time - begin_Time) / CLOCKS_PER_SEC;
        
    } else {
        
        if (rank < ((size-1)/2)) {
            cout << "Worker rank " << rank << "ready for computing ED forces";
            
            MPI_Recv();
            
            // Calculate ED forces
            edmd.calculate_ED_Forces();
            
            MPI_Send();
            
        } else {
            cout << "Worker rank " << rank << "ready for computing VDW forces";
            
            MPI_Recv();
            
            // Calculate vdw forces
            edmd.calculate_VDW_Forces();
            
            MPI_Send();
        }
    }
    
    MPI_Finalize();
    
    return 0;
}






