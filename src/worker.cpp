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
 * File:  worker.cpp
 * Brief: The implementation of the Worker class functions
 */

#include "worker.hpp"


Worker::Worker(void) {
    
    comm      = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
}



Worker::~Worker(void) {
    
    // Deallocate memory of arrays
    array.deallocate_2D_Double_Array(pair_Lists);
    array.deallocate_2D_Int_Array(ED_Index);
    array.deallocate_2D_Int_Array(NB_Index);
    array.deallocate_2D_Double_Array(NB_Forces);

    // Free the MPI Data type
    for (int i = 0; i < num_Tetrads; i++) {
        tetrad[i].deallocate_Tetrad_Arrays();
        mpi.free_MPI_ED_Forces(&(MPI_ED_Forces[i]));
    }
    mpi.free_MPI_Crds(&MPI_Crds);
    
}



void Worker::recv_Parameters(void) {
    
    double edmd_Para[11];
    
    // Receive edmd simulation parameters
    MPI_Bcast(edmd_Para, 11, MPI_DOUBLE, 0, comm);
    edmd.initialise(edmd_Para[0], edmd_Para[1], edmd_Para[2], edmd_Para[3],
                    edmd_Para[4], edmd_Para[5], edmd_Para[6], edmd_Para[7]);
    num_Tetrads = (int) edmd_Para[8];
    num_Pairs   = (int) edmd_Para[9];
    max_Atoms   = (int) edmd_Para[10];
    int * tetrad_Para = new int[2 * num_Tetrads];
    
    // Receive the tetrad parameters & initialise the tetrad array
    MPI_Bcast(tetrad_Para, 2 * num_Tetrads, MPI_INT, 0, comm);
    tetrad = new Tetrad[num_Tetrads];
    
    for (int i = 0; i < num_Tetrads; i++) {
        tetrad[i].num_Atoms = tetrad_Para[2 * i];
        tetrad[i].num_Evecs = tetrad_Para[2*i+1];
        tetrad[i].allocate_Tetrad_Arrays();
    }
    
    // Create new MPI_Datatype for ED/NB force calcation.
    MPI_ED_Forces = new MPI_Datatype [num_Tetrads]; // For every tetrad
    for (int i = 0; i < num_Tetrads; i++) {
        mpi.create_MPI_ED_Forces(&(MPI_ED_Forces[i]), &(tetrad[i]));
    }
    mpi.create_MPI_Crds(&MPI_Crds, num_Tetrads, tetrad); // For all tetrads
    
    // The parameters for pair lists & workeload
    num_Pairs  = num_Tetrads * (num_Tetrads - 1) / 2;
    pair_Lists = array.allocate_2D_Double_Array(num_Pairs, 2);
    ED_Index   = array.allocate_2D_Int_Array(size - 1, 2);
    NB_Index   = array.allocate_2D_Int_Array(size - 1, 2);
    NB_Forces  = array.allocate_2D_Double_Array(num_Tetrads, 3 * max_Atoms + 2);
    
    delete [] tetrad_Para;
    
}



void Worker::recv_Tetrads(void) {
    
    MPI_Datatype MPI_Tetrad;
    
    mpi.create_MPI_Tetrad(&MPI_Tetrad, num_Tetrads, tetrad);
    
    MPI_Bcast(tetrad, 1, MPI_Tetrad, 0, comm);
    
    mpi.free_MPI_Tetrad(&MPI_Tetrad);
    
}



void Worker::recv_Messages(void) {
    
    int signal = 1, flag;
    MPI_Status recv_Status;
    
    while (signal != TAG_END) {
        
        // Test if there is any message arrived
        MPI_Iprobe(0, MPI_ANY_TAG, comm, &flag, &recv_Status);
        
        if (flag && recv_Status.MPI_TAG >= TAG_PAIRS && recv_Status.MPI_TAG <= TAG_PAIRS + 2) {
            MPI_Recv(&(pair_Lists[0][0]), 2 * num_Pairs, MPI_DOUBLE, 0, TAG_PAIRS, comm, &recv_Status);
            MPI_Recv(&(NB_Index[0][0]), 2 * (size - 1), MPI_INT, 0, TAG_PAIRS + 1, comm, &recv_Status);
            MPI_Recv(&(ED_Index[0][0]), 2 * (size - 1), MPI_INT, 0, TAG_PAIRS + 2, comm, &recv_Status);
        }
        
        else if (flag && recv_Status.MPI_TAG == TAG_FORCE) {
            // Receive the force calculation single
            MPI_Recv(&flag, 1, MPI_INT, 0, TAG_FORCE, comm, &recv_Status);
            
            // Receive the coordinates of all tetrads
            MPI_Bcast(tetrad, 1, MPI_Crds, 0, comm);
            
            // Start the ED/NB force calculation
            force_Calculation();
        }
        
        else if (flag && recv_Status.MPI_TAG == TAG_END){
            MPI_Recv(&signal, 1, MPI_INT, 0, TAG_END, comm, &recv_Status);
        }

    }
    
}




void Worker::force_Calculation() {
    
    int i, j, i1, i2;
    int index = ED_Index[rank - 1][0];
    int workload = ED_Index[rank - 1][1];
    MPI_Request send_Request[workload];
    MPI_Status send_Status[workload];
    
    // Calculate the ED forces and the random terms
    for (i = index; i < index + workload; i++) {
        
        edmd.calculate_ED_Forces(&(tetrad[i]));
        edmd.calculate_Random_Terms(&(tetrad[i]), rank);
        
        MPI_Isend(&(tetrad[i]), 1, MPI_ED_Forces[i], 0, TAG_ED + i, comm, &(send_Request[i - index]));
        
    }
    
    // Calculate the NB forces
    empty_NB_Forces();
    for (i = NB_Index[rank - 1][0]; i < NB_Index[rank - 1][0] + NB_Index[rank - 1][1]; i++) {
        
        i1 = pair_Lists[i][0];
        i2 = pair_Lists[i][1];
        edmd.calculate_NB_Forces(&tetrad[i1], &tetrad[i2]);
        
        // Sum up the NB forces of the specific tetrads
        for (j = 0; j < 3 * tetrad[i1].num_Atoms + 2; j++) {
            NB_Forces[i1][j] += tetrad[i1].NB_Forces[j];
        }
        for (j = 0; j < 3 * tetrad[i2].num_Atoms + 2; j++) {
            NB_Forces[i2][j] += tetrad[i2].NB_Forces[j];
        }
    }
    
    // Reduce & sum up the NB forces to the master
    MPI_Reduce(&(NB_Forces[0][0]), &(NB_Forces[0][0]), num_Tetrads * (3 * max_Atoms + 2), MPI_DOUBLE, MPI_SUM, 0, comm);
    
    // Wait all ED forces to be received
    MPI_Waitall(workload, send_Request, send_Status);
    
}



void Worker::empty_NB_Forces(void) {
    
    for (int i = 0; i < num_Tetrads; i++) {
        for (int j = 0; j < 3 * tetrad[i].num_Atoms + 2; j++) {
            NB_Forces[i][j] = 0.0;
        }
    }
    
}




