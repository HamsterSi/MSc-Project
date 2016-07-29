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
 * File:  worker.cpp
 * Brief: Implementation of the Worker class functions
 */

#include "worker.hpp"


Worker::Worker(void) {
    
    comm = MPI_COMM_WORLD;
    
}



Worker::~Worker(void) {
    
    // Deallocate arrays
    for (int i = 0; i < num_Tetrads; i++) {
        tetrad[i].deallocate_Tetrad_Arrays();
    }
    
}



void Worker::recv_Parameters(void) {
    
    double edmd_Para[9];
    MPI_Status status;
    
    // Receive edmd simulation parameters
    MPI_Bcast(edmd_Para, 9, MPI_DOUBLE, 0, comm);
    edmd.initialise(edmd_Para[0], edmd_Para[1], edmd_Para[2], edmd_Para[3],
                    edmd_Para[4], edmd_Para[5], edmd_Para[6], edmd_Para[7]);
    num_Tetrads = (int) edmd_Para[8];
    int * tetrad_Para = new int[2 * num_Tetrads];
    
    // Receive the number of atoms & evecs of all tetrads from master
    MPI_Bcast(tetrad_Para, 2 * num_Tetrads, MPI_INT, 0, comm);
    
    // Initialise tetrad, allocate memory & assignments
    tetrad = new Tetrad[num_Tetrads];
    for (int i = 0; i < num_Tetrads; i++) {
        tetrad[i].num_Atoms = tetrad_Para[2 * i];
        tetrad[i].num_Evecs = tetrad_Para[2*i+1];
        
        tetrad[i].allocate_Tetrad_Arrays();
    }
    
    delete [] tetrad_Para;
    
}



void Worker::recv_Tetrads(void) {
    
    MPI_Datatype MPI_Tetrad;
    
    // Create MPI_Datatype "MPI_Tetrad" for all tetrads
    MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, num_Tetrads, tetrad);
    
    // Receive tetrads from the master
    MPI_Bcast(tetrad, 1, MPI_Tetrad, 0, comm);
    
    // Free the MPI_Datatype "MPI_Tetrad"
    MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);
    
}




void Worker::ED_Calculation(int index[], MPI_Request send_Rqt[], MPI_Request recv_Rqt[]) {
    
    MPI_Status recv_Status;
    
    // Receive coordinates
    MPI_Irecv(tetrad[index[0]].coordinates, 3 * tetrad[index[0]].num_Atoms,
              MPI_DOUBLE, 0, TAG_ED, comm, &recv_Rqt[1]);
    
    MPI_Wait(&recv_Rqt[1], &recv_Status);
    
    // Calculate ED forces & energy
    edmd.calculate_ED_Forces(&tetrad[index[0]]);
    
    // send ED forces, energy back
    MPI_Isend(index, 3, MPI_INT, 0, TAG_INDEX, comm, &send_Rqt[0]);
    
    MPI_Isend(tetrad[index[0]].ED_Forces,   3 * tetrad[index[0]].num_Atoms + 1,
              MPI_DOUBLE, 0, TAG_ED + 1, comm, &send_Rqt[1]);
    MPI_Isend(tetrad[index[0]].coordinates, 3 * tetrad[index[0]].num_Atoms,
              MPI_DOUBLE, 0, TAG_ED + 2, comm, &send_Rqt[2]);
    
}



void Worker::NB_Calculation(int index[], MPI_Request send_Rqt[], MPI_Request recv_Rqt[]) {
    
    MPI_Status recv_Status;
    
    // Receive coordinates
    MPI_Irecv(tetrad[index[0]].coordinates, 3 * tetrad[index[0]].num_Atoms,
              MPI_DOUBLE, 0, TAG_NB + 1, comm, &recv_Rqt[1]);
    MPI_Irecv(tetrad[index[1]].coordinates, 3 * tetrad[index[1]].num_Atoms,
              MPI_DOUBLE, 0, TAG_NB + 2, comm, &recv_Rqt[2]);
    
    MPI_Wait(&recv_Rqt[1], &recv_Status);
    MPI_Wait(&recv_Rqt[2], &recv_Status);
    
    // Calculate NB forces, NB energy & Electrostatic Energy
    edmd.calculate_NB_Forces(&tetrad[index[0]], &tetrad[index[1]]);
    
    // send NB forces, energy back
    MPI_Isend(index, 3, MPI_INT, 0, TAG_INDEX, comm, &send_Rqt[0]);
    
    MPI_Isend(tetrad[index[0]].NB_Forces, 3 * tetrad[index[0]].num_Atoms + 2,
              MPI_DOUBLE, 0, TAG_NB + 1, comm, &send_Rqt[1]);
    MPI_Isend(tetrad[index[1]].NB_Forces, 3 * tetrad[index[1]].num_Atoms + 2,
              MPI_DOUBLE, 0, TAG_NB + 2, comm, &send_Rqt[2]);
    
}



void Worker::force_Calculation(void) {
    
    int signal = 1, index[3];
    
    MPI_Request send_Rqt[3], recv_Rqt[3];
    MPI_Status  send_Status, recv_Status;
    
    // Receive the tetard indexes for ED/NB force calcualtion
    MPI_Irecv(index, 3, MPI_INT, 0, TAG_INDEX, comm, &recv_Rqt[0]);
    MPI_Wait(&recv_Rqt[0], &recv_Status);
    
    switch (index[2]) {
        case TAG_ED: ED_Calculation(index, send_Rqt, recv_Rqt); break;
        case TAG_NB: NB_Calculation(index, send_Rqt, recv_Rqt); break;
        case TAG_SIGNAL: signal = 0; break;
        default: break;
    }
    
    // The loop to take new jobs & send results back to the master
    while (signal == 1) {
        
        MPI_Irecv(index, 3, MPI_INT, 0, TAG_INDEX, comm, &recv_Rqt[0]);
        MPI_Wait(&recv_Rqt[0], &recv_Status);
        
        MPI_Wait(&send_Rqt[0], &send_Status);
        MPI_Wait(&send_Rqt[1], &send_Status);
        MPI_Wait(&send_Rqt[2], &send_Status);
        
        switch (index[2]) {
            case TAG_ED: ED_Calculation(index, send_Rqt, recv_Rqt); break;
            case TAG_NB: NB_Calculation(index, send_Rqt, recv_Rqt); break;
            case TAG_SIGNAL: signal = 0; break;
            default: break;
        }
        
    }
    
}


