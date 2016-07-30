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

    // Deallocate tetrad memeory
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




void Worker::ED_Calculation(int index[], MPI_Request send_Request[]) {
    
    // Calculate the ED forces & ED energy
    edmd.calculate_ED_Forces(&tetrad[index[0]]);
    
    // Send tetrad index back to master to tell it the calculation's finished
    MPI_Isend(index, 3, MPI_INT, 0, TAG_INDEX, comm, &send_Request[0]);
    
    // Send tetrad index, ED forces & coordinates to rank 1
    MPI_Isend(index, 3, MPI_INT, 1, TAG_INDEX, comm, &send_Request[1]);
    MPI_Isend(tetrad[index[0]].ED_Forces, 3 * tetrad[index[0]].num_Atoms + 1,
              MPI_DOUBLE, 1, TAG_ED + 1, comm, &send_Request[2]);
    MPI_Isend(tetrad[index[0]].coordinates,   3 * tetrad[index[0]].num_Atoms,
              MPI_DOUBLE, 1, TAG_ED + 2, comm, &send_Request[3]);
    
}



void Worker::NB_Calculation(int index[], MPI_Request send_Request[]) {
    
    // Calculate the NB forces, NB energy & Electrostatic energy
    edmd.calculate_NB_Forces(&tetrad[index[0]], &tetrad[index[1]]);
    
    // Send tetrad indexes back to master to tell it the calculation's finished
    MPI_Isend(index, 3, MPI_INT, 0, TAG_INDEX, comm, &send_Request[0]);
    
    // Send tetrad indexes, two NB forces to rank 1 to sum up
    MPI_Isend(index, 3, MPI_INT, 1, TAG_INDEX, comm, &send_Request[1]);
    MPI_Isend(tetrad[index[0]].NB_Forces, 3 * tetrad[index[0]].num_Atoms + 2,
              MPI_DOUBLE, 1, TAG_NB + 1, comm, &send_Request[2]);
    MPI_Isend(tetrad[index[1]].NB_Forces, 3 * tetrad[index[1]].num_Atoms + 2,
              MPI_DOUBLE, 1, TAG_NB + 2, comm, &send_Request[3]);
    
}



void Worker::force_Calculation(void) {
    
    int signal = 1, index[3];
    MPI_Request send_Request[4];
    MPI_Status send_Status, recv_Status;
    
    // Receive the tetrad index(es) & force tag from master
    MPI_Recv(index, 3, MPI_INT, 0, TAG_INDEX, comm, &recv_Status);
    
    // Calculation performed according to the force tag
    switch (index[2]) {
        case TAG_ED: ED_Calculation(index, send_Request); break;
        case TAG_NB: NB_Calculation(index, send_Request); break;
        case TAG_SIGNAL: signal = 0; break;
        default: break;
    }
    
    // The loop to receive new ED/NB force calcualtion task & calculation
    while (signal == 1) {
        
        // Receive new tetrad index(es) & force tag from master
        MPI_Recv(index, 3, MPI_INT, 0, TAG_INDEX, comm, &recv_Status);
        
        // Wait the previous send request to be finished.
        MPI_Wait(&(send_Request[0]), &send_Status);
        MPI_Wait(&(send_Request[1]), &send_Status);
        MPI_Wait(&(send_Request[2]), &send_Status);
        MPI_Wait(&(send_Request[3]), &send_Status);
        
        // ED/NB force calculation according to the the force tag
        switch (index[2]) {
            case TAG_ED: ED_Calculation(index, send_Request); break;
            case TAG_NB: NB_Calculation(index, send_Request); break;
            case TAG_SIGNAL: signal = 0; break;
            default: break;
        }
        
    }
    
}



void Worker::recv_ED_Forces(int index, int source) {
    
    MPI_Status recv_Status;
    
    // Receive the ED forces & coordinates from other workers.
    MPI_Recv(tetrad[index].ED_Forces, 3 * tetrad[index].num_Atoms + 1,
             MPI_DOUBLE, source, TAG_ED + 1, comm, &recv_Status);
    MPI_Recv(tetrad[index].coordinates,   3 * tetrad[index].num_Atoms,
             MPI_DOUBLE, source, TAG_ED + 2, comm, &recv_Status);
    
}



void Worker::empty_NB_Forces(void) {
    
    for (int i = 0; i < num_Tetrads; i++) {
        for (int j = 0; j < 3 * tetrad[i].num_Atoms + 2; j++) {
            tetrad[0].NB_Forces[j] = 0.0;
        }
    }
    
}



void Worker::sum_NB_Forces(int index[], int source, double* NB_Forces1, double* NB_Forces2) {

    MPI_Request recv_Request[2];
    MPI_Status  recv_Status;
    
    // Receive the two arrays of NB forces form other workers
    MPI_Irecv(NB_Forces1, 3 * tetrad[index[0]].num_Atoms + 2, MPI_DOUBLE,
              source, TAG_NB + 1, comm, &recv_Request[0]);
    MPI_Irecv(NB_Forces2, 3 * tetrad[index[1]].num_Atoms + 2, MPI_DOUBLE,
              source, TAG_NB + 2, comm, &recv_Request[1]);
    
    // Wait one NB forces array to arrive & sum up
    MPI_Wait(&recv_Request[0], &recv_Status);
    for (int i = 0; i < 3 * tetrad[index[0]].num_Atoms + 2; i++) {
        tetrad[index[0]].NB_Forces[i] += NB_Forces1[i];
    }
    
    // Wait one NB forces array to arrive & sum up
    MPI_Wait(&recv_Request[1], &recv_Status);
    for (int i = 0; i < 3 * tetrad[index[1]].num_Atoms + 2; i++) {
        tetrad[index[1]].NB_Forces[i] += NB_Forces2[i];
    }
    
}



void Worker::clip_NB_Forces(void) {
    
    double max_Forces = 1.0;
    
    // Clip all NB forces between -1.0 and 1.0
    for (int i  = 0; i < num_Tetrads; i++) {
        
        for (int j = 0; j < 3 * tetrad[i].num_Atoms; j++) {
            tetrad[i].NB_Forces[j] = min( max_Forces, tetrad[i].NB_Forces[j]);
            tetrad[i].NB_Forces[j] = max(-max_Forces, tetrad[i].NB_Forces[j]);
        }
        
    }
    
}



void Worker::force_Processing(void) {
    
    int index[3], signal = 1, max_Atoms = 0;
    
    // Get the maximum number of atoms in tetrads to allocate MPI buffer
    for (int i = 0; i < num_Tetrads; i++) {
        max_Atoms = max(max_Atoms, tetrad[i].num_Atoms);
    }
    
    // Allocate the MPi buffer for receivin Nb forces
    double * NB_Forces1 = new double [3 * max_Atoms + 2];
    double * NB_Forces2 = new double [3 * max_Atoms + 2];
    
    MPI_Status recv_Status;
    MPI_Datatype MPI_Force;
    MPI_Library::create_MPI_Force(&MPI_Force, num_Tetrads, tetrad);
    
    double temp[3] = {0.0, 0.0, 0.0};
    
    while (signal == 1) {
        
        MPI_Recv(index, 3, MPI_INT, MPI_ANY_SOURCE, TAG_INDEX, comm, &recv_Status);
        
        switch (index[2]) {
            case TAG_ED    : recv_ED_Forces(index[0], recv_Status.MPI_SOURCE);
                             edmd.calculate_Random_Forces(&tetrad[index[0]]); break;
            case TAG_CLEAN : empty_NB_Forces(); break;
            case TAG_NB    : sum_NB_Forces(index, recv_Status.MPI_SOURCE, NB_Forces1, NB_Forces2); break;
            case TAG_FORCE : clip_NB_Forces(); MPI_Send(tetrad, 1, MPI_Force, 0, TAG_FORCE, comm); break;
            case TAG_SIGNAL: signal = 0; break;
            default: break;
        }
        
    }
    
    delete [] NB_Forces1;
    delete [] NB_Forces2;
    MPI_Library::free_MPI_Force(&MPI_Force);
    
}


