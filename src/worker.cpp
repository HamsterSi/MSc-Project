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
    array.deallocate_2D_Array(send_Buf);
    array.deallocate_2D_Array(recv_Buf);
    
    for (int i = 0; i < num_Tetrads; i++) {
        array.deallocate_Tetrad_Arrays(&tetrad[i]);
    }
    
}



void Worker::recv_Parameters(void) {
    
    int i, parameters[2], signal = 1;
    double edmd_Para[8];
    MPI_Status status;
    
    // Receive parameters from the master process
    MPI_Recv(parameters, 2, MPI_INT, 0, TAG_DATA, comm, &status);
    num_Tetrads = parameters[0]; // The number of tetrads
    max_Atoms   = parameters[1]; // The maximum number of atoms in tetrads
    
    send_Buf = array.allocate_2D_Array(2, 3 * max_Atoms + 2);
    recv_Buf = array.allocate_2D_Array(2, 3 * max_Atoms + 2);
    int * tetrad_Para = new int[2 * num_Tetrads];
    
    // Receive edmd parameters & assignment
    MPI_Recv(edmd_Para, 8, MPI_DOUBLE, 0, TAG_DATA, comm, &status);
    edmd.initialise(edmd_Para[0], edmd_Para[1], edmd_Para[2], edmd_Para[3],
                    edmd_Para[4], edmd_Para[5], edmd_Para[6], edmd_Para[7]);
    
    // Receive the number of atoms & evecs in tetrads
    MPI_Recv(tetrad_Para, 2 * num_Tetrads, MPI_INT, 0, TAG_DATA, comm, &status);
    
    // Assign values & allocate memory for all tetrads
    tetrad = new Tetrad[num_Tetrads];
    for (i = 0; i < num_Tetrads; i++) {
        tetrad[i].num_Atoms = tetrad_Para[2 * i];
        tetrad[i].num_Evecs = tetrad_Para[2*i+1];
        
        array.allocate_Tetrad_Arrays(&tetrad[i]);
    }
    
    delete [] tetrad_Para;
    
    // Send feedback to master that has received all parameters.
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_DATA, comm);
    
}



void Worker::recv_Tetrads(void) {
    
    int i, signal = 1;
    MPI_Datatype MPI_Tetrad;
    MPI_Status status;
    
    // Receive all tetrads parameters from the master process
    for (i = 0; i < num_Tetrads; i++) {
        
        MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &tetrad[i]);
        MPI_Recv(&tetrad[i], 1, MPI_Tetrad, 0, TAG_TETRAD+i, comm, &status);
        MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);
        
    }
    
    // Send feedback to master that has received all tetrads
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_TETRAD, comm);
    
}



void Worker::ED_Calculation(void) {
    
    // Assign values for tetrad indexes and coordinates
    int index = send_Buf[0][3 * max_Atoms + 1] = (int) recv_Buf[0][3 * max_Atoms + 1];
    
    // Calculate ED forces (ED energy), Parameters:
    // PARAMETERS: Tetrad* tetrad, double* old_Crds, double* ED_Forces, double* new_Crds, int atoms
    edmd.calculate_ED_Forces(&tetrad[index], recv_Buf[0], send_Buf[0], send_Buf[1], 3 * max_Atoms);
    
}



void Worker::NB_Calculation(void) {
    
    // Assign tetrad indexes from recv buffer to send buffer
    int i = send_Buf[0][3 * max_Atoms + 1] = (int) recv_Buf[0][3 * max_Atoms + 1];
    int j = send_Buf[1][3 * max_Atoms + 1] = (int) recv_Buf[1][3 * max_Atoms + 1];
    
    // Calculate NB forces, NB energy & Electrostatic Energy
    edmd.calculate_NB_Forces(&tetrad[i], &tetrad[j], recv_Buf, send_Buf, 3 * max_Atoms);
    
}



void Worker::force_Calculation(void) {
    
    int signal = 1, num_Buf = 2 * (3 * max_Atoms + 2);
    MPI_Status send_Status, recv_Status;
    MPI_Request send_Request, recv_Request;
    
    // Receive the first tetrad index(es) & coordinates from master for ED/NB calculation
    MPI_Recv(&(recv_Buf[0][0]), num_Buf, MPI_DOUBLE, 0, MPI_ANY_TAG, comm, &status);
    
    switch (recv_Status.MPI_TAG) {
        case TAG_ED: ED_Calculation(); // ED force calculation & send ED forces, energy back
            MPI_Isend(&(send_Buf[0][0]), num_Buf, MPI_DOUBLE, 0, TAG_ED, comm, &send_Request);
            break;
        case TAG_NB: NB_Calculation(); // NB force calculation & send NB forces, energy back
            MPI_Isend(&(send_Buf[0][0]), num_Buf, MPI_DOUBLE, 0, TAG_NB, comm, &send_Request);
            break;
        default: break;
    }
    
    // The loop to take new jobs & send results back to the master
    while (signal == 1) {
        
        MPI_Irecv(&(recv_Buf[0][0]), num_Buf, MPI_DOUBLE, 0, MPI_ANY_TAG, comm, &recv_Request);
        
        MPI_Wait(&recv_Request, &recv_Status);
        MPI_Wait(&send_Request, &send_Status);
        
        switch (recv_Status.MPI_TAG) {
            case TAG_ED: ED_Calculation();
                MPI_Isend(&(send_Buf[0][0]), num_Buf, MPI_DOUBLE, 0, TAG_ED, comm, &send_Request);
                break;
            case TAG_NB: NB_Calculation();
                MPI_Isend(&(send_Buf[0][0]), num_Buf, MPI_DOUBLE, 0, TAG_NB, comm, &send_Request);
                break;
            case TAG_SIGNAL: signal = 0; break;
            default: break;
        }
        
    }
    
}


