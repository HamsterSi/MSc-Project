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
    
    array.deallocate_2D_Array(send_Buffer);
    array.deallocate_2D_Array(recv_Buffer);
    
    // Deallocate memory spaces for all arrays in tetrads
    for (int i = 0; i < num_Tetrads; i++) {
        array.deallocate_Tetrad_Arrays(&tetrad[i]);
    }
    
}



void Worker::recv_Parameters(void) {
    
    int i, parameters[2], signal = 1;
    double edmd_Para[8];
    
    // Receive parameters from the master process
    MPI_Recv(parameters, 2, MPI_INT, 0, TAG_DATA, comm, &status);
    num_Tetrads = parameters[0]; // The number of tetrads
    max_Atoms   = parameters[1]; // The maximum number of atoms in tetrads
    
    send_Buffer = array.allocate_2D_Array(2, 3 * max_Atoms + 2);
    recv_Buffer = array.allocate_2D_Array(2, 3 * max_Atoms + 2);
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
        
        tetrad[i].ED_Energy = tetrad[i].NB_Energy   = 0.0;
        tetrad[i].EL_Energy = tetrad[i].temperature = 0.0;
    }
    
    delete [] tetrad_Para;
    
    // Send feedback to master that has received all parameters.
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_DATA, comm);
    
}



void Worker::recv_Tetrads(void) {
    
    int i, signal = 1;
    MPI_Datatype MPI_Tetrad;
    
    // Receive all tetrads parameters from the master process
    for (i = 0; i < num_Tetrads; i++) {
        
        MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &tetrad[i]);
        MPI_Recv(&tetrad[i], 1, MPI_Tetrad, 0, TAG_TETRAD+i, comm, &status);
        MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);
        
    }
    
    // Send feedback to master that has received all tetrads
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_TETRAD, comm);
    
}



void Worker::ED_Calculation(MPI_Request* request) {
    
    // Assign values for tetrad indexes and coordinates
    int index = (int) recv_Buffer[0][3 * max_Atoms + 1];
    send_Buffer[0][3 * max_Atoms + 1] = recv_Buffer[0][3 * max_Atoms + 1];

    // Calculate ED forces (ED energy)
    edmd.calculate_ED_Forces(&tetrad[index], recv_Buffer[0], send_Buffer[0], send_Buffer[1], 3 * max_Atoms);
    
    // Send the calculated ED forces, ED energy index back
    MPI_Isend(&(send_Buffer[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, 0, TAG_ED, comm, request);
    
}



void Worker::NB_Calculation(MPI_Request* request) {
    
    // Assign values for tetrad indexes and coordinates
    int i = (int) recv_Buffer[0][3 * max_Atoms + 1];
    send_Buffer[0][3 * max_Atoms + 1] = recv_Buffer[0][3 * max_Atoms + 1];
    int j = (int) recv_Buffer[1][3 * max_Atoms + 1];
    send_Buffer[1][3 * max_Atoms + 1] = recv_Buffer[1][3 * max_Atoms + 1];
    
    // Calculate NB forces, NB energy & Electrostatic Energy
    edmd.calculate_NB_Forces(&tetrad[i], &tetrad[j], recv_Buffer, send_Buffer, 3 * max_Atoms);
    
    //cout << "worker: " << recv_Buffer[0][3 * max_Atoms + 1] << " " << recv_Buffer[1][3 * max_Atoms + 1] << " " << send_Buffer[0][3*max_Atoms] << " " << send_Buffer[0][3*max_Atoms] << endl;
    
    // Send NB forces, energies & indexes back to master
    MPI_Isend(&(send_Buffer[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, 0, TAG_NB, comm, request);
    
}



int Worker::terminate(void) {
    
    int signal;
    
    // Receive the terminate signal from master
    MPI_Recv(&signal, 1, MPI_INT, 0, TAG_SIGNAL, comm, &status);
    return signal;
}



