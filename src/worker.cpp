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
    
    array.deallocate_2D_Array(buffer);
    
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
    
    buffer = array.allocate_2D_Array(2, 3 * max_Atoms + 2);
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



void Worker::ED_Calculation(void) {
    
    // Receive the tetrad index from the master process
    //MPI_Recv(&(buffer[0][0]), 3 * max_Atoms + 2, MPI_DOUBLE, 0, TAG_ED, comm, &status);
    
    // Assign values for tetrad indexes and coordinates
    int index = (int) buffer[0][3 * max_Atoms + 1];

    // Calculate ED forces (ED energy)
    edmd.calculate_ED_Forces(&tetrad[index], buffer[0]);
    
    // Assign ED forces & random Forces to the 2D array for sending once
    array.assignment(3 * tetrad[index].num_Atoms, tetrad[index].ED_Forces  , buffer[0]);
    array.assignment(3 * tetrad[index].num_Atoms, tetrad[index].coordinates, buffer[1]);
    
    // Need to send the ED Energy back
    buffer[0][3 * max_Atoms] = tetrad[index].ED_Energy;
    
    // Send the calculated ED forces, ED energy, random forces & index back
    MPI_Send(&(buffer[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, 0, TAG_ED, comm);
    
}



void Worker::NB_Calculation(void) {

    // Receive tetrad indexes for NB forces calculation
    //MPI_Recv(&(buffer[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, 0, TAG_NB, comm, &status);
    
    // Assign values for tetrad indexes and coordinates
    int idx1 = (int) buffer[0][3 * max_Atoms + 1];
    int idx2 = (int) buffer[1][3 * max_Atoms + 1];
    
    // Calculate NB forces, NB energy & Electrostatic Energy
    edmd.calculate_NB_Forces(&tetrad[idx1], &tetrad[idx2], buffer[0], buffer[1]);
    
    // Assign NB forces to the 2D array to send back once
    array.assignment(3 * tetrad[idx1].num_Atoms, tetrad[idx1].NB_Forces, buffer[0]);
    array.assignment(3 * tetrad[idx2].num_Atoms, tetrad[idx2].NB_Forces, buffer[1]);
    
    // Need to send NB Energy & Electrostatic Energy back (bacause the energies in both tetrads are the same when calculation, so only need to send one set)
    buffer[0][3 * max_Atoms] = tetrad[idx1].NB_Energy;
    buffer[1][3 * max_Atoms] = tetrad[idx1].EL_Energy;
    
    // Send NB forces, energies & indexes back to master
    MPI_Send(&(buffer[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, 0, TAG_NB, comm);
    
}



int Worker::terminate(void) {
    
    int signal;
    
    // Receive the terminate signal from master
    MPI_Recv(&signal, 1, MPI_INT, 0, TAG_SIGNAL, comm, &status);
    return signal;
}



