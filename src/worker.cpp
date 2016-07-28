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
    
    int i, parameters[2], signal = 1;
    double edmd_Para[8];
    MPI_Status status;
    
    // Receive parameters from the master process
    MPI_Recv(parameters, 2, MPI_INT, 0, TAG_DATA, comm, &status);
    num_Tetrads = parameters[0]; // The number of tetrads
    max_Atoms   = parameters[1]; // The maximum number of atoms in tetrads
    
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
        
        tetrad[i].allocate_Tetrad_Arrays();
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




void Worker::ED_Calculation(int index[], int num_Buf, double** recv_Buf, MPI_Request send_Rqt[]) {
    
    // Calculate ED forces & energy
    index[0] = (int) recv_Buf[0][num_Buf];
    index[1] = 0;
    index[2] = TAG_ED;
    edmd.calculate_ED_Forces(&tetrad[index[0]], recv_Buf[0], num_Buf);
    
    // send ED forces, energy back
    MPI_Isend(index, 3, MPI_INT, 0, TAG_INDEX, comm, &send_Rqt[0]);
    MPI_Isend(tetrad[index[0]].ED_Forces,  num_Buf + 1, MPI_DOUBLE, 0, TAG_ED + 1, comm, &send_Rqt[1]);
    MPI_Isend(tetrad[index[0]].coordinates, num_Buf,    MPI_DOUBLE, 0, TAG_ED + 2, comm, &send_Rqt[2]);
    
}



void Worker::NB_Calculation(int index[], int num_Buf, double** recv_Buf, MPI_Request send_Rqt[]) {
    
    // Calculate NB forces, NB energy & Electrostatic Energy
    index[0] = (int) recv_Buf[0][num_Buf];
    index[1] = (int) recv_Buf[1][num_Buf];
    index[2] = TAG_NB;
    edmd.calculate_NB_Forces(&tetrad[index[0]], &tetrad[index[1]], recv_Buf, num_Buf);
    
    // send NB forces, energy back
    MPI_Isend(index, 3, MPI_INT, 0, TAG_INDEX, comm, &send_Rqt[0]);
    MPI_Isend(tetrad[index[0]].NB_Forces, num_Buf + 2, MPI_DOUBLE, 0, TAG_NB + 1, comm, &send_Rqt[1]);
    MPI_Isend(tetrad[index[1]].NB_Forces, num_Buf + 2, MPI_DOUBLE, 0, TAG_NB + 2, comm, &send_Rqt[2]);
    
}



void Worker::force_Calculation(void) {
    
    int signal = 1, index[3], num_Buf = 3 * max_Atoms;
    double ** recv_Buf = Array::allocate_2D_Array(2, num_Buf + 1);
    
    MPI_Request send_Rqt[3], recv_Rqt;
    MPI_Status send_Status, recv_Status;
    
    // Receive the first tetrad index(es) & coordinates from master for ED/NB calculation
    MPI_Recv(&(recv_Buf[0][0]), 2 * (num_Buf + 1), MPI_DOUBLE, 0, MPI_ANY_TAG, comm, &recv_Status);
    
    switch (recv_Status.MPI_TAG) {
        case TAG_ED: ED_Calculation(index, num_Buf, recv_Buf, send_Rqt); break;
        case TAG_NB: NB_Calculation(index, num_Buf, recv_Buf, send_Rqt); break;
        default: break;
    }
    
    // The loop to take new jobs & send results back to the master
    while (signal == 1) {
        
        MPI_Irecv(&(recv_Buf[0][0]), 2 * (num_Buf + 1), MPI_DOUBLE, 0, MPI_ANY_TAG, comm, &recv_Rqt);
        
        MPI_Wait(&recv_Rqt,    &recv_Status);
        MPI_Wait(&send_Rqt[0], &send_Status);
        MPI_Wait(&send_Rqt[1], &send_Status);
        MPI_Wait(&send_Rqt[2], &send_Status);
        
        switch (recv_Status.MPI_TAG) {
            case TAG_ED: ED_Calculation(index, num_Buf, recv_Buf, send_Rqt); break;
            case TAG_NB: NB_Calculation(index, num_Buf, recv_Buf, send_Rqt); break;
            case TAG_SIGNAL: signal = 0; break;
            default: break;
        }
        
        //cout << " W: " << index[0] << " " << index[1] << " " << index[2] << endl;
    
    }

    Array::deallocate_2D_Array(recv_Buf);
    
}


