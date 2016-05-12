//
//  projectTest.cpp
//  
//
//  Created by Zhuowei Si on 04/04/2016.
//
//

#include "mpilibrary.hpp"

/* Initialise the MPI library */
void MPI_Library::initialize(int* argc, char **argv[]){
    MPI_Init(argc, argv);
}

int MPI_Library::get_Size(int* size){
    
    MPI_Comm_size(MPI_COMM_WORLD, size);
    return *size;
}

int MPI_Library::get_Rank(int* rank){
    
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    return *rank;
}

void MPI_Library::create_MPI_Tetrad(MPI_Datatype MPI_Tetrad, int num_Atoms_In_Tetrad, int num_Evecs) {
    int num_Atoms = 3 * num_Atoms_In_Tetrad;
    
    // For 2D array "eigenvectors"
    MPI_Datatype type_2D;
    MPI_Type_contiguous(num_Atoms, MPI_FLOAT, &type_2D);
    
    //MPI_Datatype MPI_Tetrad;
    int counts[9] = {1, 1, num_Atoms, num_Atoms, num_Atoms, num_Evecs, num_Evecs, num_Atoms, num_Atoms};
    MPI_Datatype types[9] = {MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, type_2D, MPI_FLOAT, MPI_FLOAT};
    MPI_Aint displs[9] = { offsetof(Tetrad, num_Atoms_In_Tetrad),
        offsetof(Tetrad, num_Evecs),    offsetof(Tetrad, avg_Structure),
        offsetof(Tetrad, masses),       offsetof(Tetrad, abq),
        offsetof(Tetrad, eigenvalues),  offsetof(Tetrad, eigenvectors),
        offsetof(Tetrad, coordinates),  offsetof(Tetrad, velocities) };
    
    MPI_Type_create_struct(9, counts, displs, types, &MPI_Tetrad);
    MPI_Type_commit(&MPI_Tetrad);
    
    MPI_Type_free(&type_2D);
}

void MPI_Library::free_MPI_Tetrad(MPI_Datatype MPI_Tetrad) {
    MPI_Type_free(&MPI_Tetrad);
}

/* Finalize the MPI library */
void MPI_Library::finalize(void) {
    MPI_Finalize();
}