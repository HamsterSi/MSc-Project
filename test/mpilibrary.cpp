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

/* Finalize the MPI library */
void MPI_Library::finalize(void) {
    MPI_Finalize();
}