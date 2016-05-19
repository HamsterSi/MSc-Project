//
//  projectTest.cpp
//  
//
//  Created by Zhuowei Si on 04/04/2016.
//
//

#include "mpilibrary.hpp"

void MPI_Library::create_MPI_Tetrad(MPI_Datatype* MPI_Tetrad, Tetrad* tetrad) {
    int num_Atoms = 3 * tetrad->num_Atoms_In_Tetrad;
    
    // For 2D array "eigenvectors"
    MPI_Datatype type_2D;
    MPI_Type_contiguous(num_Atoms, MPI_FLOAT, &type_2D);
    
    //MPI_Datatype MPI_Tetrad;
    /*
    int counts[9] = {1, 1, num_Atoms, num_Atoms, num_Atoms, tetrad->num_Evecs, tetrad->num_Evecs, num_Atoms, num_Atoms};
    MPI_Datatype types[9] = {MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, type_2D, MPI_FLOAT, MPI_FLOAT};
    
    MPI_Aint base, displs[9];
    
    MPI_Address(tetrad, &base);
    MPI_Address(&(tetrad->num_Atoms_In_Tetrad), &displs[0]);
    MPI_Address(&(tetrad->num_Evecs), &displs[1]);
    MPI_Address(&(tetrad->avg_Structure), &displs[2]);
    MPI_Address(&(tetrad->masses), &displs[3]);
    MPI_Address(&(tetrad->abq), &displs[4]);
    MPI_Address(&(tetrad->eigenvalues), &displs[5]);
    MPI_Address(&(tetrad->eigenvectors), &displs[6]);
    MPI_Address(&(tetrad->coordinates), &displs[7]);
    MPI_Address(&(tetrad->velocities), &displs[8]);
    for (int i = 8; i >= 0; i--) {
        displs[i] -= base;
    }
    
    MPI_Aint displs[9] = {0, sizeof(int), 2*sizeof(int),
        2*sizeof(int) +   num_Atoms*sizeof(float),
        2*sizeof(int) + 2*num_Atoms*sizeof(float),
        2*sizeof(int) + 3*num_Atoms*sizeof(float),
        2*sizeof(int) + 3*num_Atoms*sizeof(float) + tetrad->num_Evecs*sizeof(float),
        2*sizeof(int) + 3*num_Atoms*sizeof(float) + tetrad->num_Evecs*sizeof(float) + tetrad->num_Evecs*num_Atoms*sizeof(float),
        2*sizeof(int) + 4*num_Atoms*sizeof(float) + tetrad->num_Evecs*sizeof(float) + tetrad->num_Evecs*num_Atoms*sizeof(float),
    };
    
    MPI_Aint displs[9] = { offsetof(Tetrad, num_Atoms_In_Tetrad),
        offsetof(Tetrad, num_Evecs),    offsetof(Tetrad, avg_Structure),
        offsetof(Tetrad, masses),       offsetof(Tetrad, abq),
        offsetof(Tetrad, eigenvalues),  offsetof(Tetrad, eigenvectors),
        offsetof(Tetrad, coordinates),  offsetof(Tetrad, velocities)
    };
    
    MPI_Type_create_struct(9, counts, displs, types, MPI_Tetrad); */
    
    int counts[6] = {1, 1, num_Atoms, num_Atoms, num_Atoms, tetrad->num_Evecs};
    MPI_Datatype types[6] = {MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
}

void MPI_Library::free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad) {
    MPI_Type_free(MPI_Tetrad);
}
