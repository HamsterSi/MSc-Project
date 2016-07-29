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
 * File:  mpilib.cpp
 * Brief: Implementation of class functions for the MPI library for this code
 */

#include "mpilib.hpp"


void MPI_Library::create_MPI_Tetrad(MPI_Datatype* MPI_Tetrad, int num_Tetrads, Tetrad* tetrad) {
    
    int i, j, * counts = new int [5 * num_Tetrads];
    MPI_Datatype * old_Types = new MPI_Datatype [5 * num_Tetrads];
    MPI_Aint base,  * displs = new MPI_Aint [5 * num_Tetrads];
    
    for (i = 0; i < num_Tetrads; i++) {
        
        // The number of elements in every array
        counts[5 * i] = counts[5*i+1] = counts[5*i+2] = 3 * tetrad[i].num_Atoms;
        counts[5*i+3] = tetrad[i].num_Evecs;
        counts[5*i+4] = tetrad[i].num_Evecs * (3 * tetrad[i].num_Atoms);
        
        // The old data type of arrays
        for (j = 0; j < 5; j++) { old_Types[5*i+j] = MPI_DOUBLE; }
        
        // Get the memory address of every elements in tetrad
        MPI_Get_address(&(tetrad[i].avg_Structure[0]),   &displs[5 * i]);
        MPI_Get_address(&(tetrad[i].masses[0]),          &displs[5*i+1]);
        MPI_Get_address(&(tetrad[i].abq[0]),             &displs[5*i+2]);
        MPI_Get_address(&(tetrad[i].eigenvalues[0]),     &displs[5*i+3]);
        MPI_Get_address(&(tetrad[i].eigenvectors[0][0]), &displs[5*i+4]);
        
    }
    
    // Calculate the displacements
    MPI_Get_address(&(tetrad[0]), &base);
    for (i = 5 * num_Tetrads - 1; i >= 0; i--) { displs[i] -= base; }
    
    // Create the MPI data type "MPI_Tetrad".
    MPI_Type_create_struct(5 * num_Tetrads, counts, displs, old_Types, MPI_Tetrad);
    MPI_Type_commit(MPI_Tetrad);

}



void MPI_Library::free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad) {
    
    MPI_Type_free(MPI_Tetrad);
    
}
