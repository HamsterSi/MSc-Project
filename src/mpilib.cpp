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
    
    int i, j, * counts = new int [6 * num_Tetrads];
    MPI_Datatype * old_Types = new MPI_Datatype [6 * num_Tetrads];
    MPI_Aint base,  * displs = new MPI_Aint [6 * num_Tetrads];
    
    for (i = 0; i < num_Tetrads; i++) {
        
        // The number of elements in every array
        counts[6 * i] = counts[6*i+1] = 3 * tetrad[i].num_Atoms;
        counts[6*i+2] = counts[6*i+5] = 3 * tetrad[i].num_Atoms;
        counts[6*i+3] = tetrad[i].num_Evecs;
        counts[6*i+4] = tetrad[i].num_Evecs * (3 * tetrad[i].num_Atoms);
        
        // The old data type of arrays
        for (j = 0; j < 6; j++) { old_Types[6*i+j] = MPI_DOUBLE; }
        
        // Get the memory address of every elements in tetrad
        MPI_Get_address(&(tetrad[i].avg[0]),             &displs[6 * i]);
        MPI_Get_address(&(tetrad[i].masses[0]),          &displs[6*i+1]);
        MPI_Get_address(&(tetrad[i].abq[0]),             &displs[6*i+2]);
        MPI_Get_address(&(tetrad[i].eigenvalues[0]),     &displs[6*i+3]);
        MPI_Get_address(&(tetrad[i].eigenvectors[0][0]), &displs[6*i+4]);
        MPI_Get_address(&(tetrad[i].coordinates[0]),     &displs[6*i+5]);
        
    }
    
    // Calculate the displacements
    MPI_Get_address(&(tetrad[0]), &base);
    for (i = 6 * num_Tetrads - 1; i >= 0; i--) { displs[i] -= base; }
    
    // Create the MPI data type "MPI_Tetrad".
    MPI_Type_create_struct(6 * num_Tetrads, counts, displs, old_Types, MPI_Tetrad);
    MPI_Type_commit(MPI_Tetrad);
    
    delete [] counts;
    delete [] old_Types;
    delete [] displs;

}



void MPI_Library::free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad) {
    
    MPI_Type_free(MPI_Tetrad);
    
}



void MPI_Library::create_MPI_Force(MPI_Datatype* MPI_Force, int num_Tetrads, Tetrad* tetrad) {
    
    int i, j, * counts = new int [4 * num_Tetrads];
    MPI_Datatype * old_Types = new MPI_Datatype [4 * num_Tetrads];
    MPI_Aint base,  * displs = new MPI_Aint [4 * num_Tetrads];
    
    for (i = 0; i < num_Tetrads; i++) {
        
        // The number of elements in every array
        counts[4 * i] = 3 * tetrad[i].num_Atoms + 1; // ED
        counts[4*i+1] = 3 * tetrad[i].num_Atoms;     // Crds
        counts[4*i+2] = 3 * tetrad[i].num_Atoms;     // random
        counts[4*i+3] = 3 * tetrad[i].num_Atoms + 2; // NB
        
        // The old data type of NB forces
        for (j = 0; j < 4; j++) { old_Types[4*i+j] = MPI_DOUBLE; }
        
        // Get the memory address of every elements in tetrad
        MPI_Get_address(&(tetrad[i].ED_Forces[0]),     &displs[4 * i]);
        MPI_Get_address(&(tetrad[i].coordinates[0]),   &displs[4*i+1]);
        MPI_Get_address(&(tetrad[i].random_Forces[0]), &displs[4*i+2]);
        MPI_Get_address(&(tetrad[i].NB_Forces[0]),     &displs[4*i+3]);
        
    }
    
    // Calculate the displacements
    MPI_Get_address(&(tetrad[0]), &base);
    for (i = 4 * num_Tetrads - 1; i >= 0; i--) { displs[i] -= base; }
    
    // Create the MPI data type "MPI_Tetrad".
    MPI_Type_create_struct(4 * num_Tetrads, counts, displs, old_Types, MPI_Force);
    MPI_Type_commit(MPI_Force);
    
    delete [] counts;
    delete [] old_Types;
    delete [] displs;
    
}



void MPI_Library::free_MPI_Force(MPI_Datatype* MPI_Force) {
    
    MPI_Type_free(MPI_Force);
    
}




