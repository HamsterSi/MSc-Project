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


void MPI_Library::create_MPI_Tetrad(MPI_Datatype* MPI_Tetrad, Tetrad* tetrad) {
  
    
    int i, counts[5], num_Atoms = 3 * tetrad->num_Atoms;
    MPI_Datatype old_Types[5];
    MPI_Aint  base, displs[5];
    
    // Assign values for counts & old types
    for (i = 0; i < 5; i++) {
        if (i < 3) counts[i] = num_Atoms;
        else if (i == 3) counts[i] = tetrad->num_Evecs;
        else counts[i] = tetrad->num_Evecs * num_Atoms;
        
        old_Types[i] = MPI_DOUBLE;
    }
    
    // Get the memory address of every elements in tetrad
    MPI_Get_address(tetrad, &base);
    MPI_Get_address(tetrad->avg_Structure, &displs[0]);
    MPI_Get_address(tetrad->masses,        &displs[1]);
    MPI_Get_address(tetrad->abq,           &displs[2]);
    MPI_Get_address(tetrad->eigenvalues,   &displs[3]);
    MPI_Get_address(&(tetrad->eigenvectors[0][0]),  &displs[4]);
    
    // Calculate the displacements
    for (i = 4; i >= 0; i--) { displs[i] -= base; }
    
    // Create the MPI data type "MPI_Tetrad".
    MPI_Type_create_struct(5, counts, displs, old_Types, MPI_Tetrad);
    MPI_Type_commit(MPI_Tetrad);

}



void MPI_Library::free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad) {
    
    MPI_Type_free(MPI_Tetrad);
    
}
