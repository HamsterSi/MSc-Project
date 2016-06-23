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
    
    int i, num_Atoms = 3 * tetrad->num_Atoms;
    
    // Allocate the arrays for counts, old types & displacements of parameters of tetrads
    int * counts = new int [8 + tetrad->num_Evecs];
    MPI_Datatype * old_Types = new MPI_Datatype [8 + tetrad->num_Evecs];
    MPI_Aint base, * displs = new MPI_Aint [8 + tetrad->num_Evecs];
    
    // Assign values for counts & old types
    counts[0] = counts[1] = 1;  old_Types[0] = old_Types[1] = MPI_INT;
    for (i = 0; i < 6 + tetrad->num_Evecs; i++) {
        if (i + 2 == 5) counts[i + 2] = tetrad->num_Evecs;
        else counts[i + 2] = num_Atoms;
        
        old_Types[i + 2] = MPI_DOUBLE;
    }
    
    // Get the memory address of every elements in tetrad
    MPI_Get_address(tetrad, &base);
    MPI_Get_address(&(tetrad->num_Atoms),  &displs[0]);
    MPI_Get_address(&(tetrad->num_Evecs),  &displs[1]);
    MPI_Get_address(tetrad->avg_Structure, &displs[2]);
    MPI_Get_address(tetrad->masses,        &displs[3]);
    MPI_Get_address(tetrad->abq,           &displs[4]);
    MPI_Get_address(tetrad->eigenvalues,   &displs[5]);
    for (i = 0; i < tetrad->num_Evecs; i++) {
        MPI_Get_address(tetrad->eigenvectors[i], &displs[6 + i]);
    }
    MPI_Get_address(tetrad->coordinates,   &displs[8 + tetrad->num_Evecs - 2]);
    MPI_Get_address(tetrad->velocities,    &displs[8 + tetrad->num_Evecs - 1]);
    
    // Calculate the displacements
    for (i = 8 + tetrad->num_Evecs - 1; i >= 0; i--) {
        displs[i] -= base;
    }
    
    // Create the MPI data type "MPI_Tetrad".
    MPI_Type_create_struct(8 + tetrad->num_Evecs, counts, displs, old_Types, MPI_Tetrad);
    MPI_Type_commit(MPI_Tetrad);

}



void MPI_Library::free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad) {
    
    // Free the MPI data type
    MPI_Type_free(MPI_Tetrad);
}
