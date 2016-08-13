/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                  MSc in High Performance Computing, EPCC                     *
 *                       The University of Edinburgh                            *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  mpilib.cpp
 * Brief: The implementation of the MPI_Lib class functions
 */

#include "mpilib.hpp"


void MPI_Lib::create_MPI_Tetrad(MPI_Datatype* MPI_Tetrad, int num_Tetrads, Tetrad* tetrad) {
    
    int i, j, * counts = new int [5 * num_Tetrads];
    MPI_Datatype * old_Types = new MPI_Datatype [5 * num_Tetrads];
    MPI_Aint base,  * displs = new MPI_Aint [5 * num_Tetrads];
    
    for (i = 0; i < num_Tetrads; i++) {
        
        // The number of elements of each array
        counts[5 * i] = 3 * tetrad[i].num_Atoms;
        counts[5*i+1] = 3 * tetrad[i].num_Atoms;
        counts[5*i+2] = 3 * tetrad[i].num_Atoms;
        counts[5*i+3] = tetrad[i].num_Evecs;
        counts[5*i+4] = tetrad[i].num_Evecs * (3 * tetrad[i].num_Atoms);
        
        // The original data type of the arrays
        for (j = 0; j < 5; j++) { old_Types[5*i+j] = MPI_DOUBLE; }
        
        // Get the memory address of every elements in tetrad
        MPI_Get_address(&(tetrad[i].avg[0]),             &displs[5 * i]);
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
    
    delete [] counts;
    delete [] old_Types;
    delete [] displs;

}



void MPI_Lib::free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad) {
    
    MPI_Type_free(MPI_Tetrad);
    
}



void MPI_Lib::create_MPI_ED_Forces(MPI_Datatype* MPI_ED_Forces, Tetrad* tetrad) {
    
    int i, counts[2];
    MPI_Datatype old_Types[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint base, displs[2];
    
    counts[0] = 3 * tetrad->num_Atoms + 1; // ED
    counts[1] = 3 * tetrad->num_Atoms;     // random
    
    MPI_Get_address(tetrad, &base);
    MPI_Get_address(&(tetrad->ED_Forces[0]),    &displs[0]);
    MPI_Get_address(&(tetrad->random_Terms[0]), &displs[1]);

    displs[0] -= base;
    displs[1] -= base;

    MPI_Type_create_struct(2, counts, displs, old_Types, MPI_ED_Forces);
    MPI_Type_commit(MPI_ED_Forces);
    
}



void MPI_Lib::free_MPI_ED_Forces(MPI_Datatype* MPI_ED_Forces) {
    
    MPI_Type_free(MPI_ED_Forces);
    
}




void MPI_Lib::create_MPI_Crds(MPI_Datatype* MPI_Crds, int num_Tetrads, Tetrad* tetrad) {
    
    int i, * counts = new int [num_Tetrads];
    MPI_Datatype * old_Types = new MPI_Datatype [num_Tetrads];
    MPI_Aint base,  * displs = new MPI_Aint [num_Tetrads];
    
    for (i = 0; i < num_Tetrads; i++) {
        
        counts[i] = 3 * tetrad[i].num_Atoms;
        old_Types[i] = MPI_DOUBLE;
        MPI_Get_address(&(tetrad[i].coordinates[0]), &displs[i]);
        
    }
    
    MPI_Get_address(&(tetrad[0]), &base);
    for (i = num_Tetrads - 1; i >= 0; i--) { displs[i] -= base; }
    
    MPI_Type_create_struct(num_Tetrads, counts, displs, old_Types, MPI_Crds);
    MPI_Type_commit(MPI_Crds);
    
    delete [] counts;
    delete [] old_Types;
    delete [] displs;
    
}



void MPI_Lib::free_MPI_Crds(MPI_Datatype* MPI_Crds) {
    
    MPI_Type_free(MPI_Crds);
    
}





