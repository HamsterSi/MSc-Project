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
 * File:  array.cpp
 * Brief: Implementation of functions for 2D array allocation and deallocation
 */

#include "array.hpp"


double** Array::allocate_2D_Array(int rows, int cols) {
    
    double ** array = new double * [rows];
    double * sub_Array = new double [rows * cols];
    
    for (int i = 0; i < rows; i++) {
        array[i] = sub_Array; sub_Array += cols;
    }
    
    return array;
}



void Array::deallocate_2D_Array(double** array) {
    
    delete [] array[0];
    delete [] array;
    
}



void Array::allocate_Tetrad_Arrays(Tetrad* tetrad) {
    
    tetrad->avg_Structure = new double[3 * tetrad->num_Atoms];
    tetrad->masses        = new double[3 * tetrad->num_Atoms];
    tetrad->abq           = new double[3 * tetrad->num_Atoms];
    tetrad->eigenvalues   = new double[tetrad->num_Evecs];
    tetrad->eigenvectors  = allocate_2D_Array(tetrad->num_Evecs, 3 * tetrad->num_Atoms);
    tetrad->velocities    = new double[3 * tetrad->num_Atoms];
    tetrad->coordinates   = new double[3 * tetrad->num_Atoms];
    tetrad->ED_Forces     = new double[3 * tetrad->num_Atoms];
    tetrad->random_Forces = new double[3 * tetrad->num_Atoms];
    tetrad->NB_Forces     = new double[3 * tetrad->num_Atoms];
    
}



void Array::deallocate_Tetrad_Arrays(Tetrad* tetrad) {
    
    delete [] tetrad->avg_Structure;
    delete [] tetrad->masses;
    delete [] tetrad->abq;
    delete [] tetrad->eigenvalues;
    deallocate_2D_Array(tetrad->eigenvectors);
    delete [] tetrad->velocities;
    delete [] tetrad->coordinates;
    delete [] tetrad->ED_Forces;
    delete [] tetrad->random_Forces;
    delete [] tetrad->NB_Forces;
    
}




