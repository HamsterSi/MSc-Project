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



void Array::assignment(int num, double* source, double* dest) {
    
    for (int i = 0; i < num; i++) {
        dest[i] = source[i];
    }
    
}



