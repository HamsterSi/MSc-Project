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
 * File:  array.hpp
 * Brief: Declaration of functions for 2D array allocation and deallocation
 */

#ifndef arrays_hpp
#define arrays_hpp

#include <iostream>

#include "tetrad.hpp"

/**
 * Brief: A class for 2D array allocation and deallocation.
 */
class Array{
    
public:
    
    /**
     * Function:  Create a 2D array with continguous memory space
     *
     * Parameter: int rows -> The number of rows
     *            int cols -> The number of columns
     *
     * Return:    A 2D array
     */
    static double** allocate_2D_Array(int rows, int cols);
    
    /**
     * Function:  free the memory space of 2D array
     *
     * Parameter: double** array -> The 2D array to be freed
     *
     * Return:    None
     */
    static void deallocate_2D_Array(double** array);
    
    /**
     * Function:  Allocate all arrays in tetrad
     *
     * Parameter: Tetrad* tetrad -> The tetrad whose arrays to be allocated
     *
     * Return:    None
     */
    static void allocate_Tetrad_Arrays(Tetrad* tetrad);
    
    /**
     * Function:  Deallocate all arrays in tetrad
     *
     * Parameter: Tetrad* tetrad -> The tetrad whose arrays to be deallocated
     *
     * Return:    None
     */
    static void deallocate_Tetrad_Arrays(Tetrad* tetrad);
    
    /**
     * Function:  Assign values from source array to dest array (for 1D array)
     *
     * Parameter: int num        -> The number of elements to be assigned
     *            double* source -> The source array
     *            double* dest   -> The dest array
     *
     * Return:    None
     */
    static void assignment(int num, double* source, double* dest);
    
};

#endif /* arrays_hpp */
