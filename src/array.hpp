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

#include <stdio.h>

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
    static void allocate_2D_Array(double** array);
    
};

#endif /* arrays_hpp */
