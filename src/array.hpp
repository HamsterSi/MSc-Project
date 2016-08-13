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
 * File:  array.hpp
 * Brief: The declaration of a Array class for 2D array allocation and deallocation
 */

#ifndef arrays_hpp
#define arrays_hpp

#include <iostream>

#include "tetrad.hpp"

/**
 * Brief: The Array class for 2D (double and integer) array allocation and deallocation.
 */
class Array{
    
public:
    
    /**
     * Function:  Create a 2D double array within a continguous memory space
     *
     * Parameter: int rows -> The number of rows
     *            int cols -> The number of columns
     *
     * Return:    The 2D double array
     */
    static double** allocate_2D_Double_Array(int rows, int cols);
    
    /**
     * Function:  free the memory space of the 2D double array
     *
     * Parameter: double** array -> The 2D double array to be freed
     *
     * Return:    None
     */
    static void deallocate_2D_Double_Array(double** array);
    
    /**
     * Function:  Create a 2D int array within continguous memory space
     *
     * Parameter: int rows -> The number of rows
     *            int cols -> The number of columns
     *
     * Return:    A 2D int array
     */
    static int** allocate_2D_Int_Array(int rows, int cols);
    
    /**
     * Function:  free the memory space of 2D array
     *
     * Parameter: int** array -> The 2D int array to be freed
     *
     * Return:    None
     */
    static void deallocate_2D_Int_Array(int** array);
    
};

#endif /* arrays_hpp */
