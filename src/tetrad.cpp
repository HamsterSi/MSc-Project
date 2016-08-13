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
 * File:  tetrad.cpp
 * Brief: The implementation of the Tetrad class functions
 */

#include "tetrad.hpp"


void Tetrad::allocate_Tetrad_Arrays(void) {
    
    avg          = new double[3 * num_Atoms];
    masses       = new double[3 * num_Atoms];
    abq          = new double[3 * num_Atoms];
    eigenvalues  = new double[num_Evecs];
    eigenvectors = Array::allocate_2D_Double_Array(num_Evecs, 3 * num_Atoms);
    velocities   = new double[3 * num_Atoms];
    coordinates  = new double[3 * num_Atoms];
    ED_Forces    = new double[3 * num_Atoms + 1];
    random_Terms = new double[3 * num_Atoms];
    NB_Forces    = new double[3 * num_Atoms + 2];
    
}



void Tetrad::deallocate_Tetrad_Arrays(void) {
    
    delete [] avg;
    delete [] masses;
    delete [] abq;
    delete [] eigenvalues;
    Array::deallocate_2D_Double_Array(eigenvectors);
    delete [] velocities;
    delete [] coordinates;
    delete [] ED_Forces;
    delete [] random_Terms;
    delete [] NB_Forces;
    
}




