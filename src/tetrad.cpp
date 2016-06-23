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
 * File:  tetrad.cpp
 * Brief: Implementation of class functions for DNA tetrads
 */

#include "tetrad.hpp"

/*
Tetrad::Tetrad(void) {
    
    num_Atoms = 0;
    num_Evecs = 0;
    
}


Tetrad::~Tetrad(void) {
   
    delete []avg_Structure;
    delete []masses;
    delete []abq;
    delete []eigenvalues;
    
    for (int i = 0; i < num_Evecs; i++) {
        delete []eigenvectors[i];
    }
    delete []eigenvectors;
    
    delete []velocities;
    delete []coordinates;
    
    delete []ED_Forces;
    delete []random_Forces;
    delete []NB_Forces;
    
}*/

