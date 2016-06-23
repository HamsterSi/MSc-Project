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
 * File:  tetrad.hpp
 * Brief: Declaration of a class for DNA tetrads
 */

#ifndef tetrad_hpp
#define tetrad_hpp

#include <iostream>

using namespace std;

/**
 * Brief: A class which contains all essential parameters for tetrads
 *
 * The "Tetrad" is the basic class in the code, whose parameters are either obtained
 * from input files or calculated from simualtion.
 */
class Tetrad {
    
public:
    
    int num_Atoms;         // The number of atoms in tetrads
    
    int num_Evecs;         // The number of eigenvectors & eigenvalues
    
    double *avg_Structure; // The average structure of tetrad
    
    double *masses;        // The masses of atoms in tetrad
    
    double *abq;           // The abq
    
    double *eigenvalues;   // The eigenvalues
    
    double **eigenvectors; // The eigenvectors
    
    double *velocities;    // The velocities of tetrad
    
    double *coordinates;   // The coordinates of tetrad
    
    double *ED_Forces;     // The ED forces
    
    double *random_Forces; // The random forces
    
    double *NB_Forces;     // The NB forces
    
    double ED_Energy;      // The ED energy
    
    double NB_Energy;      // The NB energy
    
    double EL_Energy;      // The electrostatic energy
    
    double temperature;    // The temperature of tetrad
    
};

#endif /* tetrad_hpp */
