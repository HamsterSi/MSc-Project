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
#include "array.hpp"

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
    
    double *avg;           // The reference average structure of tetrad
    
    double *masses;        // The masses of every atom in tetrad
    
    double *abq;           // The non-bonded parameters
    
    double *eigenvalues;   // The eigenvalues
    
    double **eigenvectors; // The eigenvectors

    double *ED_Forces;     // The ED forces (Laset element is ED energy)
    
    double *random_Forces; // The random forces
    
    double *NB_Forces;     // The NB forces (Laset two elements are NB energy & Electrostatic Energy)
    
    double temperature;    // The temperature of tetrad
    
    double *velocities;    // The velocities of tetrad
    
    double *coordinates;   // The coordinates of tetrad
   
public:
    
    /**
     * Function:  Allocate memory space for all arrays in tetrad
     *
     * Parameter: None
     *
     * Return:    None
     */
    void allocate_Tetrad_Arrays(void);
    
    /**
     * Function:  Deallocate the memory space  of all arrays in tetrad
     *
     * Parameter: None
     *
     * Return:    None
     */
    void deallocate_Tetrad_Arrays(void);

};

#endif /* tetrad_hpp */
