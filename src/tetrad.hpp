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
 * File:  tetrad.hpp
 * Brief: The declaration of the Tetrad class
 */

#ifndef tetrad_hpp
#define tetrad_hpp

#include <iostream>
#include "array.hpp"

using namespace std;

/**
 * Brief: The Tetrad class that contains all the essential parameters and varialbes
 *        of tetrads for the ED/MD simulation.
 */
class Tetrad {
    
public:
    
    int num_Atoms;         // The number of atoms in tetrads
    
    int num_Evecs;         // The number of eigenvectors & eigenvalues
    
    double * avg;          // The reference average structure
    
    double * masses;       // The masses of every atom in tetrad
    
    double * abq;          // The non-bonded parameters
    
    double * eigenvalues;  // The eigenvalues calculated from PCA
    
    double** eigenvectors; // The eigenvectors obtained from PCA

    double * ED_Forces;    // The ED forces (Laset element: ED energy)
    
    double * random_Terms; // The random terms for Langevin dynamics
    
    double * NB_Forces;    // The NB forces (Laset two elements: NB energy & Electrostatic Energy)
    
    double temperature;    // The temperature of tetrad
    
    double * velocities;   // The velocities of tetrad
    
    double * coordinates;  // The coordinates of tetrad
   
public:
    
    /**
     * Function:  Allocate memory space for all the arrays in tetrad
     *
     * Parameter: None
     *
     * Return:    None
     */
    void allocate_Tetrad_Arrays(void);
    
    /**
     * Function:  Deallocate the memory space of the arrays in tetrad
     *
     * Parameter: None
     *
     * Return:    None
     */
    void deallocate_Tetrad_Arrays(void);

};

#endif /* tetrad_hpp */
