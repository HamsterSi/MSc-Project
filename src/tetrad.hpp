/*
 * This the header file of the Tetrad class, it contains 
 * some basic charasteristics of tetrads.
 */

#ifndef tetrad_hpp
#define tetrad_hpp

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

/* Initial parameters of tetrads are read from "prm" file */
    double *NB_Forces;
class Tetrad {
    
public:
    
    int num_Atoms; // The number of atoms in tetrads
    
    int num_Evecs;           // The number of eigenvectors & eigenvalues
    
    double *avg_Structure;   // The average structure (coordinates) of tetrad
    
    double *masses;          // The masses of atoms of tetrads
    
    double *abq;             // The abq
    
    double *eigenvalues;     // The eigenvalues
    
    double **eigenvectors;   // The eigenvectors
    
    double *velocities;      // The velocities of tetrads
    
    double *coordinates;     // The coordinates of tetrads
    
    double *ED_Forces;
    
    double *random_Forces;
    
    double *NB_Forces;
    
    double ED_Energy;
    
    double NB_Energy;
    
    double EL_Energy;
    
    double temperature;    // The acerage temperature of tetrad
    
};



#endif /* tetrad_hpp */
