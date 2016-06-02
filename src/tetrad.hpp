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

// It still needs to be discussed whether it is better to use structure or class for tetrads
typedef struct _Tetrad {
    
    int num_Atoms;        // The number of atoms in tetrads
    
    int num_Evecs;        // The number of eigenvectors & eigenvalues
    
    float *avg_Structure; // The average structure (coordinates) of tetrad
    
    float *masses;        // The masses of atoms of tetrads
    
    float *abq;           // The abq
    
    float *eigenvalues;   // The eigenvalues
    
    float **eigenvectors; // The eigenvectors
    
    float *velocities;    // The velocities of tetrads

    float *coordinates;   // The coordinates of tetrads
    
    float *ED_Forces;
    
    float *random_Forces;
    
    float *NB_Forces;
    
    float energies[3];    // ED_Energy, NB_Energy, Electrostatic_Energy
    
    float temperature;    // The acerage temperature of tetrad
    
}Tetrad;

/*
class Tetrad {
    
public:
    
    int num_Atoms; // The number of atoms in tetrads
    
    int num_Evecs;           // The number of eigenvectors & eigenvalues
    
    float *avg_Structure;   // The average structure (coordinates) of tetrad
    
    float *masses;          // The masses of atoms of tetrads
    
    float *abq;             // The abq
    
    float *eigenvalues;     // The eigenvalues
    
    float **eigenvectors;   // The eigenvectors
    
    float *velocities;      // The velocities of tetrads
    
    float *coordinates;     // The coordinates of tetrads
    
    float *ED_Forces;
    
    float *random_Forces;
    
    float *NB_Forces;
    
public:
    
    Tetrad(void);
    
    ~Tetrad(void);
};*/



#endif /* tetrad_hpp */
