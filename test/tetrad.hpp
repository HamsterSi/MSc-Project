//
//  tetrad.hpp
//  
//
//  Created by Zhuowei Si on 11/04/2016.
//
//

#ifndef tetrad_hpp
#define tetrad_hpp

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "mpi.h"

/* Parameters pf tetrads form prm file */
typedef struct _Tetrad {
    int num_Atoms_In_Tetrad; // The number of atoms in tetrads
    int num_Evecs;           // The number of eigenvectors & eigenvalues
    float *avg_Structure;    // The average structure (coordinates) of tetrad
    float *masses;           // The masses of atoms of tetrads
    float *abq;              // The abq
    float *eigenvalues;      // The eigenvalues
    float **eigenvectors;    // The eigenvectors
    float *coordinates;      // The coordinates of tetrads
    float *velocities;       // The velocities of tetrads
    
    void calculate_ED_Forces(void);
    
    void calculate_VDW_Forces(void);
    
    void update_Velocities(void);
    
    void update_Coordinates(void);
}Tetrad;



#endif /* tetrad_hpp */
