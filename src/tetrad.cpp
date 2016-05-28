//
//  tetrad.cpp
//  
//
//  Created by Zhuowei Si on 11/04/2016.
//
//

#include "tetrad.hpp"


Tetrad::Tetrad(void) {
    
    num_Atoms_In_Tetrad = 0;
    num_Evecs = 0;
    
}


Tetrad::~Tetrad(void) {
    
    /*
    delete []avg_Structure;
    delete []masses;
    delete []abq;
    delete []eigenvalues;
    
    for (int i = 0; i < num_Evecs; i++) {
        delete []eigenvectors[i];
    }
    delete []eigenvectors;
    
    delete []coordinates;
    delete []velocities;
     */
}