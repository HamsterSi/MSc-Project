
#ifndef worker_hpp
#define worker_hpp

#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "mpi.h"

#include "mpilibrary.hpp"
#include "edmd.hpp"
#include "tetrad.hpp"
#include "io.hpp"

using namespace std;

class Worker_Management {
    
public:
    
    int parameters[2]; // num_Tetrads, max_Atoms_In_Tetrad
    
    int *num_Atoms_N_Evecs;
    
    Tetrad *tetrad;
    
    EDMD edmd;
    
    MPI_Status status;
    
    MPI_Comm comm;
    
public:
    
    Worker_Management(void);
    
    ~Worker_Management(void);
    
    void parameters_Receiving(void);
    
    void tetrads_Receiving(void);
    
    void ED_Calculation(void);
    
    void NB_Calculation(void);
    
};


#endif /* worker_hpp */
