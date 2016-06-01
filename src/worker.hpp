
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
    
    int rank;
    
    int num_Tetrads;
    
    int max_Atoms;
    
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
