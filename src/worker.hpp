
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
    
    int rank;           // The rank of worker process
    
    int num_Tetrads;    // The number of total tetrads
    
    int max_Atoms;      // The maximum number of atoms in tetrads
    
    Tetrad *tetrad;     // Tetrad array, used to stroe tetrads
    
    EDMD edmd;          // EDMD class, needs to call functions to calculate forces
    
    MPI_Status status;  // MPI Status
    
    MPI_Comm comm;      // MPI Communicator
    
public:
    
    Worker_Management(void);
    
    ~Worker_Management(void);
    
    void parameters_Receiving(void);
    
    void tetrads_Receiving(void);
    
    void ED_Calculation(void);
    
    void NB_Calculation(void);
    
};


#endif /* worker_hpp */
