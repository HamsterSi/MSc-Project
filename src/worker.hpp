
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

class Worker {
    
public:
    
    int rank;           // The rank of worker process
    
    int num_Tetrads;    // The number of total tetrads
    
    int max_Atoms;      // The maximum number of atoms in tetrads
    
    Tetrad *tetrad;     // Tetrad array, used to stroe tetrads
    
    EDMD edmd;          // EDMD class, needs to call functions to calculate forces
    
    MPI_Status status;  // MPI Status
    
    MPI_Comm comm;      // MPI Communicator
    
public:
    
    Worker(void);
    
    ~Worker(void);
    
    void recv_Parameters(void);
    
    void recv_Tetrads(void);
    
    void ED_Calculation(void);
    
    void NB_Calculation(void);
    
    int  terminate(void);
    
};


#endif /* worker_hpp */
