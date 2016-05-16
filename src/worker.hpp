
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
    
    int data[3]; // io.prm.num_Tetrads, io.prm.max_Atoms, io.prm.max_Evecs
    
    Tetrad *tetrad;
    
    EDMD edmd;
    
    MPI_Status status;
    
    MPI_Comm comm;
    
public:
    
    Worker_Management(void);
    
    void data_Receiving(void);
    
    void tetrad_Receiving(void);
    
    void ED_Calculation(void);
    
    void NB_Calculation(void);
    
};


#endif /* worker_hpp */
