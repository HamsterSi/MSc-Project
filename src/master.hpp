
#ifndef master_hpp
#define master_hpp

#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <string>
#include "mpi.h"

#include "mpilibrary.hpp"
#include "edmd.hpp"
#include "tetrad.hpp"
#include "io.hpp"

using namespace std;

/*
 * This is the class that will be executed by the master process.
 */
class Master {
    
public:
    
    int size;
    int max_Atoms;

    float * velocities;
    float * coordinates;
    
    EDMD edmd;
    IO io;
    
    MPI_Comm comm;
    MPI_Status status;
    
public:
    
    Master(void);
    
    ~Master(void);
    
    void initialise(void);
    
    void send_Parameters(void);
    
    void send_Tetrads(void);
    
    void send_Worker_Pairlists(int* j, int num_Pairs, int source, int pair_List[][2]);
    
    void force_Calculation(void);
    
    void cal_Velocities(void);
    
    void cal_Coordinate(void);
    
    void data_Processing(void);
    
    void write_Energy(void);
    
    void write_Forces(void);
    
    void write_Files(void);
    
    void finalise(void);
    
};

#endif /* master_hpp */
