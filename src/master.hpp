
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
class Master_Management {
    
public:
    
    int size;
    int max_Atoms;
    int * displs;
    
    float * velocities;
    float * coordinates;
    
    EDMD edmd;
    IO io;
    
    MPI_Comm comm;
    MPI_Status status;
    
public:
    
    Master_Management(void);
    
    ~Master_Management(void);
    
    
    // Functions of reading files & initialisation
    void read_Cofig(void);
    
    void read_Prm(void);
    
    void read_Crd(void);
    
    void initialise(void);
    
    
    // Functions of data passing & calculation
    void send_Parameters(void);
    
    void send_Tetrads(void);
    
    void force_Passing(void);
    
    void cal_Velocities(void);
    
    void cal_Coordinate(void);
    
    
    // Functions of writing out data & finalisation
    void write_Energy(void);
    
    void write_Forces(void);
    
    void write_Trajectory(void);
    
    void update_Crd_File(void);
    
    void finalise(void);
    
};

#endif /* master_hpp */
