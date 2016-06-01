
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
    
    string prm_File;
    
    string crd_File;
    
    string output_File;
    
    EDMD edmd;
    IO io;
    
    MPI_Comm comm;
    MPI_Status status;
    
    int size;
    int max_Atoms;
    int * displs;
    
    float * whole_Velocities;
    float * whole_Coordinates;
    
public:
    
    Master_Management(void);
    
    ~Master_Management(void);
    
    void initialise(void);
    
    void parameters_Sending(void);
    
    void tetrads_Sending(void);
    
    void force_Calculation(void);
    
    void velocity_Calculation(void);
    
    void coordinate_Calculation(void);
    
    void finalise(void);
    
};

#endif /* master_hpp */
