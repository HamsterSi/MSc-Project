
#ifndef master_hpp
#define master_hpp

#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "mpi.h"

#include "pool/pool.h"
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
    
    float *masses;
    
    int start_Index, end_Index;
    float *displacement;
    
    int max_Atoms;
    
    float *ED_Forces;
    float *total_ED_Forces;
    
    float **NB_Forces;
    float *total_NB_Forces;
    
    float *noise_Factor;
    float *langevin_Forces;
    
    float *total_Velocities;
    float *total_Coordinates;
    
    
public:
    
    Master_Management(void);
    
    ~Master_Management(void);
    
    void master_ED_Forces(void);
    
    void master_NB_Forces(void);
    
    void master_LV_Forces(void);
    
    void master_Total_Forces(void);
    
    void master_Velocities(void);
    
    void master_Coordinates(void);
    
};

#endif /* master_hpp */
