
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
    
    int size;
    int max_Atoms;
    
    float *masses;
    float *displacement;
    
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
    
    void initialise(void);
    
    void parameters_Sending(void);
    
    void tetrads_Sending(void);
    
    void force_Passing(void);
    
    //void ED_Forces(void);
    
    //void NB_Forces(void);
    
    void LV_Forces(void);
    
    void total_Forces(void);
    
    void velocities(void);
    
    void coordinates(void);
    
    void finalise(void);
    
};

#endif /* master_hpp */
