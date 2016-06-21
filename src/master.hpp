
#ifndef master_hpp
#define master_hpp

#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <string>
#include "mpi.h"

#include "mpilib.hpp"
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
    
    int num_Pairs;
    int effective_Pairs;
    int * pair_List;

    double * velocities;
    double * coordinates;
    
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
    
    void send_Vels_n_Crds(void);
    
    void generate_Pair_Lists(void);
    
    void send_Tetrad_Index(int* i, int* j, int source);
    
    void cal_Forces(void);
    
    void cal_Velocities(void);
    
    void cal_Coordinate(void);
    
    void data_Processing(void);
    
    void write_Energy(int istep);
    
    void write_Forces(void);
    
    void write_Trajectories(void);
    
    void finalise(void);
    
};

#endif /* master_hpp */
