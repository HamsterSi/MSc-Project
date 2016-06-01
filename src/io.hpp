/*
 * This is the IO of the program, reading from files and write out data
 * are all done by the functions in the class IO.
 */

#ifndef io_hpp
#define io_hpp

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "mpi.h"

#include "tetrad.hpp"

using namespace std;

/* Parameters read from crd file */
typedef struct _Crd {
    
    int num_BP;           // Number of DNA Base Pairs
    
    int *num_Atoms_In_BP; // Number of atoms in each base pair
    
    int total_Atoms;      // Total number of atoms in DNA base pairs
    
    float *ini_BP_Crds;   // Initial BP coordinates;
    
}Crd;


/* Parameters read from prm file
 * Other parameters are stored into class "Tetrad" */
typedef struct _Prm {
    
    int num_Tetrads; // Total numbers of overlapped tetrads
    
}Prm;


/*
 * The IO class can read data from files ("crd" and "prm") and wrtie 
 * results out to new files (To be discussed and implemented).
 */
class IO {
    
public:
    
    Crd crd;
    
    Prm prm;
    
    Tetrad *tetrad;
    
public:
    
    void read_Prm(string prm_File);
    
    void read_Crd(string crd_File, bool redundant);
    
    void read_Initial_Crds(void);
    
    void write_Results(string output_File, float* velocities, float* coordinates);
    
};

#endif /* io_hpp */
