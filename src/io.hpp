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
#include <sstream>
#include <cstdio>
#include "mpi.h"

#include "edmd.hpp"
#include "tetrad.hpp"

using namespace std;

/* Parameters read from crd file */
typedef struct _Crd {
    
    int num_BP;           // Number of DNA Base Pairs
    
    int * num_BP_Atoms;   // Number of atoms in each base pair
    
    int total_Atoms;      // Total number of atoms in DNA
    
    double * ini_BP_Vels; // Initial BP velocities;
    
    double * ini_BP_Crds; // Initial BP coordinates;
    
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
    
    int irest;     // Indicates to read which crd file
    int nsteps;    // The total steps of iterations
    int ntsync;    // The frequency of snyc
    int ntwt;      // The frequency of writing of traj
    int ntpr;      // The frequency of writing of info (Energies & temperature)
    int ncycs;     // Total cycles
    
    bool circular; // The Shape of DNA
    
    int * displs;  // The displacement of Base Pairs
    
    string prm_File;
    string crd_File;
    string energy_File;
    string forces_File;
    string trj_File;
    string new_Crd_File;
    
    Crd crd;
    
    Prm prm;
    
    Tetrad *tetrad;
    
public:
    
    IO(void);
    
    ~IO(void);
    
    void read_Cofig(EDMD* edmd);
    
    void read_Prm(void);
    
    void read_Crd(void);
    
    void generate_Displacements(void);
    
    void initialise_Tetrad_Crds(void);
    
    void write_Template(ofstream* fout, double* data);
    
    void write_Energies(int istep, double energies[]);
    
    void write_Forces(double* ED_Forces, double* random_Forces, double* NB_Forces);
    
    void write_Trajectory(int istep, int total_Atoms, int index, double* coordinates);
    
    void update_Crd(double* velocities, double* coordinates);
    
};

#endif /* io_hpp */
