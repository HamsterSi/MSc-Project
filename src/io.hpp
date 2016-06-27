/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                 MSc in High Performance Computing, EPCC                      *
 *                      The University of Edinburgh                             *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  io.hpp
 * Brief: Declaration of the IO class for input/output
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

#include "array.hpp"
#include "edmd.hpp"
#include "tetrad.hpp"

using namespace std;

/**
 * Brief: The data in crd file will be read into these variables & arrays
 */
typedef struct _Crd {
    
    int num_BP;       // Number of DNA Base Pairs
    
    int * BP_Atoms;   // Number of atoms in each base pair
    
    int total_Atoms;  // Total number of atoms in DNA
    
    double * BP_Vels; // Initial BP velocities;
    
    double * BP_Crds; // Initial BP coordinates;
    
}Crd;


/**
 * Brief: Parameters read from prm file, other data are read into the tetrad class array
 */
typedef struct _Prm {
    
    int num_Tetrads;  // Total numbers of overlapped tetrads
    
}Prm;


/**
 * Brief: The IO class reads data from files (mainly the "crd" and "prm" file) and wrtie
 *        results into new files. All input/output functions are in the IO class.
 */
class IO {
    
public:
    
    Crd crd;
    
    Prm prm;
    
    Tetrad *tetrad; // The array of tetrad class
    
    Array array;    // The array class for array operation
    
    int * displs;   // The displacements of base pairs
    
    int irest;      // Indicates to read which crd file
    int nsteps;     // The total steps of iterations
    int ntsync;     // The frequency of snyc
    int ntwt;       // The frequency of writing of trajectories
    int ntpr;       // The frequency of writing of info (Energies & temperature, crd file)
    int ncycs;      // Total cycles
    
    bool circular;  // The Shape of DNA
    
    // The strings used to store file paths
    string prm_File;
    string crd_File;
    string energy_File;
    string forces_File;
    string trj_File;
    string new_Crd_File;
    
public:
    
    /**
     * Function:  The constructor of IO class. Set default parameters.
     *
     * Parameter: None
     *
     * Return:    None
     */
    IO(void);
    
    /**
     * Function:  The destructor of IO class. Deallocate memorys.
     *
     * Parameter: None
     *
     * Return:    None
     */
    ~IO(void);
    
    /**
     * Function:   Read the config file & initialisation before the simualtion starts.
     *
     * Parameters: EDMD* edmd -> Some data will be read into the instance of EDMD class.
     *
     * Returns:    None.
     */
    void read_Cofig(EDMD* edmd);
    
    /**
     * Function:   Read the prm file. The DNA structure (linear/circular) needs to be
     *             taken into account in the prm-file generation.
     *
     * Parameters: None.
     *
     * Returns:    None.
     */
    void read_Prm(void);
    
    /**
     * Funtion:    Read a coordinates (.crd) file.
     *
     * Parameters: None
     *
     * Returns:    None.
     */
    void read_Crd(void);
    
    /**
     * Function:   Generate the displacements of base pairs.
     *             Set up some book-keeping stuff. "displs" store the displacements of base pairs.
     *             If 1st BP is 0, then 2nd BP displacement is 0+(3 * number of atoms in 1st BP)
     *
     * Parameters: None.
     *
     * Returns:    None.
     */
    void generate_Displacements(void);
    
    /**
     * Function:   Initialise velocities, coordinates and other parameters for tetrads from crd
     *
     * Parameters: None.
     *
     * Returns:    None.
     */
    void initialise_Tetrad_Crds(void);
    
    /**
     * Function:   Write out some array in a certain form
     *
     * Parameters: ofstream* fout -> The file stream pointer
     *             double* data   -> The data needs to write out
     *
     * Returns:    None.
     */
    void write_Template(ofstream* fout, double* data);
    
    void open_File(ofstream* fout, string file_Path);

    void write_Energy(ofstream* fout, int istep, double energies[]);
    
    void write_Trajectories(ofstream* fout, int istep, int index, double* coordinates);
    
    void close_File(ofstream* fout);
    
    /**
     * Function:   Write out energies & temperature of all tetrads.
     *
     * Parameters: int istep         -> The iterations of simulation
     *             double energies[] -> Energies & temperature of tetrads
     *
     * Returns:    None.
     */
    void write_Energies(int istep, double energies[]);
    
    /**
     * Function:   Write out all ED forces, random forces & NB forces
     *
     * Parameters: double* ED_Forces     -> The total ED forces of all DNA base pairs
     *             double* random_Forces -> The total random forces of all DNA base pairs
     *             double* NB_Forces     -> The total NB forces of all DNA base pairs
     *
     * Returns:    None.
     */
    void write_Forces(double* ED_Forces, double* random_Forces, double* NB_Forces);
    
    /**
     * Function:   Write out trajectories evry certain iterations
     *
     * Parameters: int istep           -> The iterations of simulation
     *             int total_Atoms     -> Total atoms of the DNA
     *             int index           -> The index of writing trajectory to which place
     *             double* coordinates -> The coordinates of all atoms in DNA
     *
     * Returns:    None.
     */
    void write_Trajectory(int istep, int index, double* coordinates);
    
    /**
     * Function:   Update the crd file for next simulation, write out velocities & coordinates
     *
     * Parameters: double* velocities  -> The velocities  of all atoms in DNA
     *             double* coordinates -> The coordinates of all atoms in DNA
     *
     * Returns:    None.
     */
    void update_Crd(double* velocities, double* coordinates);
    
};

#endif /* io_hpp */

