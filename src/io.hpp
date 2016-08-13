/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                  MSc in High Performance Computing, EPCC                     *
 *                       The University of Edinburgh                            *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  io.hpp
 * Brief: The declaration of the IO class for input/output
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
 * Brief: The variables and arrays to be initialised from the coordinate file.
 *        The coordinate file contains the number of atoms, coordinates (and 
 *        velocities) of base pairs.
 */
typedef struct _Crd {
    
    int num_BP;       // The number of DNA Base Pairs
    
    int * BP_Atoms;   // The number of atoms in each base pair
    
    int total_Atoms;  // The total number of atoms in DNA
    
    double * BP_Vels; // Initial velocities of base pairs
    
    double * BP_Crds; // Initial coordinates of base pairs
    
}Crd;


/**
 * Brief: Parameters read from tetrad parameter file, 
 *        The other data from the Prm file are read into the tetrad class array
 */
typedef struct _Prm {
    
    int num_Tetrads;  // Total numbers of the overlapped tetrads
    
}Prm;


/**
 * Brief: The IO class with the input/output functions to read data from the 
 *        coordinate file and the tetrad parameter file, and to wrtie the trajectories,
 *        the enegies, etc. to new files.
 */
class IO {
    
public:
    
    Crd crd;        // For reading the coordinate file
    
    Prm prm;        // For reading the tetrad parameter file
    
    Tetrad *tetrad; // The tetrad array
    
    int * displs;   // The displacements of base pairs
    
    int irest;      // Indicates to read which crd file
    int nsteps;     // The total steps of iterations
    int ntsync;     // The frequency of synchronization
    int ntwt;       // The frequency of writing of energy & trajectory
    int ntpr;       // The frequency of updating the crd file

    // The strings of the input/output file paths
    string prm_File;
    string crd_File;
    string energy_File;
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
     * Function:  The destructor of IO class. Deallocate memory space of dynamical arrays.
     *
     * Parameter: None
     *
     * Return:    None
     */
    ~IO(void);
    
    /**
     * Function:   Read the config file & initialisation before the simualtion starts.
     *
     * Parameters: EDMD* edmd -> Some data is read into the instance of the EDMD class.
     *
     * Returns:    None.
     */
    void read_Cofig(EDMD* edmd);
    
    /**
     * Function:   Read the tetrad parameter file (mainly read into the tetrad array).
     *
     * Parameters: None.
     *
     * Returns:    None.
     */
    void read_Prm(void);
    
    /**
     * Funtion:    Read the coordinate file of the DNA base pairs.
     *
     * Parameters: None
     *
     * Returns:    None.
     */
    void read_Crd(void);
    
    /**
     * Function:   Generate the displacements of base pairs.
     *             Set up some book-keeping stuff. "displs" store the displacements of base pairs.
     *             If 1st BP is 0, then 2nd BP displacement is 0+(3 * number of atoms of 1st BP)
     *
     * Parameters: None.
     *
     * Returns:    None.
     */
    void generate_Displacements(void);
    
    /**
     * Function:   Initialise the coordinates and other variables of tetrads according to the 
     *             data read from the coordinate file.
     *             Firstly check whether the data from the crd file and prm file matches or not.
     *
     * Parameters: None.
     *
     * Returns:    None.
     */
    void initialise_Tetrad_Crds(void);
    
    /**
     * Function:   Write out the energies & temperature of the DNA.
     *             All the energies & temperature of tetrads should be summed up before calling
     *             this function.
     *
     * Parameters: int istep         -> The iterations of the ED/MD simulation
     *             double energies[] -> Energies & temperature of the DNA
     *
     * Returns:    None.
     */
    void write_Energies(int istep, double energies[]);
    
    /**
     * Function:   Write out the trajectories of the DNA.
     *
     * Parameters: int istep           -> The iterations of the ED/MD simulation
     *             int index           -> The index of writing which parts of the trajectory
     *             double* coordinates -> The coordinates of the DNA
     *
     * Returns:    None.
     */
    void write_Trajectory(int istep, int index, double* coordinates);
    
    /**
     * Function:   Update the coordinate file.
     *             Write out the number of atoms, the velocities & coordinates of base pairs
     *
     * Parameters: double* velocities  -> The velocities  of the DNA
     *             double* coordinates -> The coordinates of the DNA
     *
     * Returns:    None.
     */
    void update_Crd_File(double* velocities, double* coordinates);
    
};

#endif /* io_hpp */

