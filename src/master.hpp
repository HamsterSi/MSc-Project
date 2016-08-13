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
 * File:  master.hpp
 * Brief: The declaration of the Master class
 */

#ifndef master_hpp
#define master_hpp

#include <iostream>
#include "mpi.h"

#include "array.hpp"
#include "mpilib.hpp"
#include "edmd.hpp"
#include "tetrad.hpp"
#include "io.hpp"

using namespace std;


/**
 * Brief: The Master class is responsible for initialsation and finalisation of the EDMD
 *        simulation, while it also calculates the velocities and coordinates of tetrads
 *        The IO is also called in the Master class.
 *        Message passing required with the workers.
 */
class Master {
    
public:
    
    int size;             // The number of MPI processes
   
    EDMD edmd;            // For updating velocities and coordiantes of tetrads
    
    IO io;                // The IO class for inputs and outputs
    
    Array array;          // For 2D array operations
    
    MPI_Lib mpi;          // For creating MPI_Datatype
    
    int      max_Atoms;   // The maximum number of atoms in tetrads
    
    int      num_Pairs;   // The number of non-bonded pairs
    
    double ** pair_Lists; // The 2D array of NB pair lists
    
    int    ** NB_Index;   // The workload distribution of NB force caulcaiton
    
    int    ** ED_Index;   // The workload distribution of ED force caulcaiton
    
    double ** NB_Forces;  // The 2D array to store the NB forces
    
    double * velocities;  // The velocities of the DNA
    
    double * coordinates; // The coordinates of the DNA
    
    MPI_Comm comm;        // The MPI communicator
    
    MPI_Datatype * MPI_ED_Forces; // For receiving ED forces & random terms
    
    MPI_Datatype   MPI_Crds;      // For sending the coordinates of all tetrads

    
    
public:
    
    /**
     * Function:  The constructor of Master class.
     *
     * Parameter: None
     *
     * Return:    None
     */
    Master(void);
    
    /**
     * Function:  The destructor of Master class. Deallocate memory of arrays.
     *
     * Parameter: None
     *
     * Return:    None
     */
    ~Master(void);
    
    /**
     * Function:  Master initialises the simulation.
     *            Read-in tetrads parameters, memory allocation and other initialisations.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void initialise(void);
    
    /**
     * Function:  Master sends the EDMD simualtion parameters, the number of atoms
     *            and number of evecs in every tetrad to worker processes
     *
     * Parameter: None
     *
     * Return:    None
     */
    void send_Parameters(void);
    
    /**
     * Function:  Master sends all tetrads to workers
     *
     * Parameter: None
     *
     * Return:    None
     */
    void send_Tetrads(void);
    
    /**
     * Function:  Calculate the center of mass (actually, centre of geom)
     *
     * Parameter: double** com -> Center of mass
     *
     * Return:    None
     */
    void cal_Centre_of_Mass(double** com);
    
    /**
     * Function:  Generate the pair lists of tetrads for non-bonded forces calculation
     *
     * Parameter: None
     *
     * Return:    None
     */
    void generate_Pair_Lists(void);
    
    /**
     * Function:  Master divides the workload (according to the pair lists) of the ED and 
     *            NB force calculation into similar size chunks.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void generate_Indexes(void);
    
    /**
     * Function:  Master sends the pair lsit, workload index to all workers
     *
     * Parameter: None
     *
     * Return:    None
     */
    void send_Workload_Indexes(void);

    /**
     * Function:  Master send the calculationg signal & the coordinates of all tetrads
     *            to all workers, and the workers then can start the ED/NB force calculation 
     *            according to the workload index sent before. 
     *            The master then receive the ED forces & sum up the NB forces with
     *            the MPI_Reduce operation.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void calculate_Forces(void);
    
    /**
     * Function:  Clip the NB forces into range (-1.0, 1.0) & assign NB forces to tetrads
     *
     * Parameter: None
     *
     * Return:    None
     */
    void process_NB_Forces(void);
    
    /**
     * Function:  Master calculates the velocities of all tetrads
     *
     * Parameter: None
     *
     * Return:    None
     */
    void update_Velocity(void);
    
    /**
     * Function:  Master calculates the coordinates of all tetrads
     *
     * Parameter: None
     *
     * Return:    None
     */
    void update_Coordinate(void);
    
    /**
     * Function:  Master merges the velocities & coordinates of tetrad together
     *            & divide them by 4 (as they are fourfold overlapped)
     *
     * Parameter: None
     *
     * Return:    None
     */
    void merge_Vels_n_Crds(void);
    
    /**
     * Function:  Master writes out the energies, temperature and trajectories.
     *
     * Parameter: int istep -> The iterations of simulation
     *
     * Return:    None
     */
    void write_Info(int istep);
    
    /**
     * Function:  Master updates the new coordinate file at a certain frequency.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void write_Crds(void);
    
    /**
     * Function:  Master terminates all workers & finalisation the simualtion
     *
     * Parameter: None
     *
     * Return:    None
     */
    void finalise(void);
    
    
};

#endif /* master_hpp */
