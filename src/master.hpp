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
 * File:  master.hpp
 * Brief: Declaration of the Master class
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
 * Brief: The Master class is responsible for initialsation and finalisation the simulation 
 *        as well as doing the simualtion on velocities, coordinates and other calculations.
 *        It requires message pssing between worker processes.
 */
class Master {
    
public:
    
    int size;    // The number of MPI processes
   
    EDMD edmd;   // The EDMD class for EDMD calculation
    
    IO io;       // The IO class for inputs and outputs
    
    Array array; // For 2D array allocation & deallocation
    
    MPI_Lib mpi; // For creating MPI_Datatype
    
    int       max_Atoms;  // The maximum number of atoms in tetrads
    
    int       num_Pairs;  // The number of non-bonded pairs
    
    double ** pair_Lists; // The 2D array of NB pair lists
    
    int    ** NB_Index;   // The workload distribution of NB force caulcaiton
    
    int    ** ED_Index;   // The workload distribution of ED force caulcaiton
    
    double ** NB_Forces;  // The 2D array to store the NB forces
    
    double * velocities;  // The velocities of the DNA
    
    double * coordinates; // The coordinates of the DNA
    
    MPI_Comm comm;        // The MPI communicator
    
    MPI_Datatype * MPI_ED_Forces; // For receiving ED forces, the random terms & coordinates
    
    MPI_Datatype   MPI_Crds;      // For sending all coordinates of tetrads

    
    
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
     * Function:  Master divides the workload (according to the pair lists) of the NB force
     *            calculation into similar size chunk.
     *            It is actually the displacements of the pair lists.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void generate_Indexes(void);
    
    /**
     * Function:  Master broadcast the NB parameters (pair lsit, pair index) to all workers
     *
     * Parameter: None
     *
     * Return:    None
     */
    void send_Workload_Indexes(void);

    /**
     * Function:  Master distrubutes the ED force calculation among worker processes, send
     *            tetrad index & coordinates to workers and receive forces & energy back.
     * Function:  Master distrubutes the NB force calculation among all MPI processes
     *            , send the coordinates of all tetrads to
     *            workers at first, then after the NB force calculation, a reduction
     *            operation is made to sum up all NB forces.
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
     * Function:  Master merge the velocities & coordinates of tetrad together
     *            & divide them by 4 (as they are fourfold overlapped)
     *
     * Parameter: None
     *
     * Return:    None
     */
    void merge_Vels_n_Crds(void);
    
    /**
     * Function:  Master writes the energies, temperature and trajectories.
     *
     * Parameter: int istep -> The iterations of simulation
     *
     * Return:    None
     */
    void write_Info(int istep);
    
    /**
     * Function:  Master writes a new coordinate file & keeps writing out
     *            on this file at certain frequency.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void write_Crds(void);
    
    /**
     * Function:  Master terminates workers & finalisation the simualtion
     *
     * Parameter: None
     *
     * Return:    None
     */
    void finalise(void);
    
    
};

#endif /* master_hpp */
