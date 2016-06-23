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
    
    EDMD edmd;            // The EDMD class for calculation
    
    IO io;                // The IO class for input and output
    
    int size;             // The size of MPI processes
    
    int max_Atoms;        // The maximum number of atoms of tetrads
    
    int num_Pairs;        // The total number of pairs for NB forces
    
    int effective_Pairs;  // The number of pairs that have NB forces to be calculated
    
    int * pair_Lists;     // The array to store pair lists

    double * velocities;  // Store the velocities of all atoms in DNA
    
    double * coordinates; // Store the coordinates of all atoms in DNA
    
    MPI_Comm comm;        // The MPI communicator
    
    MPI_Status status;    // MPI status
    
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
    
    /*
     * Function:  Master initialises the simulation.
     *            Read-in tetrads parameters, memory allocation and other initialisations.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void initialise(void);
    
    /*
     * Function:  Master sends the number of tetrads, edmd parameters, number of atoms
     *            and number of evecs in every tetrad to worker processes
     *
     * Parameter: None
     *
     * Return:    None
     */
    void send_Parameters(void);
    
    /*
     * Function:  Master sends tetrads to workers
     *
     * Parameter: None
     *
     * Return:    None
     */
    void send_Tetrads(void);
    
    /*
     * Function:  Master sends the velocities and coordinates of tetrads to workers
     *
     * Parameter: None
     *
     * Return:    None
     */
    void send_Vels_n_Crds(void);
    
    /*
     * Function:  Generate the pair lists of tetrads for non-bonded forces calculation
     *
     * Parameter: None
     *
     * Return:    None
     */
    void generate_Pair_Lists(void);
    
    /*
     * Function:  Master send the index of Tetrads to workers for ED & NB forces calculation.
     *            It also calculate the random forces.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void send_Tetrad_Index(int* i, int* j, int source);
    
    /*
     * Function:  Master distrubute ED/NB forces calculation among worker processes,
     *            send tetrad indexes to workers and receive forces & energies back.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void cal_Forces(void);
    
    /*
     * Function:  Master calculates velocities of tetrads
     *
     * Parameter: None
     *
     * Return:    None
     */
    void cal_Velocities(void);
    
    /*
     * Function:  Master calculates coordinates of tetrads
     *
     * Parameter: None
     *
     * Return:    None
     */
    void cal_Coordinate(void);
    
    /*
     * Function:  Master processes the velocities & coordinates of DNA
     *
     * Parameter: None
     *
     * Return:    None
     */
    void data_Processing(void);
    
    /*
     * Function:  Master writes out energies of all tetrads
     *
     * Parameter: int istep -> the iterations of simulation
     *
     * Return:    None
     */
    void write_Energy(int istep);
    
    /*
     * Function:  Master writes out forces of all atoms in DNA
     *
     * Parameter: None
     *
     * Return:    None
     */
    void write_Forces(void);
    
    /*
     * Function:  Master writes the trajectories.
     *
     * Parameter: int istep -> The iterations of simulation
     *
     * Return:    None
     */
    void write_Trajectories(int istep);
    
    /*
     * Function:  Master writes a new crd file.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void write_Crds(void);
    
    /*
     * Function:  Master terminates worker processes
     *
     * Parameter: None
     *
     * Return:    None
     */
    void finalise(void);
    
};

#endif /* master_hpp */
