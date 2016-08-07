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
 * File:  worker.hpp
 * Brief: Declaration of the Worker class
 */

#ifndef worker_hpp
#define worker_hpp

#include <iostream>
#include <ctime>
#include "mpi.h"

#include "array.hpp"
#include "mpilib.hpp"
#include "edmd.hpp"
#include "tetrad.hpp"
#include "io.hpp"

using namespace std;

/**
 * Brief: The Worker class is mainly used to calculate ED/NB forces
 */
class Worker {
    
public:
    
    int rank, size;  // The MPI rank and size

    int num_Tetrads; // The number of tetrads
    
    Tetrad *tetrad;  // The tetrads array
    
    EDMD edmd;       // The EDMD class for EDMD calculation
    
    Array array;     // For 2D array allocation & deallocation
    
    MPI_Lib mpi;     // For creating MPI_Datatype
    
    int max_Atoms;   // The maximum number of atoms in tetrads
    
    int num_Pairs;   // The number of non-bonded pairs
    
    double ** pair_Lists; // The 2D array of NB pair lists
    
    int    ** NB_Index;   // The workload distribution of NB force caulcaiton
    
    int    ** ED_Index;   // The workload distribution of ED force caulcaiton
    
    double ** NB_Forces;  // The 2D array to store the NB forces
    
    MPI_Comm comm;        // The MPI communicator
    
    MPI_Datatype * MPI_ED_Forces; // For receiving ED forces, the random terms & coordinates
    
    MPI_Datatype   MPI_Crds;      // For sending all coordinates of tetrads
    
public:
    
    /**
     * Function:  The constructor of Worker class.
     *
     * Parameter: None
     *
     * Return:    None
     */
    Worker(void);
    
    /**
     * Function:  The destructor to deallocate memory & free MPI_Datatype
     *
     * Parameter: None
     *
     * Return:    None
     */
    ~Worker(void);
    
    /**
     * Function:  Workers receive the EDMD simualtion parameters, number of atoms
     *            and number of evecs in every tetrad from master
     *
     * Parameter: None
     *
     * Return:    None
     */
    void recv_Parameters(void);
    
    /**
     * Function:  Workers receive all tetrads from master
     *
     * Parameter: None
     *
     * Return:    None
     */
    void recv_Tetrads(void);
    
    /**
     * Function:  Receive the coordinates of all tetrads from master, compute NB 
     *            forces of specified tetrads, sum them up and reduce them to master
     *
     * Parameter: None
     *
     * Return:    None
     */
    void EDNB_Calculation();
    
    /**
     * Function:  Receive new messages from master & calculate the ED/NB forces
     *
     * Parameter: None
     *
     * Return:    None
     */
    void force_Calculation(void);
    
    /**
     * Function:  Set the NB forces of tetrads to 0.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void empty_NB_Forces(void);

    
};

#endif /* worker_hpp */

