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
    
private:

    int num_Tetrads;    // The number of total tetrads
    
    Tetrad *tetrad;     // The Tetrad array to store all tetrads
    
    EDMD edmd;          // EDMD class to calculate ED/NB forces
    
    MPI_Comm comm;      // MPI Communicator
    
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
     * Function:  The destructor of Worker class. Deallocate memory, etc.
     *
     * Parameter: None
     *
     * Return:    None
     */
    ~Worker(void);
    
    /**
     * Function:  Workers receive the number of tetrads, edmd parameters, number of atoms
     *            and number of evecs in every tetrad from the master process
     *
     * Parameter: None
     *
     * Return:    None
     */
    void recv_Parameters(void);
    
    /**
     * Function:  Workers receive tetrads from master
     *
     * Parameter: None
     *
     * Return:    None
     */
    void recv_Tetrads(void);
    
    /**
     * Function:  Compute ED forces of tetrads.
     *
     * Parameter: int index[] -> The tetrad index & the index of force type
     *            MPI_Request send_Request[] -> The MPI send request
     *
     * Return:    None
     */
    void ED_Calculation(int index[], MPI_Request send_Request[]);
    
    /**
     * Function:  Compute NB forces of tetrads.
     *
     * Parameter: int index[] -> The tetrad index & the index of force type
     *            MPI_Request send_Request[] -> The MPI send request
     *
     * Return:    None
     */
    void NB_Calculation(int index[], MPI_Request send_Request[]);
    
    /**
     * Function:  Function for workers whose ranks are >= 2. Responsible for receive new tasks
     *            from master to calculate ED/NB forces and send forces & energies to rank 1.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void force_Calculation(void);
    
    /**
     * Function:  Function for the worker (rank 1). Receive the ED forces &
     *            coordinates from other workers.
     *
     * Parameter: int index  -> The index of the tetrad
     *            int source -> The source where the ED forces & coordinates from.
     *
     * Return:    None
     */
    void recv_ED_Forces(int index, int source);
    
    /**
     * Function:  Function for the worker (rank 1). Set all array elements of 
     *            the NB forces of tetrads to 0.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void empty_NB_Forces(void);
    
    /**
     * Function: Function for the worker (rank 1). Sum up the NB forces & NB 
     *           energy, electrostatic energy received from other workers.
     *
     * Parameter: int index[]        -> The indexes of tetrads
     *            int source         -> The source where the NB forces from.
     *            double* NB_Forces1 -> The NB forces of one tetrad
     *            double* NB_Forces2 -> The NB forces of one tetrad
     *
     * Return:    None
     */
    void sum_NB_Forces(int index[], int source, double* NB_Forces1, double* NB_Forces2);
    
    /**
     * Function:  Function for the worker (rank 1). Clip the NB forces into range
     *            (-1.0, 1.0) after received all NB forces & summed up.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void clip_NB_Forces(void);
    
    /**
     * Function:  Function for the worker (rank 1). It received singals from both the
     *            master & other workers to process the ED/NB forces.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void force_Processing(void);
    
};

#endif /* worker_hpp */

