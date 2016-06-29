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
    
    int max_Atoms;      // The maximum number of atoms in tetrads
    
    Tetrad *tetrad;     // Tetrad array, used to store tetrads
    
    double ** buffer;   // The buffer for send & receive data
    
    EDMD edmd;          // EDMD class, functions will be called to calculate forces
    
    Array array;        // The array class for array operation
    
    MPI_Status status;  // MPI Status
    
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
     * Function:  Compute ED forces of tetrads and send the ED forces & energy back to master.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void ED_Calculation(void);
    
    /**
     * Function:  Compute NB forces of tetrads. Send the NB forces & energies back to master.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void NB_Calculation(void);
    
    /**
     * Function:  Receive the terminate signal from master & terminate work
     *
     * Parameter: None
     *
     * Return:    None
     */
    int terminate(void);
    
};

#endif /* worker_hpp */

