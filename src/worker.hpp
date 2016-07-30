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
    
    Tetrad *tetrad;     // Tetrad array, used to store tetrads
    
    EDMD edmd;          // EDMD class, functions will be called to calculate forces
    
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
     * Parameter: int index[]            -> The tetrad index & the index of force type
     *            MPI_Request send_Rqt[] -> The MPI send request
     *
     * Return:    None
     */
    void ED_Calculation(int index[], MPI_Request send_Rqt[]);
    
    /**
     * Function:  Compute NB forces of tetrads.
     *
     * Parameter: int index[]            -> The tetrad index & the index of force type
     *            MPI_Request send_Rqt[] -> The MPI send request
     *
     * Return:    None
     */
    void NB_Calculation(int index[], MPI_Request send_Rqt[]);
    
    /**
     * Function:  Responsible for receive instructions from master to calculate ED/NB forces and
     *            send forces & energies back to master.
     *
     * Parameter: None
     *
     * Return:    None
     */
    void force_Calculation(void);
    
    void recv_ED_Forces(int index, int source);
    
    void empty_NB_Forces(void);
    
    void sum_NB_Forces(int index[], int source, double* NB_Forces1, double* NB_Forces2);
    
    /**
     * Function:  Clip the NB forces into range (-1.0, 1.0)
     *
     * Parameter: None
     *
     * Return:    None
     */
    void clip_NB_Forces(void);
    
    void NB_Force_Processing(void);
    
};

#endif /* worker_hpp */

