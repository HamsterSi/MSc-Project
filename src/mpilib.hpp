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
 * File:  mpilib.hpp
 * Brief: Declaration of a class with some specific MPI funcitons
 */

#ifndef projectTest_hpp
#define projectTest_hpp

#include <iostream>
#include <cstddef>
#include "mpi.h"

#include "tetrad.hpp"

// Define the MPI tag for message passing
#define TAG_DATA   1    // For passing parameters (number of tetrads, etc.)
#define TAG_ED     2    // For message passing of ED forces calculation
#define TAG_NB     3    // For message passing of NB forces calculation
#define TAG_DEATH  4    // For terminating the simulation of workers
#define TAG_TETRAD 5    // For passing tetrads between master and workers

using namespace std;


/**
 * Brief: A class which contains two functions for creating and freeing the MPI data type
 *        for tetrads.
 */
class MPI_Library{
    
public:
    
    /**
     * Function:  Create the MPI_Datatype for tetrads
     *
     * Parameter: MPI_Datatype* MPI_Tetrad -> The MPI data type for tetrad
     *            Tetrad* tetrad           -> The instance of a tetrad
     *
     * Return:    None
     */
    static void create_MPI_Tetrad(MPI_Datatype* MPI_Tetrad, Tetrad* tetrad);
    
    /**
     * Function:  Free the MPI_Datatype of tetrads
     *
     * Parameter: MPI_Datatype* MPI_Tetrad -> The MPI data type of tetrads
     *
     * Return:    None
     */
    static void free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad);
    
};


#endif /* projectTest_hpp */

