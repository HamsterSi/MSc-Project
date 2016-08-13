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
 * File:  mpilib.hpp
 * Brief: The declaration of the MPI_Lib class with user-defined functions
 */

#ifndef projectTest_hpp
#define projectTest_hpp

#include <iostream>
#include <cstddef>
#include "mpi.h"

#include "tetrad.hpp"


// MPI tags for message passing
#define TAG_DATA   1  // For passing the EDMD simulation parameters
#define TAG_TETRAD 2  // For passing the tetrads from master to workers
#define TAG_PAIRS  3  // For passing the pair list and displacements of workload
#define TAG_END    6  // For terminating simualtion
#define TAG_FORCE  7  // For ED and NB force calculation
#define TAG_NB     8  // For NB force calculation
#define TAG_ED     9  // For ED force calculation


using namespace std;


/**
 * Brief: The MPI_Lib class with functions to create and free various MPI data type
 *        for message passing between the master and the workers.
 */
class MPI_Lib{
    
public:
    
    /**
     * Function:  Create the MPI_Datatype for passing the tetrad parameters
     *
     * Parameter: MPI_Datatype* MPI_Tetrad -> The MPI data type for tetrads
     *            int num_Tetrads          -> The number of tetrads
     *            Tetrad* tetrad           -> The tetrad array
     *
     * Return:    None
     */
    static void create_MPI_Tetrad(MPI_Datatype* MPI_Tetrad, int num_Tetrads, Tetrad* tetrad);
    
    /**
     * Function:  Free the MPI_Datatype "MPI_Tetrad"
     *
     * Parameter: MPI_Datatype* MPI_Tetrad -> The MPI data type of tetrads
     *
     * Return:    None
     */
    static void free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad);
    
    /**
     * Function:  Create the MPI_Datatype for passing the ED forces
     *
     * Parameter: MPI_Datatype* MPI_Forces -> The MPI data type of the ED forces
     *            Tetrad* tetrad           -> The tetrad whose ED forces to be passed
     *
     * Return:    None
     */
    static void create_MPI_ED_Forces(MPI_Datatype* MPI_ED_Forces, Tetrad* tetrad);
    
    /**
     * Function:  Free the MPI_Datatype "MPI_ED_Forces"
     *
     * Parameter: MPI_Datatype* MPI_ED_Forces -> The MPI data type of ED forces
     *
     * Return:    None
     */
    static void free_MPI_ED_Forces(MPI_Datatype* MPI_ED_Forces);
    
    /**
     * Function:  Create the MPI_Datatype for passing the coordinates of all tetrads
     *
     * Parameter: MPI_Datatype* MPI_Crds -> The MPI data type of coordinates
     *            int num_Tetrads        -> The number of tetrads
     *            Tetrad* tetrad         -> he tetrad array
     *
     * Return:    None
     */
    static void create_MPI_Crds(MPI_Datatype* MPI_Crds, int num_Tetrads, Tetrad* tetrad);
    
    /**
     * Function:  Free the MPI_Datatype "MPI_Crds"
     *
     * Parameter: MPI_Datatype* MPI_Crds -> The MPI data type of coordinates
     *
     * Return:    None
     */
    static void free_MPI_Crds(MPI_Datatype* MPI_Crds);
    
};


#endif /* projectTest_hpp */

