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
 * File:  simulation.hpp
 * Brief: Contains two functions for the master and the workers respectively.
 */

#ifndef simulation_hpp
#define simulation_hpp

#include <iostream>
#include <ctime>
#include "mpi.h"

#include "mpilib.hpp"
#include "master.hpp"
#include "worker.hpp"

using namespace std;

/**
 * Function:  
 Define a Master class and manage the master working progress.
 *            It contorls the iteration structres of the simulation.
 *
 * Parameter: None
 *
 * Return:    None
 */
void master_Code(void);

/**
 * Function:  Define a Worker class and manage the worker working progress.
 *            If there is a message arrived, then it performs some functions according
 *            to the message tag.
 *
 * Parameter: None
 *
 * Return:    None
 */
void worker_Code(void);

#endif /* simulation_hpp */
