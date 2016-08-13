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
 * File:  edmd.hpp
 * Brief: The declaration of the EDMD class for ED/MD simulation
 */

#ifndef parameters_hpp
#define parameters_hpp

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "mpi.h"

#include "./qcprot/qcprot.h"
#include "array.hpp"
#include "tetrad.hpp"

using namespace std;


/*
 * Brief: Constants used in the DNA ED/MD simulations.
 *        Merge with the global constants class, unify units
 */
typedef struct _Constants {
    
    double Boltzmann; // Boltzmann constant in internal units
    
    double timefac;   // Factor to convert time-realted parameter to internal units

}Constants;



/**
 * Brief: The EDMD class that implements the ED/NB force calculation, as well as
 *        the velocity and cooridnate calculation.
 */
class EDMD {
    
public:
    
    Constants constants; // The DNA constants
    
    double dt;           // Timestep, in ps
    
    double gamma;        // The friction coefficient, in ps⁻¹
    
    double tautp;        // Berendsen temperature coupling parameter
    
    double temperature;  // Temperature, in Kelvin
    
    double scaled;       // Scale factor to scale ED forces
    
    double mole_Cutoff;  // Molecular level cutoffs, used between tetrads
    
    double atom_Cutoff;  // Atomic level cutoffs, used between the atoms of two tetrads
    
    double mole_Least;   // Molecules farther than the mole_Least won't have NB forces
    
public:
    
    /**
     * Function:  The constructor of the EDMD class
     *
     * Parameter: None
     *
     * Return:    None
     */
    EDMD(void);
    
    /**
     * Function:  The destructor of the EDMD class
     *
     * Parameter: None
     *
     * Return:    None
     */
    ~EDMD(void);
    
    /**
     * Function:  Assign parameters for EDMD class
     *
     * Parameter: double _dt          -> Timestep, in ps
     *            double _gamma       -> The friction coefficient, in ps⁻¹
     *            double _tautp       -> Berendsen temperature coupling parameter
     *            double _temperature -> Temperature, in Kelvin
     *            double _scaled      -> Scale factor to scale ED forces
     *            double _mole_Cutoff -> Molecular cutoffs
     *            double _atom_Cutoff -> Atomic cutoffs
     *            double _mole_Least  -> Molecules less than NB_Cutoff won't have NB ints.
     *
     * Return:    None
     */
    void initialise(double _dt, double _gamma, double _tautp, double _temperature, double _scaled, double _mole_Cutoff, double _atom_Cutoff, double _mole_Least);
    
    /**
     * Function:  Calculate ED forces of tetrad
     *
     * Parameter: Tetrad* tetrad -> The tetrad whose ED forces to be calculated
     *
     * Return:    None, the ED forces are stored in the tetrad itself
     */
    void calculate_ED_Forces(Tetrad* tetrad);
    
    /**
     * Function:  Generate the Gaussian stochastic term. Assuming unitless.
     *
     *            Part of code is adapted from the following Fortran 77 code
     *            !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
     *            !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
     *            !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
     *
     *            !  The function returns a normally distributed pseudo-random
     *            !  number with zero mean and unit variance.
     *
     *            !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
     *            !  and J.F. Monahan augmented with quadratic bounding curves.
     *
     * Parameter: Tetrad* tetrad -> The tetrad whose random terms to be calculated
     *            int rank       -> The rank of worker used as the random seed
     *
     * Return:    None
     */
    void calculate_Random_Terms(Tetrad* tetrad, int rank);
    
    /**
     * Function:  Calculate the NB forces between two interacting tetrads
     *
     * Parameter: Tetrad* t1 -> The tetrad whose NB forces to be calculated
     *            Tetrad* t2 -> The tetrad whose NB forces to be calculated
     *
     * Return:    None, the NB forces are stored in two tetrads
     */
    void calculate_NB_Forces(Tetrad* t1, Tetrad* t2);
    
    /**
     * Function:  Update the velocities of tetrad (Berendsen temperature control applied)
     *
     * Parameter: Tetrad* tetrad -> The tetrad whose velocities to be calculated
     *
     * Return:    None
     */
    void update_Velocities(Tetrad* tetrad); 
    
    /**
     * Function:  Update Coordinates of tetrad
     *
     * Parameter: Tetrad* tetrad -> The tetrad whose cooridnates to be calculated
     *
     * Return:    None
     */
    void update_Coordinates(Tetrad* tetrad);
    
};

#endif /* parameters_hpp */

