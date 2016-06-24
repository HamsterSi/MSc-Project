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
 * File:  edmd.hpp
 * Brief: Declaration of a class for ED/MD simulation
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
 * Brief: Constants used in the DNA simulations.
 *        Merge with the global constants class, unify units
 */
typedef struct _Constants {
    
    double Boltzmann; // Boltzmann constant in internal units
    
    double timefac;   // Factor to convert time-realted parameter to internal units

}Constants;



/**
 * Brief: The EDMD class that implements the ED & NB forces calculation and other 
 *        calculations on velocities and coordinates.
 */
class EDMD {
    
public:
    
    Constants constants; // The DNA constants
    
    double dt;           // Timestep, in ps
    
    double gamma;        // The friction coefficient, in ps⁻¹
    
    double tautp;        // Berendsen temperature coupling parameter
    
    double temperature;  // Temperature, in Kelvin
    
    double scaled;       // Scale factor to scale ED forces
    
    double mole_Cutoff;  // Molecular cutoffs
    
    double atom_Cutoff;  // Atomic cutoffs
    
    double mole_Least;   // Molecules less than NB_Cutoff won't have NB ints.
    
public:
    
    /**
     * Function:  The constructor of EDMD class
     *
     * Parameter: None
     *
     * Return:    None
     */
    EDMD(void);
    
    /**
     * Function:  The destructor of EDMD class
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
     * Function:  Calculate ED forces for every tetrad, the results are stored in the Tetrad class
     *
     * Parameter: Tetrad* tetrad -> The instance of Tetrad cleass
     *
     * Return:    None
     */
    void calculate_ED_Forces(Tetrad* tetrad);
    
    /**
     * Function:  Calculate the LV random forces.
     *            Generate the Gaussian stochastic term. Assuming unitless.
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
     * Parameter: Tetrad* tetrad -> The instance of Tetrad cleass
     *
     * Return:    None
     */
    void calculate_Random_Forces(Tetrad* tetrad);
    
    /**
     * Function:  Calculate NB forces, results stored in Tetrad class
     *
     * Parameter: Tetrad* t1 -> The instance of Tetrad cleass
     *            Tetrad* t2 -> The instance of Tetrad cleass
     *
     * Return:    None
     */
    void calculate_NB_Forces(Tetrad* t1, Tetrad* t2);
    
    /**
     * Function:  Update velocities & Berendsen temperature control
     *
     * Parameter: Tetrad* tetrad -> The instance of Tetrad cleass
     *
     * Return:    None
     */
    void update_Velocities(Tetrad* tetrad); 
    
    /**
     * Function:  Update Coordinates of tetrads
     *
     * Parameter: Tetrad* tetrad -> The instance of Tetrad cleass
     *
     * Return:    None
     */
    void update_Coordinates(Tetrad* tetrad);
    
};

#endif /* parameters_hpp */

