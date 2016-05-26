
#ifndef parameters_hpp
#define parameters_hpp

#include <iostream>
#include <cmath>

#include "mpi.h"
#include "tetrad.hpp"
#include "./qcprot/qcprot.h"

using namespace std;

/*
 * Constants used in the DNA simulations.
 * TODO: merge with the global constants class, unify units
 *
 * Amber uses internally:
 *   - Lengths in Angstroms
 *   - Masses in atomic mass units,
 *   - Energies in kcal/mol
 *   - This means that the unit of time is 1/20.455 ps
 *   - VDW parameters: R* is in Angstroms, epsilon in kcal/mol.
 *   See http://ambermd.org/Questions/units.html
 */
typedef struct _Constants {
    
    float Boltzmann; // Boltzmann constant in internal units
    
    float timefac;   // Factor to convert time-realted parameter to internal units

}Constants;


/*
 * The ED/MD simulation class
 */
class EDMD {
    
public:
    
    Constants  constants;
    
    bool  circular;     // Whether the sistem has a circular or linear topology
    
    float RNG_Seed;    // Seed for the random number generator
    
    float dt;          // Timestep, in ps
    float gamma;       // The friction coefficient, in ps⁻¹
    float tautp;       // Berendsen temperature coupling parameter
    
    float temperature; // Temperature, in Kelvin
    float scaled;      // Scale factor to scale ED forces
    
public:
    
    EDMD(void);
    
    void calculate_ED_Forces(Tetrad* tetrad, float* forces, float scaled, int ED_Energy);
    
    void calculate_NB_Forces(Tetrad* tetrad1, Tetrad* tetrad2, float** NB_Forces, int NB_Energy, int Electrostatic_Energy);
    
    float generate_Stochastic_Term(float tetrad_ID);
    
    void generate_Pair_Lists(int pair_List[][2], int num_Tetrads, Tetrad* tetrad);
    
    void update_Velocities(float* velocities, float* ED_Forces, float* NB_Forces, float* masses, int total_Atoms);
    
    void update_Coordinates(float* coordinates, float* velocities, int total_Atoms);
    
};



#endif /* parameters_hpp */
