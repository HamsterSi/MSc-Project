
#ifndef parameters_hpp
#define parameters_hpp

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

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
    
    double Boltzmann; // Boltzmann constant in internal units
    
    double timefac;   // Factor to convert time-realted parameter to internal units

}Constants;


/*
 * The ED/MD simulation class
 */
class EDMD {
    
public:
    
    Constants  constants;
    
    double dt;          // Timestep, in ps
    double gamma;       // The friction coefficient, in ps⁻¹
    double tautp;       // Berendsen temperature coupling parameter
    
    double temperature; // Temperature, in Kelvin
    double scaled;      // Scale factor to scale ED forces
    
    double mole_Cutoff; // Molecular cutoffs
    double atom_Cutoff; // Atomic cutoffs
    double mole_Least;  // Molecules less than NB_Cutoff won't have NB ints.
    
public:
    
    EDMD(void);
    
    void initialise(double _dt, double _gamma, double _tautp, double _temperature, double _scaled, double _mole_Cutoff, double _atom_Cutoff, double _mole_Least);
    
    void calculate_ED_Forces(Tetrad* tetrad);
    
    void calculate_Random_Forces(Tetrad* tetrad);
    
    void generate_Pair_Lists(int pair_List[][2], int* effective_Pairs, int num_Tetrads, Tetrad* tetrad);
    
    void calculate_NB_Forces(Tetrad* t1, Tetrad* t2);
    
    void update_Velocities(Tetrad* tetrad, int index);
    
    void update_Coordinates(Tetrad* tetrad);
    
};



#endif /* parameters_hpp */
