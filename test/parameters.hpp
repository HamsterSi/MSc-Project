//
//  parameters.hpp
//  
//
//  Created by Zhuowei Si on 06/04/2016.
//
//

#ifndef parameters_hpp
#define parameters_hpp

#include <iostream>
#include <string>
#include <ctime>

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

typedef struct _Energies {
    float eed; // ED energy
    float evdw;
    float eel;
    float ke;  // Current kinetic energy
    float ke0; // Target KE for Berendsen
}Energies;

/*
 * Parameters used in EDMD simulations.
 */
typedef struct _Parameters {
    
    bool  circular;    // Whether the sistem has a circular or linear topology
    float dt;          // Timestep, in ps
    float gamma;       // The friction coefficient, in ps⁻¹
    float tautp;       // Berendsen temperature coupling parameter
    float RNG_Seed;    // Seed for the random number generator
    float temperature; // Temperature, in Kelvin
    float scaled;      // Scale factor to scale ED forces
    float gamfac;      // Velocity scale factor
    float mole_Cutoff; // Molecular cutoffs
    float atom_Cutoff; // Atomic cutoffs
    float NB_Cutoff;   // Molecules less than NB_Cutoff won't have NB ints.
    
}Parameters;


/*
 * The ED/MD simulation class
 */
class EDMD {
    
public:
    
    Constants  DNA_Constants;
    
    Energies   EDMD_Energies;
    
    Parameters EDMD_Parameters;
    
public:
    
    EDMD(void);
    
    ~EDMD(void);
    
    void generate_Stochastic_Term(void);
    
    void parameters_Setting(void);

};



#endif /* parameters_hpp */
