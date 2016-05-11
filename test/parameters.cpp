//
//  parameters.cpp
//  
//
//  Created by Zhuowei Si on 07/04/2016.
//
//

#include "parameters.hpp"

EDMD::EDMD(void) {
    
    DNA_Constants.Boltzmann = 0.002;
    DNA_Constants.timefac = 20.455;
    
    EDMD_Parameters.circular = true;
    EDMD_Parameters.dt = 0.002;
    EDMD_Parameters.gamma = 0.4;
    EDMD_Parameters.tautp = 0.2;
    EDMD_Parameters.RNG_Seed = 0.0;
    EDMD_Parameters.temperature = 300.0;
    EDMD_Parameters.scaled = 0.0;
    EDMD_Parameters.gamfac = 0.0;
    EDMD_Parameters.mole_Cutoff = 40.0;
    EDMD_Parameters.atom_Cutoff = 10.0;
    EDMD_Parameters.NB_Cutoff = 4.0;
    
}

EDMD::~EDMD(void) {
    
}

void EDMD::parameters_Setting(void) {
    
    // Convert time-related parameters to internal units
    EDMD_Parameters.dt    *= DNA_Constants.timefac;
    EDMD_Parameters.gamma /= DNA_Constants.timefac;
    EDMD_Parameters.tautp *= DNA_Constants.timefac;
    
    // Scale factor
    EDMD_Parameters.scaled = DNA_Constants.Boltzmann * EDMD_Parameters.temperature;
    
    // Random number generator
    
    // Velocity scale factor
    EDMD_Parameters.gamfac = 1.0 / (1.0 + EDMD_Parameters.gamma * EDMD_Parameters.dt);
}
















