#include "edmd.hpp"


/*
 * Function:  The constructor of EDMD class
 *
 * Parameter: None
 *
 * Return:    None
 */
EDMD::EDMD(void) {
    
    constants.Boltzmann = 0.002;
    constants.timefac  = 20.455;
    
    dt  = 0.002;
    gamma = 2.0;
    tautp = 0.2;
    
    // Convert time-related parameters to internal units
    dt    *= constants.timefac;
    gamma /= constants.timefac;
    tautp *= constants.timefac;
    
    temperature = 300.0;
    scaled = constants.Boltzmann * temperature;
    
    mole_Cutoff = 30.0;
    atom_Cutoff = 10.0;
    mole_Least  =  5.0;
}




/*
 * Function:  Assign EDMD parameters
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::initialise(float _dt, float _gamma, float _tautp, float _temperature, float _scaled, float _mole_Cutoff, float _atom_Cutoff, float _mole_Least) {
    dt    =    _dt;
    gamma = _gamma;
    tautp = _tautp;
    temperature = _temperature;
    scaled = _scaled;
    mole_Cutoff = _mole_Cutoff;
    atom_Cutoff = _atom_Cutoff;
    mole_Least  = _mole_Least ;
}




/*
 * Function:  Calculate ED forces
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::calculate_ED_Forces(Tetrad* tetrad, float scaled) {
    
    int i, j;
    double rotmat[9], rmsd, temp = 0;
    
    // Allocate memory for temp arrays
    double * temp_Crds   = new double[3 * tetrad->num_Atoms];
    double * temp_Forces = new double[3 * tetrad->num_Atoms];
    double * proj        = new double[tetrad->num_Evecs];
    
    // Two arrays used in QCP rotation calculation
    double **avg_Crds  = new double *[3];
    double **crds      = new double *[3];
    for (i = 0; i < 3; i++) {
        avg_Crds[i]  = new double [tetrad->num_Atoms];
        crds[i]      = new double [tetrad->num_Atoms];
    }
    
    // Copy data from tetrads
    for (i = 0, j = 0; i < 3 * tetrad->num_Atoms && j < tetrad->num_Atoms; j++) {
        avg_Crds[0][j] = (double) tetrad->avg_Structure[ i ];
        avg_Crds[1][j] = (double) tetrad->avg_Structure[i+1];
        avg_Crds[2][j] = (double) tetrad->avg_Structure[i+2];
        
        crds[0][j] = (double) tetrad->coordinates[ i ];
        crds[1][j] = (double) tetrad->coordinates[i+1];
        crds[2][j] = (double) tetrad->coordinates[i+2];
        
        i += 3;
    }
    
    // Call QCP functions
    rmsd = CalcRMSDRotationalMatrix((double **) avg_Crds, (double **) crds, tetrad->num_Atoms, rotmat, NULL);
    
    // Transfer data to tetrads
    for (i = 0, j = 0; i < 3 * tetrad->num_Atoms && j < tetrad->num_Atoms; j++) {
        tetrad->avg_Structure[ i ] = (float) avg_Crds[0][j];
        tetrad->avg_Structure[i+1] = (float) avg_Crds[1][j];
        tetrad->avg_Structure[i+2] = (float) avg_Crds[2][j];
        
        tetrad->coordinates[ i ] = (float) crds[0][j];
        tetrad->coordinates[i+1] = (float) crds[1][j];
        tetrad->coordinates[i+2] = (float) crds[2][j];
        
        i += 3;
    }
    
    // Step 1: rotate x into the pcz frame of reference & remove average structure
    for (i = 0; i < 3 * tetrad->num_Atoms;) {
        temp_Crds[ i ] = rotmat[0]*tetrad->coordinates[i] + rotmat[1]*tetrad->coordinates[i+1] + rotmat[2]*tetrad->coordinates[i+2] - tetrad->avg_Structure[ i ];
        temp_Crds[i+1] = rotmat[3]*tetrad->coordinates[i] + rotmat[4]*tetrad->coordinates[i+1] + rotmat[5]*tetrad->coordinates[i+2] - tetrad->avg_Structure[i+1];
        temp_Crds[i+2] = rotmat[6]*tetrad->coordinates[i] + rotmat[7]*tetrad->coordinates[i+1] + rotmat[8]*tetrad->coordinates[i+2] - tetrad->avg_Structure[i+2];
        
        i += 3;
    }
    
    // Step 2: calculate projections
    for(i = 0; i < tetrad->num_Evecs; i++) {
        for(proj[i] = 0.0, j = 0; j < 3 * tetrad->num_Atoms; j++) {
            proj[i] += tetrad->eigenvectors[i][j] * temp_Crds[j];
        }
    }
    
    // Step 3 & Step 4
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        temp_Crds[i] = tetrad->avg_Structure[i];
        tetrad->ED_Forces[i] = 0.0;
    }
    for (i = 0; i < tetrad->num_Evecs; i++) {
        for (j = 0; j < 3 * tetrad->num_Atoms; j++) {
            // Step 3: re-embed the input coordinates in PC space - a sort of 'shake' procedure.
            //         Ideally this step is not needed, as stuff above should ensure all moves remain in PC subspace...
            temp_Crds[j] += tetrad->eigenvectors[i][j] * proj[i];
            
            // Step 4: calculate ED forces
            tetrad->ED_Forces[j] -= (tetrad->eigenvectors[i][j] * proj[i] * scaled / tetrad->eigenvalues[i]);
        }
    }
    
    // Step 5 & Step 6
    for (i = 0; i < 3 * tetrad->num_Atoms;) {
        // Step 5: rotate 'shaken' coordinates back into right frame
        tetrad->coordinates[ i ] = rotmat[0]*temp_Crds[i] + rotmat[1]*temp_Crds[i+1] + rotmat[2]*temp_Crds[i+2];
        tetrad->coordinates[i+1] = rotmat[3]*temp_Crds[i] + rotmat[4]*temp_Crds[i+1] + rotmat[5]*temp_Crds[i+2];
        tetrad->coordinates[i+2] = rotmat[6]*temp_Crds[i] + rotmat[7]*temp_Crds[i+1] + rotmat[8]*temp_Crds[i+2];
        
        // Step 6: rotate forces back to original orientation of coordinates
        temp_Forces[ i ] = rotmat[0]*tetrad->ED_Forces[i] + rotmat[1]*tetrad->ED_Forces[i+1] + rotmat[2]*tetrad->ED_Forces[i+2];
        temp_Forces[i+1] = rotmat[3]*tetrad->ED_Forces[i] + rotmat[4]*tetrad->ED_Forces[i+1] + rotmat[5]*tetrad->ED_Forces[i+2];
        temp_Forces[i+2] = rotmat[6]*tetrad->ED_Forces[i] + rotmat[7]*tetrad->ED_Forces[i+1] + rotmat[8]*tetrad->ED_Forces[i+2];
        
        i += 3;
    }

    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        tetrad->ED_Forces[i] = temp_Forces[i];
    }
    
    // Step 7: calculate the 'potential energy' (in units of kT)
    for (tetrad->energies[0] = 0.0, i = 0; i < tetrad->num_Evecs; i++) {
        tetrad->energies[0] += proj[i] * proj[i] / tetrad->eigenvalues[i];
    }
    tetrad->energies[0] *= 0.5 * scaled; // ED Energy
    
    // Deallocate memory
    for (i = 0; i < 3; i++) {
        delete [] avg_Crds[i]; delete [] crds[i];
    }
    delete [] avg_Crds; delete [] crds;
    delete [] temp_Forces; delete[] temp_Crds; delete [] proj;
    
}




/*
 * Function:  Calculate the LV random forces.
 *            Generate the Gaussian stochastic term. Assuming unitless.
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::calculate_Random_Forces(Tetrad* tetrad, int rank) {
    
    int i, j, RNG_Seed = 13579;
    float random, s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472;
    float half = 0.5, r1 = 0.27597, r2 = 0.27846, u, v, x, y, q;
    float* noise_Factor = new float[3 * tetrad->num_Atoms];

    srand((unsigned)RNG_Seed + rank);
    
    // Noise factors, sum(noise_Factor) = 3594.75 when gmma = 2.0
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        noise_Factor[i]  = sqrt(2.0 * gamma * scaled * tetrad->masses[i] / dt);
        tetrad->random_Forces[i] = 0.0;
    }
    
    // Calculate random forces;
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        /*
         ! Adapted from the following Fortran 77 code
         !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
         !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
         !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
         
         !  The function random_normal() returns a normally distributed pseudo-random
         !  number with zero mean and unit variance.
         
         !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
         !  and J.F. Monahan augmented with quadratic bounding curves.
         */
        {
            // Generate P = (u,v) uniform in rectangle enclosing acceptance region
            while (1) {
                
                u = (float)(rand()/(float)RAND_MAX);
                v = (float)(rand()/(float)RAND_MAX);
                v = 1.7156 * (v - half);
                
                x = u - s;
                y = abs(v) - t;
                q = x*x + y * (a*y - b*x);
                
                // Accept P if inside inner ellipse
                if (q < r1) { break; }
                // Reject P if outside outer ellipse
                if (q > r2) { continue; }
                // Reject P if outside acceptance region
                if (v*v < -4.0 * log(u) * (u*u)) { break; }
            }
            
            // Return ratio of P's coordinates as the normal deviate
            random = v/u;
        }
        
        tetrad->random_Forces[i] = random * noise_Factor[i];
    }
    
    delete []noise_Factor;
}




/*
 * Function:  Generate the pair list of tetrads.
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::generate_Pair_Lists(int pair_List[][2], int* effective_Pairs, int num_Tetrads, Tetrad* tetrad) {
    
    int i, j, k, num_Pairs;
    float r, ** com = new float * [num_Tetrads];
    for (i = 0; i < num_Tetrads; i++) { com[i] = new float[3]; }
    
    // The centre of mass (actually, centre of geom)
    // com(1) = sum(x(1:(3*natoms-2):3))/natoms
    // com(2) = sum(x(2:(3*natoms-1):3))/natoms
    // com(3) = sum(x(3:(3*natoms):3))  /natoms
    for (i = 0; i < num_Tetrads; i++) {
        
        com[i][0] = com[i][1] = com[i][2] = 0.0;
        for (j = 0; j < 3 * tetrad[i].num_Atoms; ) {
            com[i][0] += tetrad[i].coordinates[ j ];
            com[i][1] += tetrad[i].coordinates[j+1];
            com[i][2] += tetrad[i].coordinates[j+2];
            j += 3;
        }
        com[i][0] /= tetrad[i].num_Atoms;
        com[i][1] /= tetrad[i].num_Atoms;
        com[i][2] /= tetrad[i].num_Atoms;
    }
    
    // Loop to generate pairlists
    num_Pairs = 0; (* effective_Pairs) = 0;
    for (i = 0; i < num_Tetrads; i++) {
        for (j = i + 1; j < num_Tetrads; j++, num_Pairs++) {

            // If r exceeds mole_Cutoff then no interaction between these two mols
            // r = sum( (com(:,i)-com(:,j)) * (com(:,i)-com(:,j)) )
            for (r = 0.0, k = 0; k < 3; k++) {
                r += (com[i][k] - com[j][k]) * (com[i][k] - com[j][k]);
            }
            
            if ((r < (mole_Cutoff * mole_Cutoff)) && (abs(i - j) > mole_Least) &&
                (abs(i - j) < (num_Tetrads - mole_Least))) {
                
                pair_List[num_Pairs][0] = i;
                pair_List[num_Pairs][1] = j;
                
                (* effective_Pairs)++;
                
            } else {
                pair_List[num_Pairs][0] = -1;
                pair_List[num_Pairs][1] = -1;
            }
            
        }
    }
    
    for (i = 0; i < num_Tetrads; i++) { delete [] com[i]; }
    delete [] com;
    
}




/*
 * Function:  Calculate NB forces
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::calculate_NB_Forces(Tetrad* t1, Tetrad* t2) {
    
    int i, j;
    float dx, dy, dz, sqdist;
    float a, pair_Force;
    float krep = 100.0; // krep: soft repulsion constant
    float q;            // q: num_atoms vectors of charges
    float qfac = 0.0;   // qfac: electrostatics factor, set up for dd-dielectric constant of 4r, qfac=332.064/4.0, no electrostatics...
    
    // Initialise energies & NB forces
    t1->energies[1] = t1->energies[2] = 0.0;
    for (i = 0 ; i < t1->num_Atoms; i++) {
        t1->NB_Forces[i] = 0.0;
    }
    t2->energies[1] = t2->energies[2] = 0.0;
    for (i = 0 ; i < t2->num_Atoms; i++) {
        t2->NB_Forces[i] = 0.0;
    }
    
    for (i = 0; i < t1->num_Atoms; i++) {
        for (j = 0; j < t2->num_Atoms; j++) {
            
            dx = t1->coordinates[3 * i] - t2->coordinates[3 * j];
            dy = t1->coordinates[3*i+1] - t2->coordinates[3*j+1];
            dz = t1->coordinates[3*i+2] - t2->coordinates[3*j+2];
            
            // Avoid div0 (full atom overlap, almost impossible)
            sqdist = max(dx*dx + dy*dy + dz*dz, (float) 1e-9);
            
            // If this is the 1st cycle, cull the pairlist
            if (sqdist < (atom_Cutoff * atom_Cutoff)) {

                a = max(0.0, (2.0 - sqdist));
                q = t1->abq[3*i+2] * t2->abq[3*j+2];
                
                // NB Energy & Electrostatic Energy
                t1->energies[1] += 0.25 * krep * a * a;
                t1->energies[2] += 0.5 * qfac * q * sqdist;
                
                // NB forces
                pair_Force = -2.0 * krep * a - 2.0 * qfac * q / (sqdist * sqdist);
                t1->NB_Forces[3 * i] -= dx * pair_Force;
                t1->NB_Forces[3*i+1] -= dy * pair_Force;
                t1->NB_Forces[3*i+2] -= dz * pair_Force;
                
                t2->NB_Forces[3 * j] += dx * pair_Force;
                t2->NB_Forces[3*j+1] += dy * pair_Force;
                t2->NB_Forces[3*j+2] += dz * pair_Force;
            }
            
        }
    }
    
}




/*
 * Function:  Update velocities & Berendsen temperature control
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::update_Velocities(Tetrad* tetrad, int index) {
    
    int i;
    float kentical_Energy = 0.0;
    float target_KE;
    float tscal;  // Berendsen T-coupling factor
    float gamfac; // Velocity scale factor
    float max_Velocity = 10.0;
    
    gamfac = 1.0 / (1.0 + gamma * dt);
    
    // Simple Langevin dynamics
    float v1 = 0.0, v2 = 0.0;
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        tetrad->velocities[i] = (tetrad->velocities[i] + tetrad->ED_Forces[i] * dt + (tetrad->random_Forces[i] + tetrad->NB_Forces[i]) * dt / tetrad->masses[i]) * gamfac;
        //tetrad->velocities[i] = min( max_Velocity, tetrad->velocities[i]);
        //tetrad->velocities[i] = max(-max_Velocity, tetrad->velocities[i]);
        v1 += tetrad->velocities[i];
    }
    
    // Berendsen temperature control
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        kentical_Energy += 0.5 * tetrad->masses[i] * tetrad->velocities[i] * tetrad->velocities[i];
    }
    
    target_KE = 0.5 * scaled * 3 * tetrad->num_Atoms;
    tscal = sqrt(1.0 + (dt/tautp) * ((target_KE/kentical_Energy) - 1.0));
    
    // Calculate temperature of tetrad
    tetrad->temperature = kentical_Energy * 2 / (constants.Boltzmann * 3 * tetrad->num_Atoms);
    tetrad->temperature *= tscal * tscal;
    
    // Update velocities
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        tetrad->velocities[i] *= tscal;
        v2 += tetrad->velocities[i];
    }
    
    ///////////////////////////
    
    float ed = 0.0, ran = 0.0, nb = 0.0, v = 0.0, c = 0.0;
    for (int j = 0; j < 3 * tetrad->num_Atoms; j++) {
        ed  += tetrad->ED_Forces[j];
        ran += tetrad->random_Forces[j];
        nb  += tetrad->NB_Forces[j];
    }
    cout << "Index: "   << fixed << setw(2) << index
        << ", ED: "     << fixed << setprecision(4) << ed
        << ",\tRan: "   << fixed << setprecision(3) << ran
        << ",\tNB: "    << fixed << setprecision(4) << nb
        << ",\tKE: "    << fixed << setprecision(4) << kentical_Energy
        << ",\tTscal: " << fixed << setprecision(4) << tscal
        << ",\tT: "     << fixed << setprecision(4) << tetrad->temperature
        << ",\tVold: "  << fixed << setprecision(4) << v1
        << ",\tVnew: "  << fixed << setprecision(4) << v2 << endl;

}




/*
 * Function:  Update Coordinates of tetrads
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::update_Coordinates(Tetrad* tetrad) {
    
    for (int i = 0; i < 3 * tetrad->num_Atoms; i++) {
        tetrad->coordinates[i] += tetrad->velocities[i] * dt;
    }
}

