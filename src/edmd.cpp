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
    constants.timefac = 20.455;
    
    circular = false;
    
    dt    = 0.002;
    gamma = 0.4;
    tautp = 0.2;
    
    // Convert time-related parameters to internal units
    dt    *= constants.timefac;
    gamma /= constants.timefac;
    tautp *= constants.timefac;
    
    temperature = 300.0;
    
    // Scale factor
    scaled = constants.Boltzmann * temperature;
}


/*
 * Function:  Calculate Ed forces
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::calculate_ED_Forces(Tetrad* tetrad, float* ED_Forces, float scaled, int ED_Energy) {
    
    int i, j;
    double rotmat[9], rmsd, temp = 0;
    
    // Allocate memory for temp arrays
    double *temp_Crds   = new double[3 * tetrad->num_Atoms];
    double *temp_Forces = new double[3 * tetrad->num_Atoms];
    double *proj        = new double[tetrad->num_Evecs];
    
    // Two arrays used in QCP rotation calculation
    double **avg_Crds  = new double*[3];
    double **crds      = new double*[3];
    for (i = 0; i < 3; i++) {
        avg_Crds[i]  = new double[tetrad->num_Atoms];
        crds[i]      = new double[tetrad->num_Atoms];
    }
    
    // Copy data from tetrads
    for (i = 0, j = 0; i < 3 * tetrad->num_Atoms && j < tetrad->num_Atoms; j++) {
        avg_Crds[0][j] = (double) tetrad->avg_Structure[i];
        avg_Crds[1][j] = (double) tetrad->avg_Structure[i+1];
        avg_Crds[2][j] = (double) tetrad->avg_Structure[i+2];
        
        crds[0][j]     = (double) tetrad->coordinates[i];
        crds[1][j]     = (double) tetrad->coordinates[i+1];
        crds[2][j]     = (double) tetrad->coordinates[i+2];
        
        i += 3;
    }
    
    // Call QCP functions
    rmsd = CalcRMSDRotationalMatrix((double **) avg_Crds, (double **) crds, tetrad->num_Atoms, rotmat, NULL);
    
    // Transfer data to tetrads
    for (i = 0, j = 0; i < 3 * tetrad->num_Atoms && j < tetrad->num_Atoms; j++) {
        tetrad->avg_Structure[i]   = (float) avg_Crds[0][j];
        tetrad->avg_Structure[i+1] = (float) avg_Crds[1][j];
        tetrad->avg_Structure[i+2] = (float) avg_Crds[2][j];
        
        tetrad->coordinates[i]   = (float) crds[0][j];
        tetrad->coordinates[i+1] = (float) crds[1][j];
        tetrad->coordinates[i+2] = (float) crds[2][j];
        
        i += 3;
    }
    
    // Step 1: rotate x into the pcz frame of reference & remove average structure
    for (i = 0; i < 3 * tetrad->num_Atoms;) {
        temp_Crds[i]   = rotmat[0]*tetrad->coordinates[i] + rotmat[1]*tetrad->coordinates[i+1] + rotmat[2]*tetrad->coordinates[i+2] - tetrad->avg_Structure[i]   - tetrad->avg_Structure[i];
        temp_Crds[i+1] = rotmat[3]*tetrad->coordinates[i] + rotmat[4]*tetrad->coordinates[i+1] + rotmat[5]*tetrad->coordinates[i+2] - tetrad->avg_Structure[i+1] - tetrad->avg_Structure[i+1];
        temp_Crds[i+2] = rotmat[6]*tetrad->coordinates[i] + rotmat[7]*tetrad->coordinates[i+1] + rotmat[8]*tetrad->coordinates[i+2] - tetrad->avg_Structure[i+2] - tetrad->avg_Structure[i+1];
        
        i += 3;
    }
    
    // Step 2: calculate projections
    for(i = 0; i < tetrad->num_Evecs; i++) {
        proj[i] = 0.0;
        for(j = 0; j < 3 * tetrad->num_Atoms; j++) {
            proj[i] += tetrad->eigenvectors[i][j] * temp_Crds[j];
        }
    }
    
    // Step 3 & Step 4
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        temp_Crds[i] = tetrad->avg_Structure[i];
        ED_Forces[i] = 0.0;
    }
    for (i = 0; i < tetrad->num_Evecs; i++) {
        for (j = 0; j < 3 * tetrad->num_Atoms; j++) {
            // Step 3: re-embed the input coordinates in PC space - a sort of 'shake' procedure.
            //         Ideally this step is not needed, as stuff above should ensure all moves remain in PC subspace...
            temp_Crds[j] += tetrad->eigenvectors[i][j]*proj[i];
            
            // Step 4: calculate ED forces
            ED_Forces[j] -= (tetrad->eigenvectors[i][j]*proj[i]*scaled/tetrad->eigenvalues[i]);
        }
    }
    
    // Step 5 & Step 6
    for (i = 0; i < 3 * tetrad->num_Atoms;) {
        // Step 5: rotate 'shaken' coordinates back into right frame
        tetrad->coordinates[i]   = rotmat[0]*temp_Crds[i] + rotmat[1]*temp_Crds[i+1] + rotmat[2]*temp_Crds[i+2];
        tetrad->coordinates[i+1] = rotmat[3]*temp_Crds[i] + rotmat[4]*temp_Crds[i+1] + rotmat[5]*temp_Crds[i+2];
        tetrad->coordinates[i+2] = rotmat[6]*temp_Crds[i] + rotmat[7]*temp_Crds[i+1] + rotmat[8]*temp_Crds[i+2];
        
        // Step 6: rotate forces back to original orientation of coordinates
        temp_Forces[i]   = rotmat[0]*ED_Forces[i] + rotmat[1]*ED_Forces[i+1] + rotmat[2]*ED_Forces[i+2];
        temp_Forces[i+1] = rotmat[3]*ED_Forces[i] + rotmat[4]*ED_Forces[i+1] + rotmat[5]*ED_Forces[i+2];
        temp_Forces[i+2] = rotmat[6]*ED_Forces[i] + rotmat[7]*ED_Forces[i+1] + rotmat[8]*ED_Forces[i+2];
        
        i += 3;
    }

    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        ED_Forces[i] = temp_Forces[i];
    }
    
    // Step 7: calculate the 'potential energy' (in units of kT)
    for (i = 0; i < tetrad->num_Evecs; i++) {
        temp += proj[i] * proj[i] / tetrad->eigenvalues[i];
    }
    ED_Forces[ED_Energy] = scaled * 0.5 * temp; // ED Energy
    
    // Deallocate memory
    for (i = 0; i < 3; i++) {
        delete [] avg_Crds[i];
        delete [] crds[i];
    }
    delete [] avg_Crds;
    delete [] crds;
    delete [] temp_Crds;
    delete [] temp_Forces;
    delete [] proj;
    
}


/*
 * Function:  Calculate the LV random forces.
 *            Generate the Gaussian stochastic term. Assuming unitless.
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::calculate_Random_Forces(Tetrad* tetrad, float* random_Forces) {
    
    int i, j, rank, RNG_Seed = 13579;
    float random, s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472;
    float half = 0.5, r1 = 0.27597, r2 = 0.27846, u, v, x, y, q;
    float* noise_Factor = new float[3 * tetrad->num_Atoms];
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand((unsigned)RNG_Seed + rank);
    
    // Noise factors
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        noise_Factor[i]  = sqrt(2.0 * gamma * scaled * tetrad->masses[i] / dt);
        random_Forces[i] = 0.0;
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
        random_Forces[i] = random * noise_Factor[i];
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
void EDMD::generate_Pair_Lists(int pair_List[][2], int num_Tetrads, Tetrad* tetrad) {
    
    float mole_Cutoff = 40.0; // Molecular cutoffs
    float atom_Cutoff = 10.0; // Atomic cutoffs
    float mole_Least  = 4.0;  // Molecules less than NB_Cutoff won't have NB ints.
    
    int i, j, k, num_Pairs = 0;
    float r, com[num_Tetrads][3];
    
    // The centre of mass (actually, centre of geom)
    // com(1) = sum(x(1:(3*natoms-2):3))/natoms
    // com(2) = sum(x(2:(3*natoms-1):3))/natoms
    // com(3) = sum(x(3:(3*natoms):3))  /natoms
    for (i = 0; i < num_Tetrads; i++) {
        for (j = 0; j < 3 * tetrad[i].num_Atoms; ) {
            com[i][0] += tetrad[i].coordinates[j];
            com[i][1] += tetrad[i].coordinates[j+1];
            com[i][2] += tetrad[i].coordinates[j+2];
            j = j + 3;
        }
        com[i][0] /= tetrad[i].num_Atoms;
        com[i][1] /= tetrad[i].num_Atoms;
        com[i][2] /= tetrad[i].num_Atoms;
    }
    
    // Loop to generate pairlists
    for (i = 0; i < num_Tetrads; i++) {
        for (j = i+1; j < num_Tetrads; j++) {
            r = 0.0;
            
            // If r exceeds mole_Cutoff then no interaction between these two mols
            // r = sum( (com(:,i)-com(:,j)) * (com(:,i)-com(:,j)) )
            for (k = 0; k < 3; k++) {
                r += (com[k][i] - com[k][j]) * (com[k][i] - com[k][j]);
            }
            
            if ((r < mole_Cutoff * mole_Cutoff) && (abs(i-j) > mole_Least) &&
                (abs(i-j) < num_Tetrads - mole_Least)) {
                
                pair_List[num_Pairs][0] = i;
                pair_List[num_Pairs][1] = j;
                
            } else {
                pair_List[num_Pairs][0] = -1;
                pair_List[num_Pairs][1] = -1;
            }
            
            num_Pairs++;
        }
    }
}


/*
 * Function:  Calculate NB forces
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::calculate_NB_Forces(Tetrad* tetrad1, Tetrad* tetrad2, float* NB_Forces1, float* NB_Forces2, int energy_Index) {
    
    int i, j;
    int NB_Energy = energy_Index, Electrostatic_Energy = energy_Index;
    float dx, dy, dz, sqdist;
    float a, pair_Force;
    float krep = 100.0;   // krep: soft repulsion constant
    float q;              // q: num_atoms vectors of charges
    float qfac = 332.064; //qfac: electrostatics factor
    
    for (i = 0 ; i < energy_Index + 2; i++) {
        NB_Forces1[i] = NB_Forces2[i] = 0.0;
    }
    
    for (i = 0; i < tetrad1->num_Atoms; i++) {
        for (j = 0;  j < tetrad2->num_Atoms; j++) {
            
            dx = tetrad1->coordinates[3*i]   - tetrad2->coordinates[3*i];
            dy = tetrad1->coordinates[3*i+1] - tetrad2->coordinates[3*i+1];
            dz = tetrad1->coordinates[3*i+2] - tetrad2->coordinates[3*i+2];
            sqdist = dx*dx + dy*dy + dz*dz;
            
            a = max(0.0, 2.0-sqdist);
            q = tetrad1->abq[3*i+2] * tetrad2->abq[3*j+2];
            
            // NB energies
            NB_Forces1[NB_Energy] += 0.25 * krep * a * a;
            NB_Forces2[Electrostatic_Energy] += 0.5 * qfac * q * sqdist;
            
            // NB forces
            pair_Force = -2.0 * krep * a - 2.0 * qfac * q / (sqdist * sqdist);
            NB_Forces1[3*i]   -= dx * pair_Force;
            NB_Forces1[3*i+1] -= dy * pair_Force;
            NB_Forces1[3*i+2] -= dz * pair_Force;
            
            NB_Forces2[3*j]   += dx * pair_Force;
            NB_Forces2[3*j+1] += dy * pair_Force;
            NB_Forces2[3*j+2] += dz * pair_Force;
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
void EDMD::update_Velocities(Tetrad* tetrad) {
    
    int i;
    float kentical_Energy = 0.0;
    float target_KE;
    float actual_Temperature;
    float tscal;  // Berendsen T-coupling factor
    float gamfac; // Velocity scale factor
    
    gamfac = 1.0 / (1.0 + gamma * dt);
    
    // Simple Langevin dynamics
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        tetrad->velocities[i] = (tetrad->velocities[i] + tetrad->ED_Forces[i] * dt + tetrad->NB_Forces[i] * dt / tetrad->masses[i]) * gamfac;
    }
    
    // Berendsen temperature control
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        kentical_Energy += 0.5 * tetrad->masses[i] * tetrad->velocities[i] * tetrad->velocities[i];
    }
    
    actual_Temperature = kentical_Energy * 2 / (constants.Boltzmann * 3 * tetrad->num_Atoms);
    
    target_KE = 0.5 * scaled * 3 * tetrad->num_Atoms;
    
    tscal = sqrt(1.0 + (dt/tautp) * ((target_KE/kentical_Energy) - 1.0));
    
    // Update velocities
    for (i = 0; i < 3 * tetrad->num_Atoms; i++) {
        tetrad->velocities[i] = tetrad->velocities[i] * tscal;
    }
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

