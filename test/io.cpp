//
//  io.cpp
//  
//
//  Created by Zhuowei Si on 05/04/2016.
//
//

#include "io.hpp"


/*
 * Read-in the prm file. 
 *
 * Parameters:
 *      prm_File  -> path to the file to read the prm from.
 *
 * Returns: the tetrads plus information per-atom like masses and NBI parameters.
 * N.B. Things like circular structures are already taken into account in the prm-file generation.
 */
void IO::read_Prm(string prm_File) {
    
    ifstream fin;
    fin.open(prm_File, ios_base::in);
    
    if (fin.is_open()) {
        
        // Line 1:  Number of (overlapping) tetrads in the system (= number of base pairs in a circle)
        fin >> prm.num_Tetrads;
        
        int i, j, k;
        tetrad = new Tetrad[prm.num_Tetrads];
        
        // Rest of file has data for each tetrad as follows:
        for (i = 0; i < prm.num_Tetrads; i++) {
            
            // Line 2: Number of atoms in the tetrad, and number of eigenvectors
            fin >> tetrad[i].num_Atoms_In_Tetrad;
            fin >> tetrad[i].num_Evecs;
            
            // Allocate memory for tetrads
            tetrad[i].avg_Structure = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
            tetrad[i].masses        = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
            tetrad[i].abq           = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
            
            tetrad[i].eigenvalues   = new float[tetrad[i].num_Evecs];
            tetrad[i].eigenvectors  = new float* [tetrad[i].num_Evecs];
            for (j = 0; j < tetrad[i].num_Evecs; j++) {
                tetrad[i].eigenvectors[j] = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
            }
            
            tetrad[i].coordinates   = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
            tetrad[i].velocities    = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
            
            // Line 3 onwards: Reference (average) structure for the tetrad (x1,y1,z1,x2,y2,z2, etc as in .crd file)
            for (j = 0; j < 3 * tetrad[i].num_Atoms_In_Tetrad; j++) {
                fin >> tetrad[i].avg_Structure[j];
            }
            
            // Line ? onwards: Masses for each atom (amu)
            for (j = 0; j < 3 * tetrad[i].num_Atoms_In_Tetrad; ) {
                fin >> tetrad[i].masses[j];
                // Spread masses
                tetrad[i].masses[j+2] = tetrad[i].masses[j+1] = tetrad[i].masses[j];
                j = j + 3;
            }
            
            // Line ? onwards: Non-bonded parameters. Each line contains vdW parameters A and B, and partial charge q, for two atoms (e.g. 1st line: atoms 1 and 2, next line: atoms 3 and 4, etc.).
            for (j = 0; j < 3 * tetrad[i].num_Atoms_In_Tetrad; j++) {
                fin >> tetrad[i].abq[j];
            }
            
            // Next the eigenvector and eigenvalue data for this tetrad:
            // Line ?: The eigenvalue for the 1st (largest) eigenvector
            // Line ? onwards: Coefficients of the 1st eigenvector
            // Line ?: The eigenvalue for the second eigenvector
            // Line ? onwards: Coefficients of the 2nd eigenvector
            // (etc to the last eigenvector of this tetrad)
            for (j = 0; j < tetrad[i].num_Evecs; j++) {
                fin >> tetrad[i].eigenvalues[j];
                for (k = 0; k < 3 * tetrad[i].num_Atoms_In_Tetrad; k++) {
                    fin >> tetrad[i].eigenvectors[j][k];
                }
            }
        }
        
        fin.close();
        
    } else {
        cout << ">>> ERROR: Can not open prm file!" << endl;
        exit(1);
    }
}



/*
 * Read-in a coordinates (.crd) file.
 *
 * Parameters:
 *      crd_File  -> path to the file to read the coordinates from.
 *      redundant -> account for Charlie's code-particularity
 *                   that the last 3 bps are the first 3 too.
 *
 * Returns: an array with the coordinates for each base pair
 */
void IO::read_Crd(string crd_File, bool redundant) {
    
    ifstream fin;
    fin.open(crd_File, ios_base::in);
    
    if (fin.is_open()) {
        
        // Line 1: number of DNA base pairs in the file (size of circle +3 for the overlap of the ends)
        fin >> crd.num_BP;
        
        // Line 2 - numBP+something: numbers of atoms in each base pair (in this case always 63, as either a GC base pair  or a CG base pair)
        int i = 0;
        crd.num_Atoms_In_BP = new int[crd.num_BP];
        crd.total_Atoms = 0;
        for (i = 0; i < crd.num_BP; i++) {
            fin >> crd.num_Atoms_In_BP[i];
            crd.total_Atoms += crd.num_Atoms_In_BP[i];
        }
        
        // Lines 65 - : Initial coordinates for each base pair (x1,y1,z1,x2,y2,z2...x63,y63,z63) 10 floats per line, new line before the start of each subsequent base pair, values in angstroms
        crd.ini_BP_Crds = new float[3 * crd.total_Atoms];
        for (i = 0; i < 3 * crd.total_Atoms; i++) {
            fin >> crd.ini_BP_Crds[i];
        }
        
        fin.close();
        
    } else {
        cout << ">>> ERROR: Can not open crd file!" << endl;
        exit(1);
    }
}



/*
 * Read in initial coordinates of tetrads
 */
void IO::read_Initial_Crds(void) {
    
    int i, displacement[crd.num_BP];
    int num_Atoms, start_Index, end_Index;
    int error_Code;
    
    /* next section sets up some book-keeping stuff and does some sanity checking
     * displacement stores the displacement of base pairs. If the 1st pair is 1, then the
     * 2nd pair displacement is 1 + (3 * the number of atoms in pair 1) */
    displacement[0] = 0;
    for (i = 0; i < crd.num_BP; i++) {
        displacement[i+1] = displacement[i] + 3 * crd.num_Atoms_In_BP[i];
    }
    
    for (i = 0; i < prm.num_Tetrads; i++) {
        
        // Sum the number of atoms in 4 BPs in crd.num_Atoms_In_BP
        num_Atoms = crd.num_Atoms_In_BP[i] + crd.num_Atoms_In_BP[i+1] +
            crd.num_Atoms_In_BP[i+2] + crd.num_Atoms_In_BP[i+3];
        // Get the start and end displacement of tetrads in crd.num_Atoms_In_BP
        start_Index = displacement[i];
        end_Index   = displacement[i+4] - 1;

        // Check if all data is matching
        if (num_Atoms != tetrad[i].num_Atoms_In_Tetrad) {
            cout << ">>> ERROR: The number of atoms in crd file is wrong." << endl;
            MPI_Abort(MPI_COMM_WORLD, error_Code);
        } else if ((end_Index - start_Index + 1) != (3 * num_Atoms)) {
            cout << ">>> ERROR: Index goes wrong." << endl;
            MPI_Abort(MPI_COMM_WORLD, error_Code);
        }
        
        // Read in the initial coordinates
        for (int j = 0; j < 3 * tetrad[i].num_Atoms_In_Tetrad; j++) {
            tetrad[i].coordinates[j] = crd.ini_BP_Crds[start_Index++];
        }
    }
}



/*
 * Generate the Gaussian stochastic term. Assuming unitless. 
 */
void EDMD::generate_Stochastic_Term(void) {
    
}

/*
 *
 */
void EDMD::calculate_ED_Forces(Tetrad tetrad) {
    
}

/*
 *
 */
void EDMD::calculate_VDW_Forces(void) {
    
}

/*
 *
 */
void EDMD::update_Velocities(void) {
    
}

/*
 *
 */
void EDMD::update_Coordinates(void) {
    
}





















