/*
 * The implementation of IO class functions
 */

#include "io.hpp"


/*
 *
 */
IO::IO(void) {
    
    // The number of iteration, default shape of DNA &
    // default file paths, can be changed in "Config" file
    iteration = 0;
    circular  = false;
    
    prm_File     = "./data//GC90c12.prm";
    crd_File     = "./data//GC90_6c.crd";
    energy_File  = "./results/energies.eng";
    forces_File  = "./results/forces.fcs";
    trj_File     = "./results/trajectory.trj";
    new_Crd_File = "./results/crd.crd";
    
}





/*
 *
 */
IO::~IO(void) {
    
    // Deallocate memory of crd
    delete []crd.num_Atoms_In_BP;
    delete []crd.ini_BP_Crds;
    delete []crd.ini_BP_Vels;
    
    // Deallocate memory spaces of tetrads
    for (int i = 0; i < prm.num_Tetrads; i++) {
        
        delete []tetrad[i].avg_Structure;
        delete []tetrad[i].masses;
        delete []tetrad[i].abq;
        delete []tetrad[i].eigenvalues;
        for (int j = 0; j < tetrad[i].num_Evecs; j++) {
            delete []tetrad[i].eigenvectors[j];
        }
        delete []tetrad[i].eigenvectors;
        delete []tetrad[i].velocities;
        delete []tetrad[i].coordinates;
        delete []tetrad[i].ED_Forces;
        delete []tetrad[i].random_Forces;
        delete []tetrad[i].NB_Forces;
    }
    
}





/*
 * Function:   Read-in the config file.
 *
 * Parameters: Config   -> path to the configuration fle.
 *
 * Returns:    None.
 */
void IO::read_Cofig(void) {
    
    char line[100] = {0};
    string s1, s2, s3;
    
    ifstream fin;
    fin.open("./Config.md", ios_base::in);
    
    if (fin.is_open()) {
    
        for (int i = 0; i < 8; i++) {
            fin.getline(line, sizeof(line));
            stringstream file_Path(line);
            
            switch (i) {
                case 0:
                    file_Path >> s1 >> s2 >> s3;
                    iteration = stoi(s3); break;
                case 1:
                    file_Path >> s1 >> s2 >> s3;
                    istringstream(s3) >> boolalpha >> circular; break;
                case 2:
                    file_Path >> s1 >> s2 >> prm_File;     break;
                case 3:
                    file_Path >> s1 >> s2 >> crd_File;     break;
                case 4:
                    file_Path >> s1 >> s2 >> energy_File;  break;
                case 5:
                    file_Path >> s1 >> s2 >> forces_File;  break;
                case 6:
                    file_Path >> s1 >> s2 >> trj_File;     break;
                case 7:
                    file_Path >> s1 >> s2 >> new_Crd_File; break;
            }
        }
        
        fin.close();
        
    } else {
        cout << ">>> ERROR: Can not open the Config file!" << endl;
        exit(1);
    }
    
}





/*
 * Function:   Read-in the prm file.
 *
 * Parameters: None.
 *
 * Returns:    None.
 * N.B. Things like circular structures are already taken into account in the prm-file generation.
 */
void IO::read_Prm(void) {
    
    int i, j, k;

    ifstream fin;
    fin.open(prm_File, ios_base::in);
    
    if (fin.is_open()) {
        
        // Line 1:  Number of (overlapping) tetrads in system (= number of base pairs in a circle)
        fin >> prm.num_Tetrads;
        
        tetrad = new Tetrad[prm.num_Tetrads];
        
        // Rest of file has data for each tetrad as follows:
        for (i = 0; i < prm.num_Tetrads; i++) {
            
            // Line 2: Number of atoms in the tetrad, and number of eigenvectors
            fin >> tetrad[i].num_Atoms;
            fin >> tetrad[i].num_Evecs;
            
            // Allocate memory for arrays in tetrads
            tetrad[i].avg_Structure = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].masses        = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].abq           = new float[3 * tetrad[i].num_Atoms];
            
            tetrad[i].eigenvalues   = new float[tetrad[i].num_Evecs];
            tetrad[i].eigenvectors  = new float* [tetrad[i].num_Evecs];
            for (j = 0; j < tetrad[i].num_Evecs; j++) {
                tetrad[i].eigenvectors[j] = new float[3 * tetrad[i].num_Atoms];
            }
            
            tetrad[i].velocities    = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].coordinates   = new float[3 * tetrad[i].num_Atoms];
            
            tetrad[i].ED_Forces     = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].random_Forces = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].NB_Forces     = new float[3 * tetrad[i].num_Atoms];
            
            tetrad[i].energies[0] = tetrad[i].energies[1] = tetrad[i].energies[2] = 0.0;
            tetrad[i].temperature = 0.0;
            
            // Line 3 onwards: Reference (average) structure for the tetrad (x1,y1,z1,x2,y2,z2, etc as in .crd file)
            for (j = 0; j < 3 * tetrad[i].num_Atoms; j++) {
                fin >> tetrad[i].avg_Structure[j];
            }
            
            // Line ? onwards: Masses for each atom (amu)
            for (j = 0; j < 3 * tetrad[i].num_Atoms; ) {
                fin >> tetrad[i].masses[j]; // Read & spread masses
                tetrad[i].masses[j+2] = tetrad[i].masses[j+1] = tetrad[i].masses[j];
                
                j += 3;
            }
            
            // Line ? onwards: Non-bonded parameters. Each line contains vdW parameters A and B, and partial charge q, for two atoms (e.g. 1st line: atoms 1 and 2, next line: atoms 3 and 4, etc.).
            for (j = 0; j < 3 * tetrad[i].num_Atoms; j++) {
                fin >> tetrad[i].abq[j];
            }
            
            // Next the eigenvector and eigenvalue data for this tetrad:
            // Line ?: The eigenvalue for the 1st (largest) eigenvector
            // Line ? onwards: Coefficients of the 1st eigenvector
            // Line ?: The eigenvalue for the 2nd eigenvector
            // Line ? onwards: Coefficients of the 2nd eigenvector
            // (etc to the last eigenvector of this tetrad)
            for (j = 0; j < tetrad[i].num_Evecs; j++) {
                fin >> tetrad[i].eigenvalues[j];
                for (k = 0; k < 3 * tetrad[i].num_Atoms; k++) {
                    fin >> tetrad[i].eigenvectors[j][k];
                }
            }
        }
        
        fin.close();
        
    } else {
        cout << ">>> ERROR: Can not open the prm file!" << endl;
        exit(1);
    }
    
}





/*
 * Funtion:    Read-in a coordinates (.crd) file.
 *
 * Parameters: redundant -> account for Charlie's code-particularity
 *                          that the last 3 bps are the first 3 too.
 *
 * Returns:    None.
 */
void IO::read_Crd(void) {
    
    ifstream fin;
    if (iteration == 0) fin.open(crd_File, ios_base::in);
    else fin.open(new_Crd_File, ios_base::in);
    
    if (fin.is_open()) {
        
        // Line 1: number of DNA base pairs in the file (size of circle +3 for the overlap of the ends)
        fin >> crd.num_BP;
        
        // Line 2 - numBP+something: numbers of atoms in each base pair (in this case always 63, as either a GC base pair  or a CG base pair)
        crd.num_Atoms_In_BP = new int[crd.num_BP];
        crd.total_Atoms = 0;
        
        for (int i = 0; i < crd.num_BP; i++) {
            fin >> crd.num_Atoms_In_BP[i];
            crd.total_Atoms += crd.num_Atoms_In_BP[i];
        }
        
        // Lines 65 - : Initial coordinates for each base pair (x1,y1,z1,x2,y2,z2...x63,y63,z63) 10 floats per line, new line before the start of each subsequent base pair, values in angstroms
        crd.ini_BP_Crds = new float[3 * crd.total_Atoms];
        for (int i = 0; i < 3 * crd.total_Atoms; i++) {
            fin >> crd.ini_BP_Crds[i];
        }
        
        // Lines ? - : If this is the new crd file & not the 1st time to run the simulation, then there are velocities to read
        crd.ini_BP_Vels = new float[3 * crd.total_Atoms];
        if (iteration != 0) {
            for (int i = 0; i < 3 * crd.total_Atoms; i++) {
                fin >> crd.ini_BP_Vels[i];
            }
        }
        
        fin.close();
        
    } else {
        cout << ">>> ERROR: Can not open the crd file!" << endl;
        exit(1);
    }
    
}





/*
 * Function:   Read-in initial coordinates for every tetrads from crd
 *
 * Parameters: None.
 *
 * Returns:    None.
 */
void IO::read_Initial_Crds(void) {

    int i, j, displs[crd.num_BP+1];
    int num_Atoms, start_Index, end_Index;
    int error_Code;
    
    /* next section sets up some book-keeping stuff and does some sanity checking
     * displacement stores the displacement of base pairs. If the 1st pair is 1, then the
     * 2nd pair displacement is 1 + (3 * the number of atoms in pair 1) */
    for (displs[0] = 0, i = 1; i < crd.num_BP + 1; i++) {
        displs[i] = displs[i-1] + 3 * crd.num_Atoms_In_BP[i-1];
    }

    for (i = 0; i < prm.num_Tetrads; i++) {
        
        // Sum the number of atoms in 4 BPs in crd.num_Atoms_In_BP
        num_Atoms = crd.num_Atoms_In_BP[i] + crd.num_Atoms_In_BP[i+1] +
            crd.num_Atoms_In_BP[i+2] + crd.num_Atoms_In_BP[i+3];
        
        // Get the start and end displacement of tetrads in crd.num_Atoms_In_BP
        start_Index = displs[i];
        end_Index   = displs[i+4] - 1;

        // Check if all data is matching
        if (num_Atoms != tetrad[i].num_Atoms) {
            cout << ">>> ERROR: The number of atoms in crd file is wrong." << endl;
            MPI_Abort(MPI_COMM_WORLD, error_Code);
        } else if ((end_Index - start_Index + 1) != (3 * num_Atoms)) {
            cout << ">>> ERROR: Index goes wrong." << endl;
            MPI_Abort(MPI_COMM_WORLD, error_Code);
        }
        
        // Read in the initial coordinates
        for (j = 0; j < 3 * tetrad[i].num_Atoms; j++) {
            if (iteration == 0) tetrad[i].velocities[j] = 0.0;
            else tetrad[i].velocities[j] = crd.ini_BP_Vels[start_Index];
            
            tetrad[i].coordinates[j] = crd.ini_BP_Crds[start_Index++];
            tetrad[i].ED_Forces[j] = tetrad[i].random_Forces[j] = tetrad[i].NB_Forces[j] = 0.0;
        }
    }
    
}





/*
 * Function: Write out energies & temperature of tetrads.
 *
 * Parameters: * energies -> Energies & temperature of tetrads
 *
 * Returns: None.
 */
void IO::write_Energies(float* energies) {
    
    ofstream fout;
    fout.open(energy_File, ios_base::out);
    
    if (fout.is_open()) {
        
        // Write out energies & temperature
        fout << "Energies:" << endl;
        fout << "ED Energy \tNB_Energy \tElectrostatic_Energy \tTemperature" << endl;
        for (int i = 0; i < 4 * prm.num_Tetrads; ) {
            fout << energies[i]   << "\t\t";
            fout << energies[i+1] << "\t\t";
            fout << energies[i+2] << "\t\t";
            fout << energies[i+3] << endl;
            i += 4;
        }
        fout << endl;
        
        fout.close();
        
    } else {
        cout << ">>> ERROR: Can not open the energy file!" << endl;
        exit(1);
    }

}




/*
 * Function: Write out all ED forces, random forces & NB forces
 *
 * Parameters: * ED_Forces     -> The total ED forces
 *             * random_Forces -> The total random forces
 *             * NB_Forces     -> The total NB forces
 *
 * Returns: None.
 */
void IO::write_Forces(float* ED_Forces, float* random_Forces, float* NB_Forces) {
    
    ofstream fout;
    fout.open(forces_File, ios_base::out);
    
    if (fout.is_open()) {
        
        int i, j;
        
        // Write out ED forces
        fout << "ED Forces:" << endl;
        for (i = 0; i < crd.num_BP; i++) {
            for (j = 0; j < 3 * crd.num_Atoms_In_BP[i]; j++) {
                fout << setw(10) << ED_Forces[i * 3 * crd.num_Atoms_In_BP[i] + j] << " ";
                if ((j + 1) % 10 == 0) fout << endl;
            }
            fout << endl;
        }
        
        // Write out random forces
        fout << "\n\nRandom Forces:" << endl;
        for (i = 0; i < crd.num_BP; i++) {
            for (j = 0; j < 3 * crd.num_Atoms_In_BP[i]; j++) {
                fout << setw(10) << random_Forces[i * 3 * crd.num_Atoms_In_BP[i] + j] << " ";
                if ((j + 1) % 10 == 0) fout << endl;
            }
            fout << endl;
        }
        
        // Write out NB forces
        fout << "\n\nNB Forces:" << endl;
        for (i = 0; i < crd.num_BP; i++) {
            for (j = 0; j < 3 * crd.num_Atoms_In_BP[i]; j++) {
                fout << setw(10) << NB_Forces[i * 3 * crd.num_Atoms_In_BP[i] + j] << " ";
                if ((j + 1) % 10 == 0) fout << endl;
            }
            fout << endl;
        }
        fout << endl;
        
        fout.close();
        
    } else {
        cout << ">>> ERROR: Can not open the forces file!" << endl;
        exit(1);
    }
}





/*
 * Function: Write out velocities & coordinates
 *
 * Parameters: * coordinates -> The whole coordinates of DNA
 *
 * Returns: None.
 */
void IO::write_Trajectory(float* coordinates) {
    
    ofstream fout;
    fout.open(trj_File, ios_base::out);
    
    if (fout.is_open()) {
        
        // Write out coordinates
        fout << "Coordinates:" << endl;
        
        for (int i = 0; i < crd.num_BP; i++) {
            for (int j = 0; j < 3 * crd.num_Atoms_In_BP[i]; j++) {
                fout << setw(10) << coordinates[i * 3 * crd.num_Atoms_In_BP[i] + j] << " ";
                if ((j + 1) % 10 == 0) fout << endl;
            }
            fout << endl;
        }
        fout << endl;
        
        fout.close();
        
    } else {
        cout << ">>> ERROR: Can not open the trajectory file!" << endl;
        exit(1);
    }
    
}



/*
 * Function: Update the crd file for next iteration
 *
 * Parameters: * velocities  -> The whole velocities of DNA
 *             * coordinates -> The whole coordinates of DNA
 *
 * Returns: None.
 */
void IO::update_Crd(float* velocities, float* coordinates) {
    
    ofstream fout;
    fout.open(new_Crd_File, ios_base::out);
    
    if (fout.is_open()) {
        
        int i, j;
        
        // Write out the total number of base pairs
        fout << crd.num_BP << endl;
        
        // Write out the number of atoms in every base pairs
        for (i = 0; i < crd.num_BP; i++) {
            fout << crd.num_Atoms_In_BP[i] << endl;
        }
        
        // Write out velocities
        for (i = 0; i < crd.num_BP; i++) {
            for (j = 0; j < 3 * crd.num_Atoms_In_BP[i]; j++) {
                fout << setw(10) << velocities[i * 3 * crd.num_Atoms_In_BP[i] + j] << " ";
                if ((j + 1) % 10 == 0) fout << endl;
            }
            fout << endl;
        }
        
        // Write out coordinates
        for (i = 0; i < crd.num_BP; i++) {
            for (j = 0; j < 3 * crd.num_Atoms_In_BP[i]; j++) {
                fout << setw(10) << coordinates[i * 3 * crd.num_Atoms_In_BP[i] + j] << " ";
                if ((j + 1) % 10 == 0) fout << endl;
            }
            fout << endl;
        }
        fout << endl;
        
        fout.close();
        
    } else {
        cout << ">>> ERROR: Can not open the new crd file!" << endl;
        exit(1);
    }
    
}






















