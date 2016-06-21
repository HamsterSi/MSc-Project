/*
 * The implementation of IO class functions
 */

#include "io.hpp"


/*
 * Function:  The constructor of IO class. Set default parameters.
 *
 * Parameter: None
 *
 * Return:    None
 */
IO::IO(void) {

    irest  = 0;
    nsteps = 100000;
    ntsync = 100;
    ntwt   = 100;
    ntpr   = 1000;
    ncycs  = 1000;
    
    circular  = false;
    
    prm_File     = "./test/GC90c12.prm";
    crd_File     = "./test/GC90_6c.crd";
    energy_File  = "./data/energies.eng";
    forces_File  = "./data/forces.fcs";
    trj_File     = "./data/trajectory.trj";
    new_Crd_File = "./data/crd.crd";
    
}





/*
 * Function:  The destructor of IO class. Deallocate memory of arrays.
 *
 * Parameter: None
 *
 * Return:    None
 */
IO::~IO(void) {
    
    // Deallocate memory of crd
    delete []displs;
    delete []crd.num_BP_Atoms;
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
void IO::read_Cofig(EDMD* edmd) {
    
    char line[100] = {0};
    string s1, s2, s3;
    
    ifstream fin;
    fin.open("./config.txt", ios_base::in);
    
    if (fin.is_open()) {
    
        for (int i = 1; i < 20; i++) {
            fin.getline(line, sizeof(line));
            stringstream file_Path(line);
            istringstream iss;
            
            switch (i) {
                case  1: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> irest;  break;
                case  2: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> nsteps; break;
                case  3: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> ntsync; break;
                case  4: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> ntwt;   break;
                case  5: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> ntpr;   break;
                case  6: file_Path >> s1 >> s2 >> s3; iss.str(s3);
                    iss >> boolalpha >> circular; break;
                case  7: file_Path >> s1 >> s2 >> s3; iss.str(s3);
                    iss >> edmd->dt   ; edmd->dt    *= edmd->constants.timefac; break;
                case  8: file_Path >> s1 >> s2 >> s3; iss.str(s3);
                    iss >> edmd->gamma; edmd->gamma /= edmd->constants.timefac; break;
                case  9: file_Path >> s1 >> s2 >> s3; iss.str(s3);
                    iss >> edmd->tautp; edmd->tautp *= edmd->constants.timefac; break;
                case 10: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> edmd->temperature;
                    edmd->scaled = edmd->constants.Boltzmann * edmd->temperature;   break;
                case 11: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> edmd->mole_Cutoff; break;
                case 12: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> edmd->atom_Cutoff; break;
                case 13: file_Path >> s1 >> s2 >> s3; iss.str(s3); iss >> edmd->mole_Least ; break;
                case 14: file_Path >> s1 >> s2 >> prm_File;     break;
                case 15: file_Path >> s1 >> s2 >> crd_File;     break;
                case 16: file_Path >> s1 >> s2 >> energy_File;  break;
                case 17: file_Path >> s1 >> s2 >> forces_File;  break;
                case 18: file_Path >> s1 >> s2 >> trj_File;     break;
                case 19: file_Path >> s1 >> s2 >> new_Crd_File; break;
            }
        }
        
        fin.close();
        
    } else {
        cout << ">>> ERROR: Can not open the Config file!" << endl;
        exit(1);
    }
    
    const char * file = energy_File.c_str(); remove(file);
    file = forces_File.c_str();  remove(file);
    file = trj_File.c_str();     remove(file);
    file = new_Crd_File.c_str(); remove(file);
    
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
    fin.open(prm_File.c_str(), ios_base::in);
    
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
            tetrad[i].avg_Structure = new double[3 * tetrad[i].num_Atoms];
            tetrad[i].masses        = new double[3 * tetrad[i].num_Atoms];
            tetrad[i].abq           = new double[3 * tetrad[i].num_Atoms];
            
            tetrad[i].eigenvalues   = new double[tetrad[i].num_Evecs];
            tetrad[i].eigenvectors  = new double* [tetrad[i].num_Evecs];
            for (j = 0; j < tetrad[i].num_Evecs; j++) {
                tetrad[i].eigenvectors[j] = new double[3 * tetrad[i].num_Atoms];
            }
            
            tetrad[i].velocities    = new double[3 * tetrad[i].num_Atoms];
            tetrad[i].coordinates   = new double[3 * tetrad[i].num_Atoms];
            
            tetrad[i].ED_Forces     = new double[3 * tetrad[i].num_Atoms];
            tetrad[i].random_Forces = new double[3 * tetrad[i].num_Atoms];
            tetrad[i].NB_Forces     = new double[3 * tetrad[i].num_Atoms];
            
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
    if (irest == 0) fin.open(crd_File.c_str(), ios_base::in);
    else fin.open(new_Crd_File.c_str(), ios_base::in);
    
    if (fin.is_open()) {
        
        // Line 1: number of DNA base pairs in the file (size of circle +3 for the overlap of the ends)
        fin >> crd.num_BP;
        
        // Line 2 - numBP+something: numbers of atoms in each base pair (in this case always 63, as either a GC base pair  or a CG base pair)
        crd.num_BP_Atoms = new int[crd.num_BP];
        crd.total_Atoms = 0;
        
        for (int i = 0; i < crd.num_BP; i++) {
            fin >> crd.num_BP_Atoms[i];
            crd.total_Atoms += crd.num_BP_Atoms[i];
        }
        
        // Lines 65 - : Initial coordinates for each base pair (x1,y1,z1,x2,y2,z2...x63,y63,z63) 10 doubles per line, new line before the start of each subsequent base pair, values in angstroms
        crd.ini_BP_Crds = new double[3 * crd.total_Atoms];
        for (int i = 0; i < 3 * crd.total_Atoms; i++) {
            fin >> crd.ini_BP_Crds[i];
        }
        
        // Lines ? - : If this is the new crd file & not the 1st time to run the simulation, then there are velocities to read
        crd.ini_BP_Vels = new double[3 * crd.total_Atoms];
        if (irest != 0) {
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
 * Function:   Set up some book-keeping stuff. "displs" store the displacements of base pairs.
 *             If 1st BP is 0, then 2nd BP displacement is 0+(3 * number of atoms in 1st BP)
 *
 * Parameters: None.
 *
 * Returns:    None.
 */
void IO::generate_Displacements(void) {

    // Initialise the displacement array & generate diplacements of BP
    displs = new int[crd.num_BP + 1];
    displs[0] = 0;
    for (int i = 1; i < crd.num_BP + 1; i++) {
        displs[i] = displs[i - 1] + 3 * crd.num_BP_Atoms[i - 1];
    }
    
}





/*
 * Function:   Initial coordinates for every tetrad from crd
 *
 * Parameters: None.
 *
 * Returns:    None.
 */
void IO::initialise_Tetrad_Crds(void) {

    int i, j, error_Code = 0;
    int num_Atoms, start_Index, end_Index;

    for (i = 0; i < prm.num_Tetrads; i++) {
        
        // Sum the number of atoms in 4 BPs in crd.num_BP_Atoms
        num_Atoms = crd.num_BP_Atoms[i] + crd.num_BP_Atoms[i + 1] +
            crd.num_BP_Atoms[i + 2] + crd.num_BP_Atoms[i + 3];
        
        // Get the start & end displacement of tetrads in crd.num_BP_Atoms
        start_Index = displs[i];
        end_Index   = displs[i + 4] - 1;

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
            if (irest == 0) tetrad[i].velocities[j] = 0.0;
            else tetrad[i].velocities[j] = crd.ini_BP_Vels[start_Index];
            
            tetrad[i].coordinates[j] = crd.ini_BP_Crds[start_Index++];
            tetrad[i].ED_Forces[j] = 0.0;
            tetrad[i].random_Forces[j] = 0.0;
            tetrad[i].NB_Forces[j] = 0.0;
        }
        
        tetrad[i].ED_Energy = tetrad[i].NB_Energy   = 0.0;
        tetrad[i].EL_Energy = tetrad[i].temperature = 0.0;
    }
    
}





/*
 * Function: Write out all ED forces, random forces & NB forces
 *
 * Parameters: * fout -> The file stream pointer
 *             * data -> The data needs to write out
 *
 * Returns: None.
 */
void IO::write_Template(ofstream* fout, double* data) {
    
    int i, j, index;
    
    for (i = 0; i < crd.num_BP; i++) {
        for (index = displs[i], j = 0; j < 3 * crd.num_BP_Atoms[i]; j++) {
            
            (* fout) << fixed << setw(10) << setprecision(4) << data[index++] << " ";
            if ((j + 1) % 10 == 0) (* fout) << endl;
            
        }
        (* fout) << endl;
    }
    (* fout) << endl;
    
}






/*
 * Function: Write out energies & temperature of tetrads.
 *
 * Parameters: istep      -> The number of iterations
 *             * energies -> Energies & temperature of tetrads
 *
 * Returns: None.
 */
void IO::write_Energies(int istep, double energies[]) {
    
    ofstream fout;
    fout.open(energy_File.c_str(), ios_base::app);
    
    if (fout.is_open()) {
        
        // Write out energies & temperature
        fout << "Iteration, ED Energy, NB_Energy, ELE_Energy, Temperature: ";
        fout << istep << ", ";
        fout << setprecision(8) << energies[0] << ", " << setprecision(8) << energies[1] << ", ";
        fout << setprecision(8) << energies[2] << ", " << setprecision(8) << energies[3] << endl;
        
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
void IO::write_Forces(double* ED_Forces, double* random_Forces, double* NB_Forces) {
    
    ofstream fout;
    fout.open(forces_File.c_str(), ios_base::out);
    
    if (fout.is_open()) {
    
        // Write out ED forces
        fout << "ED Forces:" << endl;
        write_Template(&fout, ED_Forces);
        fout << endl;
        
        // Write out random forces
        fout << "Random Forces:" << endl;
        write_Template(&fout, random_Forces);
        fout << endl;
        
        // Write out NB forces
        fout << "NB Forces:" << endl;
        write_Template(&fout, NB_Forces);
        
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
void IO::write_Trajectory(double* coordinates) {
    
    ofstream fout;
    fout.open(trj_File.c_str(), ios_base::out);
    
    if (fout.is_open()) {
        
        fout << "Coordinates:" << endl;
        write_Template(&fout, coordinates);
        
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
void IO::update_Crd(double* velocities, double* coordinates) {
    
    ofstream fout;
    fout.open(new_Crd_File.c_str(), ios_base::out);
    
    if (fout.is_open()) {
        
        int i, j, index;
        
        // Write out the total number of base pairs
        fout << crd.num_BP << endl;
        
        // Write out the number of atoms in every base pairs
        for (i = 0; i < crd.num_BP; i++) {
            fout << crd.num_BP_Atoms[i] << endl;
        }
        
        // Write out coordinates & velocities
        write_Template(&fout, coordinates);
        write_Template(&fout, velocities );
        
        fout.close();
        
    } else {
        cout << ">>> ERROR: Can not open the new crd file!" << endl;
        exit(1);
    }
    
}



