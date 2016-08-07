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
 * File:  io.cpp
 * Brief: Implementation of IO class functions for input/output
 */

#include "io.hpp"


IO::IO(void) {

    // Default parameters that can be read from configuartion file
    irest  = 0;
    nsteps = 100000;
    ntsync = 100;
    ntwt   = 100;
    ntpr   = 1000;
    
    prm_File     = "./test/GC90c12.prm";
    crd_File     = "./test/GC90_6c.crd";
    energy_File  = "./data/energies.eng";
    trj_File     = "./data/trajectory.trj";
    new_Crd_File = "./data/crd.crd";
    
}



IO::~IO(void) {
    
    // Deallocate memory of crd
    delete []displs;
    delete []crd.BP_Atoms;
    delete []crd.BP_Crds;
    delete []crd.BP_Vels;
    
    // Deallocate memory spaces for all arrays in tetrads
    for (int i = 0; i < prm.num_Tetrads; i++) {
        tetrad[i].deallocate_Tetrad_Arrays();
    }
    
}



void IO::read_Cofig(EDMD* edmd) {
    
    char line[100] = {0};
    string s1, s2, s3;
    const char * file;
    
    ifstream fin;
    fin.open("./config.txt", ios_base::in);
    
    if (fin.is_open()) {
    
        // Read config data line by line and get the desired data
        for (int i = 1; i < 18; i++) {
            fin.getline(line, sizeof(line));
            stringstream data_Line(line);
            istringstream iss;
            
            switch (i) {
                case  1: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> irest;  break;
                case  2: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> nsteps; break;
                case  3: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> ntsync; break;
                case  4: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> ntwt;   break;
                case  5: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> ntpr;   break;
                case  6: data_Line >> s1 >> s2 >> s3; iss.str(s3);
                    iss >> edmd->dt   ; edmd->dt    *= edmd->constants.timefac; break;
                case  7: data_Line >> s1 >> s2 >> s3; iss.str(s3);
                    iss >> edmd->gamma; edmd->gamma /= edmd->constants.timefac; break;
                case  8: data_Line >> s1 >> s2 >> s3; iss.str(s3);
                    iss >> edmd->tautp; edmd->tautp *= edmd->constants.timefac; break;
                case  9: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> edmd->temperature;
                    edmd->scaled = edmd->constants.Boltzmann * edmd->temperature;   break;
                case 10: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> edmd->mole_Cutoff; break;
                case 11: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> edmd->atom_Cutoff; break;
                case 12: data_Line >> s1 >> s2 >> s3; iss.str(s3); iss >> edmd->mole_Least ; break;
                case 13: data_Line >> s1 >> s2 >> prm_File;     break;
                case 14: data_Line >> s1 >> s2 >> crd_File;     break;
                case 15: data_Line >> s1 >> s2 >> energy_File;  break;
                case 16: data_Line >> s1 >> s2 >> trj_File;     break;
                case 17: data_Line >> s1 >> s2 >> new_Crd_File; break;
            }
        }
        
        fin.close();
        
    } else {
        cout << ">>> ERROR: Can not open the Config file!" << endl;
        exit(1);
    }
    
    // Remove old file before new simualtion starts
    file = energy_File.c_str();  remove(file);
    file = trj_File.c_str();     remove(file);
    file = new_Crd_File.c_str(); remove(file);
    
}



void IO::read_Prm(void) {
    
    int i, j, k;
    ifstream fin;
    fin.open(prm_File.c_str(), ios_base::in);
    
    if (fin.is_open()) {
        
        // Line 1:  Number of (overlapping) tetrads in system (= number of base pairs in a circle)
        fin >> prm.num_Tetrads;
        
        tetrad = new Tetrad[prm.num_Tetrads]; // Initialise the array of tetrads
        
        // Rest of file has data for each tetrad as follows:
        for (i = 0; i < prm.num_Tetrads; i++) {
            
            // Line 2: Number of atoms in the tetrad, and number of eigenvectors
            fin >> tetrad[i].num_Atoms;
            fin >> tetrad[i].num_Evecs;
            
            // Allocate memory spaces for all arrays in tetrads
            tetrad[i].allocate_Tetrad_Arrays();
            
            // Line 3 onwards: Reference (average) structure for the tetrad (x1,y1,z1,x2,y2,z2, etc as in .crd file)
            for (j = 0; j < 3 * tetrad[i].num_Atoms; j++) {
                fin >> tetrad[i].avg[j];
            }
            
            // Line ? onwards: Masses for each atom (amu)
            for (j = 0; j < 3 * tetrad[i].num_Atoms; j += 3) {
                fin >> tetrad[i].masses[j]; // Read & spread masses
                tetrad[i].masses[j + 2] = tetrad[i].masses[j + 1] = tetrad[i].masses[j];
            }
            
            // Line ? onwards: Non-bonded parameters.
            // Each line contains vdW parameters A and B, and partial charge q, for two atoms (e.g. 1st line: atoms 1 and 2, next line: atoms 3 and 4, etc.).
            for (j = 0; j < 3 * tetrad[i].num_Atoms; j++) {
                fin >> tetrad[i].abq[j];
            }
            
            // Next the eigenvector and eigenvalue data for this tetrad:
            // Line ?: The eigenvalue for the 1st eigenvector
            // Line ? onwards: Coefficients of the 1st eigenvector
            // Line ?: The eigenvalue for the 2nd eigenvector
            // Line ? onwards: Coefficients of the 2nd eigenvector
            // (etc to the last eigenvector of this tetrad...)
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



void IO::read_Crd(void) {
    
    int i, j;
    ifstream fin;
    
    // If irest is 0 then open the initial crd file,
    // otherwise open the new crd file which is generated from simualtion.
    if (irest == 0) fin.open(crd_File.c_str(), ios_base::in);
    else fin.open(new_Crd_File.c_str(), ios_base::in);
    
    if (fin.is_open()) {
        
        // Line 1: number of DNA base pairs in the file (size of circle+3 for the overlap of the ends)
        fin >> crd.num_BP;
        
        // Line 2 - numBP + something: numbers of atoms in each base pair (in this case always 63, as either a GC base pair or a CG base pair)
        crd.BP_Atoms = new int[crd.num_BP];
        for (crd.total_Atoms = 0, i = 0; i < crd.num_BP; i++) {
            fin >> crd.BP_Atoms[i];
            crd.total_Atoms += crd.BP_Atoms[i]; // Calculate the total atoms of the DNA
        }
        
        // Allcoate memory for the whole array of velocities & coordinates.
        crd.BP_Crds = new double[3 * crd.total_Atoms];
        crd.BP_Vels = new double[3 * crd.total_Atoms];
        
        generate_Displacements(); // Generate the displacements for DNA base pairs
        
        // Lines 65 - ?: Initial coordinates (& velocities) for each base pair (x1,y1,z1,x2,y2,z2...x63,y63,z63) 10 doubles per line, new line before the start of each subsequent base pair, values in angstroms
        for (i = 0; i < crd.num_BP; i++) {
            for (j = displs[i]; j < displs[i] + 3 * crd.BP_Atoms[i]; j++) {
                fin >> crd.BP_Crds[j];
            }
            if (irest != 0) {
                for (j = displs[i]; j < displs[i] + 3 * crd.BP_Atoms[i]; j++) {
                    fin >> crd.BP_Vels[j];
                }
            }
        }
        
        fin.close();
        
    } else {
        cout << ">>> ERROR: Can not open the crd file!" << endl;
        exit(1);
    }
    
}



void IO::generate_Displacements(void) {

    // Initialise the displacement array
    displs = new int[crd.num_BP + 1];
    
    // Generate diplacements of BP
    displs[0] = 0;
    for (int i = 1; i < crd.num_BP + 1; i++) {
        displs[i] = displs[i - 1] + 3 * crd.BP_Atoms[i - 1];
        
    }
    
}



void IO::initialise_Tetrad_Crds(void) {

    int i, j, error_Code = 0;
    int num_Atoms, start_Index, end_Index;

    for (i = 0; i < prm.num_Tetrads; i++) {
        
        // Sum the number of atoms in 4 base pairs of one tetrad
        num_Atoms = crd.BP_Atoms[i] + crd.BP_Atoms[i + 1] + crd.BP_Atoms[i + 2] + crd.BP_Atoms[i + 3];
        
        // Get the start & end displacement of tetrads
        start_Index = displs[i];
        end_Index   = displs[i + 4] - 1;

        // Check if all data is matching, otherwise quit simulation
        if (num_Atoms != tetrad[i].num_Atoms) {
            cout << ">>> ERROR: The number of atoms in crd file is wrong." << endl;
            MPI_Abort(MPI_COMM_WORLD, error_Code);
        } else if ((end_Index - start_Index + 1) != (3 * num_Atoms)) {
            cout << ">>> ERROR: Index goes wrong." << endl;
            MPI_Abort(MPI_COMM_WORLD, error_Code);
        }
        
        // Read in the initial coordinates & velocities, initialise forces to 0
        for (j = 0; j < 3 * tetrad[i].num_Atoms; j++) {
            if (irest == 0) tetrad[i].velocities[j] = 0.0;
            else tetrad[i].velocities[j] = crd.BP_Vels[start_Index];
            
            tetrad[i].coordinates[j]   = crd.BP_Crds[start_Index++];
            tetrad[i].ED_Forces[j]     = 0.0;
            tetrad[i].random_Forces[j] = 0.0;
            tetrad[i].NB_Forces[j]     = 0.0;
        }
    }
    
}



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



void IO::write_Trajectory(int istep, int index, double* coordinates) {
    
    ofstream fout;
    fout.open(trj_File.c_str(), ios_base::app);
    
    if (fout.is_open()) {
        
        // Write out the total atoms in the DNA
        if (istep == 0) { fout << crd.total_Atoms << endl; }

        // Write out the coordinates
        for (int i = 0; i < index; i++) {
            fout << fixed << setw(10) << setprecision(4) << coordinates[i] << " ";
            if ((i + 1) % 10 == 0) fout << endl;
        }
        
        fout.close();
        
    } else {
        cout << ">>> ERROR: Can not open the trajectory file!" << endl;
        exit(1);
    }
    
}



void IO::update_Crd_File(double* velocities, double* coordinates) {
    
    int i, j;
    ofstream fout;
    fout.open(new_Crd_File.c_str(), ios_base::out);
    
    if (fout.is_open()) {
        
        // Write out the number of base pairs
        fout << crd.num_BP << endl;
        
        // Write out the number of atoms of every base pair
        for (i = 0; i < crd.num_BP; i++) {
            fout << crd.BP_Atoms[i] << endl;
        }
        
        // Write out coordinates & velocities
        for (i = 0; i < crd.num_BP; i++) {
            for (j = displs[i]; j < displs[i] + 3 * crd.BP_Atoms[i]; j++) {
                fout << fixed << setw(10) << setprecision(4) << coordinates[j] << " ";
                if ((j + 1) % 10 == 0) fout << endl;
            }
            fout << endl;
            
            for (j = displs[i]; j < displs[i] + 3 * crd.BP_Atoms[i]; j++) {
                fout << fixed << setw(10) << setprecision(4) << velocities[j] << " ";
                if ((j + 1) % 10 == 0) fout << endl;
            }
            fout << endl;
        }
        
        fout.close();
        
    } else {
        cout << ">>> ERROR: Can not open the new crd file!" << endl;
        exit(1);
    }
    
}



