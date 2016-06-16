
#include "master.hpp"


/*
 * Function:  The constructor of Master class.
 *            Read-in tetrads parameters, memory allocation and other initialisation
 *
 * Parameter: None
 *
 * Return:    None
 */
Master::Master(void) {
    
    max_Atoms = 0;
    
    comm      = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size); // Get size of MPI processes

}




/*
 * Function:  The destructor of Master class. 
 *            Deallocate memory of arrays.
 *
 * Parameter: None
 *
 * Return:    None
 */
Master::~Master(void) {
    
    delete [] velocities;
    delete [] coordinates;
    
}




/*
 * Function:  Master initialises the simulation
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::initialise(void) {
    
    io.read_Cofig(&edmd);

    io.read_Prm();
    
    io.read_Crd();
    
    // Initialise the displacement array & generate diplacements of BP
    io.generate_Displacements();
    
    // Initialise coordinates (velocities) of tetrads from crd file
    io.initialise_Tetrad_Crds();
    
    // Pick out the maximum number of atoms in tetrads
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        max_Atoms = max_Atoms > io.tetrad[i].num_Atoms ? max_Atoms : io.tetrad[i].num_Atoms;
    }
    
    // Allocate memory for storing velocities & coordinates of the whole DNA
    velocities  = new double [3 * io.crd.total_Atoms];
    coordinates = new double [3 * io.crd.total_Atoms];
    
    cout << endl << "Simulation starting..." << endl;
    cout << ">>> MPI Processes: " << size << endl;
    cout << ">>> DNA Shape: ";
    if (io.circular == true) cout << "Circular" << endl;
    else cout << "Linear" << endl;
    
    cout << "Reading prm & crd file...\nData reading completed." << endl << endl;
    cout << "The number of DNA Base Pairs: " << io.crd.num_BP << endl;
    cout << "The number of DNA Tetrads   : " << io.prm.num_Tetrads << endl;
    cout << "Total number of atoms in DNA: " << 3 * io.crd.total_Atoms << endl;
    
}




/*
 * Function:  Master sends the number of tetrads, number of atoms in every tetrad
 *            and number of evecs to worker processes
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::send_Parameters(void) {
    
    double edmd_Para[8] = {edmd.dt, edmd.gamma, edmd.tautp, edmd.temperature,
           edmd.scaled, edmd.mole_Cutoff, edmd.atom_Cutoff, edmd.mole_Least};
    int parameters[3] = {io.prm.num_Tetrads, max_Atoms, io.iteration};
    int i, signal, * tetrad_Para = new int[2 * io.prm.num_Tetrads];
    
    // Send parameters to worker processes
    for (i = 0; i < size - 1; i++) {
        MPI_Send(parameters, 3, MPI_INT, i + 1, TAG_DATA, comm);
    }
    
    // If it is the 1st iteration then need to send the number of atoms & evecs of tetrads
    if (io.iteration == 0) {
        
        // Assign the number of atoms & evecs of all tetrads into array
        for (i = 0; i < io.prm.num_Tetrads; i++) {
            tetrad_Para[2 * i] = io.tetrad[i].num_Atoms;
            tetrad_Para[2*i+1] = io.tetrad[i].num_Evecs;
        }
        
        // Send the edmd parameters & tetrad parameters to workers
        for (i = 0; i < size - 1; i++) {
            MPI_Send(edmd_Para, 8, MPI_DOUBLE, i + 1, TAG_DATA, comm);
            MPI_Send(tetrad_Para, 2 * io.prm.num_Tetrads, MPI_INT, i + 1, TAG_DATA, comm);
        }
    }

    // Feedback from worker processes that they have received parameters
    for (i = 0; i < size - 1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_DATA, comm, &status);
    }
    
    delete [] tetrad_Para;
}




/*
 * Function:  Master sends tetrads to worker processes
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::send_Tetrads(void) {
    
    int i, j, signal;
    MPI_Datatype MPI_Tetrad;
    
    // Send tetrads to worker. The 1st iteration needs to send all parameters
    // of tetrads to workers, but the following iterations only needs to send
    // the velocities & coordinates of tetrads to workers.
    for (i = 1; i < size; i++) {
        for (j = 0; j < io.prm.num_Tetrads; j++) {
            
            if (io.iteration == 0) {
                
                // Create MPI_Datatype for every Tetrad instance & send them to workers
                MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &io.tetrad[j]);
                MPI_Send(&io.tetrad[j], 1, MPI_Tetrad, i, TAG_TETRAD+i+j, comm);
                MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);
                
            } else {
                
                // Only send velocities & coordinates to workers.
                MPI_Send(io.tetrad[j].velocities,  3 * io.tetrad[j].num_Atoms, MPI_DOUBLE, i, TAG_TETRAD+i+j+1, comm);
                MPI_Send(io.tetrad[j].coordinates, 3 * io.tetrad[j].num_Atoms, MPI_DOUBLE, i, TAG_TETRAD+i+j+2, comm);
            }
            
        }
    }
    
    // Feedback that all worker processes have finished receiving tetrads
    for (i = 0; i < size - 1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_TETRAD, comm, &status);
    }
    
}





/*
 * Function:  Find pair list that needs to calculate Nb forces & send them to workers
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::send_Worker_Pairlists(int* j, int num_Pairs, int source, int pair_List[][2]) {
    
    while ((*j) < num_Pairs) {
        
        if (pair_List[(*j)][0] + pair_List[(*j)][1] != -2) {
            
            int indexes[2] = { pair_List[(*j)][0], pair_List[(*j)][1] };
            MPI_Send(indexes, 2, MPI_INT, source, TAG_NB, comm);
            (*j)++; break;
            
        } else { (*j)++; }
    }
}





/*
 * Function:  Master ED forces management. Distrubute ED force calculation among workers,
 *            Main job is to send tetrad parameters to workers to compute ED forces.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::force_Calculation(void) {
    
    int i, j, k, flag, index, effective_Pairs;
    int num_Pairs = io.prm.num_Tetrads * (io.prm.num_Tetrads - 1) / 2;
    double max_Forces = 1.0, temp_Forces[2][3 * max_Atoms + 2];
    
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            io.tetrad[i].ED_Forces[j]     = 0.0;
            io.tetrad[i].random_Forces[j] = 0.0;
            io.tetrad[i].NB_Forces[j]     = 0.0;
        }
        io.tetrad[i].ED_Energy = 0.0;
        io.tetrad[i].NB_Energy = 0.0;
        io.tetrad[i].EL_Energy = 0.0;
    }
    
    // Generate pair lists of tetrads for NB forces
    int pair_List[num_Pairs][2];
    edmd.generate_Pair_Lists(pair_List, &effective_Pairs, io.prm.num_Tetrads, io.tetrad);
    cout << ">>> Generating pair lists: Total Pairlists: " << num_Pairs << ", Effective pairs: " << effective_Pairs << endl;
    
    // Send tetrad indexes for ED/NB forces calculation at the beginning
    for (i = 0, j = 0; i < size - 1; i++) {
        if (i < io.prm.num_Tetrads) {
            // Send tetrad index for ED calculation & calculate random forces
            MPI_Send(&i, 1, MPI_INT, i + 1, TAG_ED, comm);
            edmd.calculate_Random_Forces(&io.tetrad[i]);
        } else {
            // i >= num_Tetrads, send tetrad indexes for NB calculation
            send_Worker_Pairlists(&j, num_Pairs, i, pair_List);
        }
    }
    
    // When there are still forces waiting for calculating,
    // receive ED/NB forces from workers & send new tetrad indexes to workers
    while (i < effective_Pairs + io.prm.num_Tetrads + size - 1 && j <= num_Pairs) {

        // Receive ED/NB forces from workers
        MPI_Recv(&(temp_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
        
        // If the MPI Tag indicates it's ED forces then store ED forces & random forces
        if (status.MPI_TAG == TAG_ED) {
            index = (int) temp_Forces[0][3 * max_Atoms + 1];
            io.tetrad[index].ED_Energy = temp_Forces[0][3 * max_Atoms];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].ED_Forces[k]   = temp_Forces[0][k];
                io.tetrad[index].coordinates[k] = temp_Forces[1][k];
            }
            
        // The MPI Tag shows it's NB forces, stroe NB forces in tetrads
        } else if (status.MPI_TAG == TAG_NB) {
            index = (int) temp_Forces[0][3 * max_Atoms + 1];
            io.tetrad[index].NB_Energy += temp_Forces[0][3 * max_Atoms];
            io.tetrad[index].EL_Energy += temp_Forces[1][3 * max_Atoms];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].NB_Forces[k] += temp_Forces[0][k];
            }
            
            index = (int) temp_Forces[1][3 * max_Atoms + 1];
            io.tetrad[index].NB_Energy += temp_Forces[0][3 * max_Atoms];
            io.tetrad[index].EL_Energy += temp_Forces[1][3 * max_Atoms];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].NB_Forces[k] += temp_Forces[1][k];
            }
            
        }
        
        // If there are some more need to be calculated, send indexes.
        if (i < io.prm.num_Tetrads) {
            // Send tetrad index for ED calculation & calculate random forces
            MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, TAG_ED, comm);
            edmd.calculate_Random_Forces(&io.tetrad[i]);
        } else {
            // i >= num_Tetrads, send tetrad indexes for NB calculation
            send_Worker_Pairlists(&j, num_Pairs, status.MPI_SOURCE, pair_List);
        }
        
        i++;
        
    }
    
    // Clip NB forces & add random forces into the NB forces
    for (i  = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            io.tetrad[i].NB_Forces[j]  = min( max_Forces, io.tetrad[i].NB_Forces[j]);
            io.tetrad[i].NB_Forces[j]  = max(-max_Forces, io.tetrad[i].NB_Forces[j]);
        }
    }
    
}





/*
 * Function:  Master calculates velocities in tetrads
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::cal_Velocities(void) {
    
    // Calculate velocities of every tetrad
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Velocities(&io.tetrad[i], i);
    }
    
}





/*
 * Function:  Master calculates coordinates in tetrads
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::cal_Coordinate(void) {
    
    // Calculate coordinates of all tetrads
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Coordinates(&io.tetrad[i]);
    }
    
}





/*
 * Function:  Master updates the whole velocities/coordinates of DNA
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::data_Processing(void) {
    
    int i, j, index;
    
    // Initialise velocities & coordinates array
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        velocities[i] = coordinates[i] = 0.0;
    }
    
    // Gather all velocities & coordinates into a single array
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = io.displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; index++, j++) {
            velocities [index] += io.tetrad[i].velocities [j];
            coordinates[index] += io.tetrad[i].coordinates[j];
        }
    }
    
    // Process the beginning & end tetrads
    index = io.displs[io.crd.num_BP - 3];
    for (i = 0; index < io.displs[io.crd.num_BP]; index++, i++) {
        velocities [index] += velocities [i];
        coordinates[index] += coordinates[i];
        velocities [i] = velocities [index];
        coordinates[i] = coordinates[index];
    }

    // Divide velocities & coordinates by 4
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        velocities[i] *= 0.25; coordinates[i] *= 0.25;
    }

    // Restore the velocities & coordinates back to tetrads
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = io.displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; index++, j++) {
            io.tetrad[i].velocities [j] = velocities [index];
            io.tetrad[i].coordinates[j] = coordinates[index];
        }
    }
    
}





/*
 * Function:  Master writes out energies of all tetrads
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::write_Energy(void) {
    
    double energies[4] = {0.0};
    
    // Gather energies & temperature of tetrads together
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        energies[0] += io.tetrad[i].ED_Energy;
        energies[1] += io.tetrad[i].NB_Energy;
        energies[2] += io.tetrad[i].EL_Energy;
        energies[3] += io.tetrad[i].temperature;
    }
    
    // Calculate the average temperature of tetrads
    energies[3] /= io.prm.num_Tetrads;
    
    // Wrtie out energies
    io.write_Energies(energies);
    
}





/*
 * Function:  Master writes out forces of all atoms in DNA
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::write_Forces(void) {
    
    int i, j, index;
    double * ED_Forces     = new double [3 * io.crd.total_Atoms];
    double * random_Forces = new double [3 * io.crd.total_Atoms];
    double * NB_Forces     = new double [3 * io.crd.total_Atoms];
    
    // Initialise arrays
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        ED_Forces[i] = random_Forces[i] = NB_Forces[i] = 0.0;
    }
    
    // Gather all forcees together into three arrays
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = io.displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; index++, j++) {
            ED_Forces[index]     += io.tetrad[i].ED_Forces[j];
            random_Forces[index] += io.tetrad[i].random_Forces[j];
            NB_Forces[index]     += io.tetrad[i].NB_Forces[j];
        }
    }
    
    // Write out all forces
    io.write_Forces(ED_Forces, random_Forces, NB_Forces);
    
    delete []ED_Forces;
    delete []random_Forces;
    delete []NB_Forces;
    
}




/*
 * Function:  Master writes file according to the frequency.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::write_Files(void) {
    
    write_Energy();
    
    write_Forces();
    
    io.write_Trajectory(coordinates);
    
    io.update_Crd(velocities, coordinates);
    
}




/*
 * Function:  Master terminates worker processes
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master::finalise(void) {
    
    cout << "Writing Energies out at     " << io.energy_File << endl;
    cout << "Writing Forces   out at     " << io.forces_File << endl;
    cout << "Writing Trajectories out at " << io.trj_File << endl;
    cout << "Writing New Crd  out at     " << io.new_Crd_File << endl << endl;

    // Send signal to stop all worker processes
    for (int signal = -1, i = 1; i < size; i++) {
        MPI_Send(&signal, 1, MPI_INT, i, TAG_DEATH, comm);
    }
    
    cout << "Simulation ended.\n" << endl;
    
}

