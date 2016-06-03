
#include "master.hpp"


/*
 * Function:  The constructor of Master_Management class.
 *            Read-in tetrads parameters, memory allocation and other initialisation
 *
 * Parameter: None
 *
 * Return:    None
 */
Master_Management::Master_Management(void) {
    
    comm     =  MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size); // Get size of MPI processes

}




/*
 * Function:  The destructor of Master_Management class. 
 *            Deallocate memory of arrays.
 *
 * Parameter: None
 *
 * Return:    None
 */
Master_Management::~Master_Management(void) {
    
    delete [] displs;
    delete [] velocities;
    delete [] coordinates;
    
}




/*
 * Function:  Master reads config file which stores the paths of other data files
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::read_Cofig(void) {
    
    io.read_Cofig();
    
    cout << "\nSimulation starting..." << endl;
    cout << "The size of MPI processes: " << size << endl;
    cout << "The iteration number of simulation: " << io.iteration << endl;
    cout << "The shape of DNA: ";
    if (io.circular == true) cout << "Circular" << endl;
    else cout << "Linear" << endl;
    cout << "\nStarting reading data..." << endl;
    
}





/*
 * Function:  Master reads the prm file into tetrads
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::read_Prm(void) {
    
    cout << ">>> Reading prm file..." << endl;
    
    io.read_Prm();
    
    cout << ">>> Read prm file finished" << endl;
    
}





/*
 * Function:  Master reads crd file
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::read_Crd(void) {
    
    cout << ">>> Reading crd file..." << endl;
    
    io.read_Crd();
    
    // Initialise the coordinates (velocities) of tetrads from the crd file
    io.read_Initial_Crds();
    
    cout << ">>> Read crd file finished\nData read completed.\n" << endl;
    
}




/*
 * Function:  Master initialises the simulation
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::initialise(void) {
    
    int i;
    
    cout << "The number of DNA Base Pairs: " << io.crd.num_BP << endl;
    cout << "The number of DNA Tetrads   : " << io.prm.num_Tetrads << endl;
    cout << "The total number of atoms in DNA: " << 3 * io.crd.total_Atoms << endl << endl;
    
    // Pick out the maximum number of atoms in tetrads
    for (max_Atoms = 0, i = 0; i < io.prm.num_Tetrads; i++) {
        max_Atoms = max_Atoms > io.tetrad[i].num_Atoms ? max_Atoms : io.tetrad[i].num_Atoms;
    }
    
    // Allocate memory & initialise arrays that store the whole DNA velocities & coordinates
    velocities  = new float [3 * io.crd.total_Atoms];
    coordinates = new float [3 * io.crd.total_Atoms];
    for (i = 0; i < io.crd.total_Atoms; i++) {
        velocities[i] = coordinates[i] = 0.0;
    }
    
    // Initialise the displacement array & generate diplacements of BP
    displs = new int[io.crd.num_BP];
    for (displs[0] = 0, i = 1; i < io.crd.num_BP; i++) {
        displs[i] = displs[i-1] + 3 * io.crd.num_Atoms_In_BP[i];
    }
    
}





/*
 * Function:  Master sends the number of tetrads, number of atoms in every tetrad
 *            and number of evecs to worker processes
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::send_Parameters(void) {
    
    int parameters[2] = {io.prm.num_Tetrads, max_Atoms};
    int i, signal, *num_Atoms_N_Evecs;
    
    cout << "Starting sending data...\n>>> Master sending parameters..." << endl;
    
    // Stor the number of atoms & evecs in tetrads together (needs to be sent to worker processes)
    num_Atoms_N_Evecs = new int[2 * io.prm.num_Tetrads];
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        num_Atoms_N_Evecs[2 * i] = io.tetrad[i].num_Atoms;
        num_Atoms_N_Evecs[2*i+1] = io.tetrad[i].num_Evecs;
    }
    
    // Send parameters to worker processes
    for (i = 0; i < size-1; i++) {
        MPI_Send(parameters, 2, MPI_INT, i+1, TAG_DATA, comm);
        MPI_Send(num_Atoms_N_Evecs, 2*io.prm.num_Tetrads, MPI_INT, i+1, TAG_DATA, comm);
    }

    // Feedback from worker processes that they have received parameters
    for (int i = 0; i < size-1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_DATA, comm, &status);
    }
    
    cout << ">>> All workers have received parameters" << endl;
    
    delete []num_Atoms_N_Evecs;
}





/*
 * Function:  Master sends tetrads to worker processes
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::send_Tetrads(void) {
    
    int i, j, signal;
    MPI_Datatype MPI_Tetrad;
    
    cout << ">>> Master sending tetrads..." << endl;
    
    // Send all tetrads to worker processes
    for (i = 1; i < size; i++) {
        for (j = 0; j < io.prm.num_Tetrads; j++) {
            
            MPI_Send(io.tetrad[j].avg_Structure, 3*io.tetrad[j].num_Atoms, MPI_FLOAT, i, TAG_TETRAD+i+j+1, comm);
            
            MPI_Send(io.tetrad[j].masses,        3*io.tetrad[j].num_Atoms, MPI_FLOAT, i, TAG_TETRAD+i+j+2, comm);
            
            MPI_Send(io.tetrad[j].abq,           3*io.tetrad[j].num_Atoms, MPI_FLOAT, i, TAG_TETRAD+i+j+3, comm);
            
            MPI_Send(io.tetrad[j].eigenvalues,   io.tetrad[j].num_Evecs,   MPI_FLOAT, i, TAG_TETRAD+i+j+4, comm);
            
            for (int k = 0; k < io.tetrad[j].num_Evecs; k++) {
                MPI_Send(io.tetrad[j].eigenvectors[k],  3*io.tetrad[j].num_Atoms, MPI_FLOAT, i, TAG_TETRAD+i+j+5+k, comm);
            }
            
            MPI_Send(io.tetrad[j].velocities,    3*io.tetrad[j].num_Atoms, MPI_FLOAT, i, TAG_TETRAD+i+j+6, comm);
            
            MPI_Send(io.tetrad[j].coordinates,   3*io.tetrad[j].num_Atoms, MPI_FLOAT, i, TAG_TETRAD+i+j+7, comm);

            /*
            MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &io.tetrad[j]);
            MPI_Send(&io.tetrad[j], 1, MPI_Tetrad, i, TAG_TETRAD+j, comm);
            MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);*/
        }
    }
    
    // Feedback that all worker processes have finished receiving tetrads
    for (int signal, i = 0; i < size-1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_TETRAD, comm, &status);
    }
    
    cout << ">>> All workers have received tetrads\nData sending completed.\n" << endl;
    
}





/*
 * Function:  Master ED forces management. Distrubute ED force calculation among workers,
 *            Main job is to send tetrad parameters to workers to compute ED forces.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::force_Passing(void) {
    
    int i, j, k, flag, index, indexes[2], effective_Pairs;
    int num_Pairs = io.prm.num_Tetrads * (io.prm.num_Tetrads - 1) / 2;
    float max_Forces = 1.0;
    float temp_Forces[2][3 * max_Atoms + 2];
    
    // Generate pair lists of tetrads for NB forces
    int pair_List[num_Pairs][2];
    edmd.generate_Pair_Lists(pair_List, &effective_Pairs, io.prm.num_Tetrads, io.tetrad);
    
    cout << "Calculation starting...\n>>> Calculating ED/NB forces..." << endl;
    cout << ">>> Generating pair lists...\n>>>    Pair list number: " << num_Pairs << endl;
    cout << ">>>    Effective pairs : " << effective_Pairs << endl;
    
    // Send tetrad indexes for ED/NB forces calculation at the beginning
    for (i = 0, j = 0; i < size - 1; i++) {
        if (i < io.prm.num_Tetrads) { // Send tetrad index for ED calculation
            MPI_Send(&i, 1, MPI_INT, i + 1, TAG_ED, comm);
            
        } else {  // i >= num_Tetrads, send tetrad indexes for NB calculation
            while (j < num_Pairs) {
                if (pair_List[j][0] + pair_List[j][1] != -2) {
                    indexes[0] = pair_List[j][0]; indexes[1] = pair_List[j++][1];
                    MPI_Send(indexes, 2, MPI_INT, i, TAG_NB, comm);
                    break;
                } else j++;
            }
        }
    }
    
    // When there are still forces need to be calculated, receive forces back from worker processes & send new tetrad indexes
    while (i < effective_Pairs + io.prm.num_Tetrads + size - 1 && j <= num_Pairs) {

        // Receive ED/NB forces from worker processes
        MPI_Recv(&(temp_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
        
        // Store ED forces & random forces
        if (status.MPI_TAG == TAG_ED) {
            index = (int) temp_Forces[0][3 * max_Atoms + 1];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].ED_Forces[k]     += temp_Forces[0][k];
                io.tetrad[index].random_Forces[k] += temp_Forces[1][k];
            }
            io.tetrad[index].energies[0] = temp_Forces[0][3 * max_Atoms];
            
        // Stroe NB forces
        } else if (status.MPI_TAG == TAG_NB) {
            index = (int) temp_Forces[0][3 * max_Atoms + 1];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].NB_Forces[k] += temp_Forces[0][k];
            }
            io.tetrad[index].energies[1] += temp_Forces[0][3 * max_Atoms];
            io.tetrad[index].energies[2] += temp_Forces[1][3 * max_Atoms];
            index = (int) temp_Forces[1][3 * max_Atoms + 1];
            
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].NB_Forces[k] += temp_Forces[1][k];
            }
            io.tetrad[index].energies[1] += temp_Forces[0][3 * max_Atoms];
            io.tetrad[index].energies[2] += temp_Forces[1][3 * max_Atoms];
            
        }
        
        // If there are some more need to be calculated, send indexes.
        if (i < io.prm.num_Tetrads) { // Send tetrad index for ED calculation
            MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, TAG_ED, comm);
            
        } else {  // i >= num_Tetrads, send tetrad indexes for NB calculation
            while (j < num_Pairs) {
                if (pair_List[j][0] + pair_List[j][1] != -2) {
                    indexes[0] = pair_List[j][0]; indexes[1] = pair_List[j++][1];
                    MPI_Send(indexes, 2, MPI_INT, status.MPI_SOURCE, TAG_NB, comm);
                    break;
                } else j++;
            }
        }
        i++;
        
    }
    
    // Clip NB forces & add random forces into the NB forces
    for (i  = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            io.tetrad[i].NB_Forces[j]  = min( max_Forces, io.tetrad[i].NB_Forces[j]);
            io.tetrad[i].NB_Forces[j]  = max(-max_Forces, io.tetrad[i].NB_Forces[j]);
            io.tetrad[i].NB_Forces[j] += io.tetrad[i].random_Forces[j];
        }
    }

}




/*
 * Function:  Master updates Velocities & Berendsen temperature control
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::cal_Velocities(void) {
    
    int i, j, index;
    
    cout << ">>> Calculating velocities..." << endl;
    
    // Calculate velocities of every tetrad
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Velocities(&io.tetrad[i]);
    }
    
    // v(posX1:posX2) = v(posX1:posX2) + xslice(1:posX2-posX1+1)
    // Aggregate all velocities together for easily wirting out
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            velocities[index++] += io.tetrad[i].velocities[j];
        }
    }
    
    // Needs to consider whether it's circular or linear
    for (i = 0; i < io.crd.num_BP; i++) {
        for (j = 0; j < 3 * io.crd.num_Atoms_In_BP[i]; j++) {
            velocities[j] /= 4;
            if (io.circular == false) {
                if (i < 3) velocities[j] /= (i + 1);
                if (i > io.crd.num_BP - 4) velocities[j] /= (io.crd.num_BP - i);
            }
        }
    }
    
}




/*
 * Function:  Master updates Coordinates
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::cal_Coordinate(void) {
    
    int i, j, index;
    
    cout << ">>> Calculating Coordinates..." << endl;
    
    // Calculate coordinates of all tetrads
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Coordinates(&io.tetrad[i]);
    }
    
    // x(posX1:posX2) = x(posX1:posX2) + xslice(1:posX2-posX1+1)
    // Aggregate all coordinates together for easily wirting out
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            coordinates[index++] += io.tetrad[i].coordinates[j];
        }
    }
    
    // Needs to consider whether it's circular or linear
    for (i = 0; i < io.crd.num_BP; i++) {
        for (j = 0; j < 3 * io.crd.num_Atoms_In_BP[i]; j++) {
            coordinates[j] /= 4;
            if (io.circular == false) {
                if (i < 3) coordinates[j] /= (i + 1);
                if (i > io.crd.num_BP - 4) coordinates[j] /= (io.crd.num_BP - i);
            }
        }
    }
    
    cout << "Calculation completed.\n" << endl;
}





/*
 * Function:  Master writes out energies of all tetrads
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::write_Energy(void) {
    
    int i, j;
    float * energies = new float [4 * io.prm.num_Tetrads];
    
    cout << "Starting writint out...\n>>> Write out energies..." << endl;
    
    // Gather all energies & temperature of tetrads together
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3; j++) {
            energies[i * 4 + j] = io.tetrad[i].energies[j];
        }
        energies[i * 4 + j] = io.tetrad[i].temperature;
    }
    
    // Wrtie out energies
    io.write_Energies(energies);
    
    delete []energies;
}





/*
 * Function:  Master writes out forces of all atoms in DNA
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::write_Forces(void) {
    
    int i, j, index;
    float * ED_Forces     = new float [3 * io.crd.total_Atoms];
    float * random_Forces = new float [3 * io.crd.total_Atoms];
    float * NB_Forces     = new float [3 * io.crd.total_Atoms];
    
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        ED_Forces[i] = random_Forces[i] = NB_Forces[i] = 0.0;
    }
    
    cout << ">>> Write out forces..." << endl;
    
    // Gather all forcees together into three arrays
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            ED_Forces[index]     += io.tetrad[i].ED_Forces[j];
            random_Forces[index] += io.tetrad[i].random_Forces[j];
            NB_Forces[index++]     += io.tetrad[i].NB_Forces[j];
        }
    }
    
    // Write out all forces
    io.write_Forces(ED_Forces, random_Forces, NB_Forces);
    
    delete []ED_Forces;
    delete []random_Forces;
    delete []NB_Forces;
    
}





/*
 * Function:  Master writes out coordinates of all atoms in DNA
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::write_Trajectory(void) {

    cout << ">>> Write out trajectories..." << endl;
    
    io.write_Trajectory(coordinates);
}





/*
 * Function:  Master updates crd file for the next simulation
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::update_Crd_File(void) {
    
    cout << ">>> Update crd file..." << endl;
    
    io.update_Crd(velocities, coordinates);
    
    cout << "Data writing completed.\n" << endl;
}





/*
 * Function:  Master terminates worker processes
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::finalise(void) {
    
    // Send signal to stop all worker processes
    for (int signal = -1, i = 1; i < size; i++) {
        MPI_Send(&signal, 1, MPI_INT, i, TAG_DEATH, comm);
    }
    
    cout << "Simulation ended.\n" << endl;
    
}