
#include "master.hpp"

/*
 * The constructor of Master_Management class
 * Function:  Read-in tetrads parameters, memory allocation and other initialisation
 *
 * Parameter: None
 *
 * Return:    None
 */
Master_Management::Master_Management(void) {
    
    prm_File = "./data//GC90c12.prm";
    crd_File = "./data//GC90_6c.crd";
    output_File =  "./result.rst";
    
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);

}


/*
 * The destructor of Master_Management class
 * Function:  Deallocate memory etc.
 *
 * Parameter: None
 *
 * Return:    None
 */
Master_Management::~Master_Management(void) {
    
    delete []displs;
    delete []whole_Velocities;
    delete []whole_Coordinates;
    
}


/*
 * Initialise the simulation
 * Function:  Master reads data from files and allocate memory for arrays;
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::initialise(void) {
    int i, j, k;
    
    cout << "\nInitialising program..." << endl;
    cout << "Size of MPI processes: " << size << endl;
    
    cout << "\nSimulation starting...\n\nData reading starting..." << endl;
    cout << ">>> Reading prm file..." << endl;
    io.read_Prm(prm_File);
    
    cout << ">>> Read prm file finished\n>>> Reading crd file..." << endl;
    io.read_Crd(crd_File, true); // "true" for circular, "flase" for linar
    
    cout << ">>> Read crd file finished" << endl;
    io.read_Initial_Crds();
    
    cout << "Data read completed.\n" << endl;
    cout << "The number of DNA Base Pairs: " << io.crd.num_BP << endl;
    cout << "The number of DNA Tetrads   : " << io.prm.num_Tetrads << "\n" << endl;
    cout << "Total atoms in DNA segment  : " << 3 * io.crd.total_Atoms << endl;
    
    // Pick out the maximum number of atoms in tetrad
    for (max_Atoms = 0, i = 0; i < io.prm.num_Tetrads; i++) {
        max_Atoms = max_Atoms > io.tetrad[i].num_Atoms_In_Tetrad ? max_Atoms : io.tetrad[i].num_Atoms_In_Tetrad;
    }
    
    whole_Velocities  = new float[3 * io.crd.total_Atoms];
    whole_Coordinates = new float[3 * io.crd.total_Atoms];
    
    // Generate diplacements of BP
    displs = new int[io.prm.num_Tetrads];
    for (displs[0] = 0, i = 1; i < io.crd.num_BP; i++) {
        displs[i] = displs[i-1] + 3 * io.crd.num_Atoms_In_BP[i];
    }
    
    /*
    cout << "displs: " << endl;
    for (i = 0; i < io.crd.num_BP; i++) cout << displs[i] << " ";
    cout << endl << endl;*/
}


/*
 * Send parameters
 * Function:  Master send the number of tetrads, number of atoms in every tetrad 
 *            and number of evecs to worker processes
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::parameters_Sending(void) {
    
    int parameters[2] = {io.prm.num_Tetrads, max_Atoms};
    int i, signal, *num_Atoms_N_Evecs;
    
    cout << "Data sending starting...\n>>> Master sending parameters..." << endl;
    
    // Stor the number of atoms & evecs in tetrads together (needs to be sent to worker processes)
    num_Atoms_N_Evecs = new int[2 * io.prm.num_Tetrads];
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        num_Atoms_N_Evecs[2 * i] = io.tetrad[i].num_Atoms_In_Tetrad;
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
 * Send tetrads
 * Function:  Master send tetrads to worker processes
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::tetrads_Sending(void) {
    
    int i, j, signal;
    MPI_Datatype MPI_Tetrad;
    
    cout << ">>> Master sending tetrads..." << endl;
    
    // Send all tetrads to worker processes
    for (i = 1; i < size; i++) {
        for (j = 0; j < io.prm.num_Tetrads; j++) {
            
            MPI_Send(io.tetrad[j].avg_Structure, 3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, i, TAG_TETRAD+i+j+1, comm);
            
            MPI_Send(io.tetrad[j].masses,        3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, i, TAG_TETRAD+i+j+2, comm);
            
            MPI_Send(io.tetrad[j].abq,           3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, i, TAG_TETRAD+i+j+3, comm);
            
            MPI_Send(io.tetrad[j].eigenvalues,   io.tetrad[j].num_Evecs,             MPI_FLOAT, i, TAG_TETRAD+i+j+4, comm);
            
            for (int k = 0; k < io.tetrad[j].num_Evecs; k++) {
                MPI_Send(io.tetrad[j].eigenvectors[k],  3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, i, TAG_TETRAD+i+j+5+k, comm);
            }
            
            MPI_Send(io.tetrad[j].velocities,    3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, i, TAG_TETRAD+i+j+6, comm);
            
            MPI_Send(io.tetrad[j].coordinates,   3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, i, TAG_TETRAD+i+j+7, comm);

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
 * Master ED forces management
 * Function:  Distrubute ED force calculation among workers,
 *            Main job is to send tetrad parameters to workers to compute ED forces.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::force_Calculation(void) {
    
    int i, j, k, flag, index, indexes[2], effective_PL = 0;
    int num_PL = io.prm.num_Tetrads * (io.prm.num_Tetrads - 1) / 2;
    float max_Forces = 1.0;
    float temp_Forces[2][3 * max_Atoms + 2];
    
    // Generate pair lists of tetrads for NB forces
    int pair_List[num_PL][2];
    edmd.generate_Pair_Lists(pair_List, io.prm.num_Tetrads, io.tetrad);
    
    cout << "Forces calculation starting..." << endl;
    cout << ">>> Pair list number: " << num_PL << endl;
    for (i = 0; i < num_PL; i++) {
        if (pair_List[i][0] + pair_List[i][1] != -2) effective_PL++;
    }
    
    // Send tetrad indexes for ED/NB forces calculation at the beginning
    for (i = 0, j = 0; i < size-1; i++) {
        if (i < io.prm.num_Tetrads) { // Send tetrad index for ED calculation
            MPI_Send(&i, 1, MPI_INT, i+1, TAG_ED, comm);
            
        } else {  // i >= num_Tetrads, send tetrad indexes for NB calculation
            while (j < num_PL) {
                if (pair_List[j][0] + pair_List[j][1] != -2) {
                    indexes[0] = pair_List[j][0]; indexes[1] = pair_List[j++][1];
                    MPI_Send(indexes, 2, MPI_INT, status.MPI_SOURCE, TAG_NB, comm);
                    break;
                } else j++;
            }
        }
    }
    
    // When there are still forces need to be calculated, receive forces back from worker processes & send new tetrad indexes
    while (i < effective_PL + io.prm.num_Tetrads + size - 1 && j <= num_PL) {

        // Receive ED/NB forces from worker processes
        MPI_Recv(&(temp_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
        
        // Store ED forces & random forces
        if (status.MPI_TAG == TAG_ED) {
            index = (int) temp_Forces[0][3 * max_Atoms + 1];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms_In_Tetrad; k++) {
                io.tetrad[index].ED_Forces[k]     += temp_Forces[0][k];
                io.tetrad[index].random_Forces[k] += temp_Forces[1][k];
            }
            
        // Stroe NB forces
        } else if (status.MPI_TAG == TAG_NB) {
            index = (int) temp_Forces[0][3 * max_Atoms + 1];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms_In_Tetrad; k++) {
                io.tetrad[index].NB_Forces[k] += temp_Forces[0][k];
            } 
            index = (int) temp_Forces[1][3 * max_Atoms + 1];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms_In_Tetrad; k++) {
                io.tetrad[index].NB_Forces[k] += temp_Forces[1][k];
            }
            
        }
        
        // If there are some more need to be calculated, send indexes.
        if (i < io.prm.num_Tetrads) { // Send tetrad index for ED calculation
            MPI_Send(&i, 1, MPI_INT, i+1, TAG_ED, comm);
            
        } else {  // i >= num_Tetrads, send tetrad indexes for NB calculation
            while (j < num_PL) {
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
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms_In_Tetrad; j++) {
            io.tetrad[i].NB_Forces[j] = min( max_Forces, io.tetrad[i].NB_Forces[j]);
            io.tetrad[i].NB_Forces[j] = max(-max_Forces, io.tetrad[i].NB_Forces[j]);
            io.tetrad[i].NB_Forces[j] += io.tetrad[i].random_Forces[j];
        }
    }
    
    /*
    for (j = 0; j < 3 * io.tetrad[89].num_Atoms_In_Tetrad; j++) {
        cout << io.tetrad[89].ED_Forces[j] << " ";
    } cout << endl;*/
    
    cout << "Forces calculation completed.\n" << endl;
    
}


/*
 * Master updates velocities
 * Function:  Update Velocities & Berendsen temperature control
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::velocity_Calculation(void) {
    
    int i, j, index;
    
    cout << "Velocities calculation starting..." << endl;
    
    // Calculate velocities of every tetrad
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Velocities(&io.tetrad[i]);
    }
    
    // v(posX1:posX2) = v(posX1:posX2) + xslice(1:posX2-posX1+1)
    // Aggregate all velocities together for easily wirting out
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms_In_Tetrad; j++) {
            whole_Velocities[index++] = io.tetrad[i].velocities[j];
        }
    }
    
    // Needs to consider whether it's circular or linear
    
    for (i = 0; i < io.crd.total_Atoms; i++) {
        whole_Velocities[i] *= 0.25;
    }
    
    /*
     for (int i = 0;  i < 3 * io.crd.total_Atoms; i++) {
     cout << whole_Velocities[i] << "\t"; if ((i+1)%10 == 0) cout << endl;
     }
     cout << endl;*/
    
    cout << "Velocities calculation complteted.\n" << endl;
    
}


/*
 * Master updates coordinates
 * Function:  Update Coordinates
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::coordinate_Calculation(void) {
    
    int i, j, index;
    
    cout << "Coordinates calculation starting..." << endl;
    
    // Calculate coordinates of all tetrads
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Coordinates(&io.tetrad[i]);
    }
    
    // x(posX1:posX2) = x(posX1:posX2)+xslice(1:posX2-posX1+1)
    // Aggregate all coordinates together for easily wirting out
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms_In_Tetrad; j++) {
            whole_Coordinates[index++] = io.tetrad[i].coordinates[j];
        }
    }
    
    // Needs to consider whether it's circular or linear
    
    for (i = 0; i < io.crd.total_Atoms; i++) {
        whole_Coordinates[i] *= 0.25;
    }
    
    /*
    for (int i = 0;  i < 3 * io.crd.total_Atoms; i++) {
        cout << whole_Coordinates[i] << "\t"; if ((i+1)%10 == 0) cout << endl;
    }
    cout << endl;*/

    cout << "Coordinates calculation complteted.\n" << endl;
}

/*
 * Master terminates all workers
 * Function:  Terminate worker processes
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
    
    // Write out the results to new file (Format to be discussed)
    io.write_Results(output_File, whole_Velocities, whole_Coordinates);
    
    cout << "Simulation ended.\n" << endl;
    
}