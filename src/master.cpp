
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
    output_File =  "./GC90_Result.rst";
    
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
    /*
    delete []masses;
    delete []displacement;
    
    delete []ED_Forces[0];
    delete []ED_Forces[1];
    delete []ED_Forces;
    delete []NB_Forces[0];
    delete []NB_Forces[1];
    delete []NB_Forces;
    
    delete []total_ED_Forces;
    delete []total_NB_Forces;
    delete []langevin_Forces;
    
    delete []total_Velocities;
    delete []total_Coordinates;
     */
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
    cout << "The number of DNA Tetrads: " << io.prm.num_Tetrads << "\n" << endl;
    
    // Allocate memory for arrays
    for (max_Atoms = 0, i = 0; i < io.prm.num_Tetrads; i++) {
        max_Atoms = max_Atoms > io.tetrad[i].num_Atoms_In_Tetrad ? max_Atoms : io.tetrad[i].num_Atoms_In_Tetrad;
    }
    ED_Forces    = new float*[2];
    ED_Forces[0] = new float [3 * max_Atoms + 2];
    ED_Forces[1] = new float [3 * max_Atoms + 2];
    NB_Forces    = new float*[2];
    NB_Forces[0] = new float [3 * max_Atoms + 2];
    NB_Forces[1] = new float [3 * max_Atoms + 2];
    
    masses            = new float[3 * io.crd.total_Atoms];
    displacement      = new float[io.crd.num_BP];
    
    total_ED_Forces   = new float[3 * io.crd.total_Atoms];
    total_NB_Forces   = new float[3 * io.crd.total_Atoms];
    langevin_Forces   = new float[3 * io.crd.total_Atoms];
    
    total_Velocities  = new float[3 * io.crd.total_Atoms];
    total_Coordinates = new float[3 * io.crd.total_Atoms];
    
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        total_ED_Forces[i]  = total_NB_Forces[i]   = 0.0;
        langevin_Forces[i]  = 0.0;
        total_Velocities[i] = total_Coordinates[i] = 0.0;
    }
    
    // Initialise the whole masses array
    for (i = 0, j = 0; j < io.prm.num_Tetrads - 1; j++) {
        for (k = 0; k < 3 * io.crd.num_Atoms_In_BP[j]; k++) {
            masses[i++] = io.tetrad[j].masses[k];
        }
    }
    for (k = 0; k < 3 * io.tetrad[j].num_Atoms_In_Tetrad; k++) {
        masses[i++] = io.tetrad[j].masses[k];
    }
    /*
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        cout << masses[i] << "\t";
        if((i+1) % 10 == 0) cout << endl;
    }
    cout << endl;*/
    
    // Generate diplacements of BP
    for (displacement[0] = 0.0, i = 0; i < io.crd.num_BP; i++) {
        displacement[i+1] = displacement[i] + 3 * io.crd.num_Atoms_In_BP[i];
    }
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
    int i, *num_Atoms_N_Evecs;
    
    num_Atoms_N_Evecs = new int[2 * io.prm.num_Tetrads];
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        num_Atoms_N_Evecs[2 * i] = io.tetrad[i].num_Atoms_In_Tetrad;
        num_Atoms_N_Evecs[2*i+1] = io.tetrad[i].num_Evecs;
    }
    
    for (i = 0; i < size-1; i++) {
        MPI_Send(parameters, 2, MPI_INT, i+1, TAG_DATA, comm);
        MPI_Send(num_Atoms_N_Evecs, 2*io.prm.num_Tetrads, MPI_INT, i+1, TAG_DATA, comm);
    }
    
    delete []num_Atoms_N_Evecs;
    cout << "Data sending starting...\n>>> Master sending parameters..." << endl;
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
    
    int signal, flag;
    MPI_Datatype MPI_Tetrad;
    MPI_Status status;
    
    for (int i = 0; i < size-1; i++) {
        
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_DATA, comm, &status);
        for (int j = 0; j < io.prm.num_Tetrads; j++) {
            
            MPI_Send(&(io.tetrad[j].avg_Structure[0]), 3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, status.MPI_SOURCE, TAG_TETRAD+j+1, comm);
            
            MPI_Send(&(io.tetrad[j].masses[0]),        3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, status.MPI_SOURCE, TAG_TETRAD+j+2, comm);
            
            MPI_Send(&(io.tetrad[j].abq[0]),           3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, status.MPI_SOURCE, TAG_TETRAD+j+3, comm);
            
            MPI_Send(&(io.tetrad[j].eigenvalues[0]),   io.tetrad[j].num_Evecs,             MPI_FLOAT, status.MPI_SOURCE, TAG_TETRAD+j+4, comm);
            
            MPI_Send(&(io.tetrad[j].eigenvectors[0][0]),  io.tetrad[j].num_Evecs*3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, status.MPI_SOURCE, TAG_TETRAD+j+5, comm);
            
            MPI_Send(&(io.tetrad[j].coordinates[0]),   3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, status.MPI_SOURCE, TAG_TETRAD+j+6, comm);
            
            MPI_Send(&(io.tetrad[j].velocities[0]),    3*io.tetrad[j].num_Atoms_In_Tetrad, MPI_FLOAT, status.MPI_SOURCE, TAG_TETRAD+j+7, comm);
            
            /*
            MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &io.tetrad[j]);
            MPI_Send(&io.tetrad[j], 1, MPI_Tetrad, status.MPI_SOURCE, TAG_TETRAD+j, comm);
            MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);*/
        }
    }
    
    cout << ">>> All workers have received parameters\n>>> Master sending tetrads..." << endl;
}


/*
 * Master receives tetrad received signal
 * Function:  Master will receive the signal from worker processes that 
 *            they have received all tetrads.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::tetrad_Received_Signal(void) {
    
    MPI_Status status;
    
    for (int signal, i = 0; i < size-1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_DATA, comm, &status);
    }
    cout << ">>> All workers have received tetrads\nData sending completed.\n" << endl;
    cout << "Forces calculation starting..." << endl;
    
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
void Master_Management::force_Passing(void) {
    
    int i, j, k, flag, num_PL;
    int tetrad_Index, pair_Index, index, indexes[3];
    MPI_Status status;

    // Generate pair lists of tetrads for NB forces
    int pair_List[(io.prm.num_Tetrads*(io.prm.num_Tetrads-1)/2)][2];
    edmd.generate_Pair_Lists(pair_List, io.prm.num_Tetrads, io.tetrad);
    cout << ">>> Pair list number: " << (io.prm.num_Tetrads*(io.prm.num_Tetrads-1)/2) << endl;
    for (num_PL = 0, i = 0; i < (io.prm.num_Tetrads*(io.prm.num_Tetrads-1)/2); i++) {
        if (pair_List[i][0] + pair_List[i][1] != -2) {
            num_PL++;
        }
    }
    
    for (i = 0, j = 0; i < size-1; i++) { // For every worker process
        
        if (i < io.prm.num_Tetrads) { // Send tetrad index for ED calculation
            MPI_Send(&i, 1, MPI_INT, i+1, TAG_ED, comm);
            
        } else { // i >= num_Tetrads, send tetrad indexes for NB calculation
            while (j < io.prm.num_Tetrads*(io.prm.num_Tetrads-1)/2) {
                if (pair_List[j][0] + pair_List[j][1] != -2) {
                    indexes[0] = j; indexes[1] = pair_List[j][0]; indexes[2] = pair_List[j][1];
                    MPI_Send(indexes, 2, MPI_INT, status.MPI_SOURCE, TAG_NB, comm);
                    j++; break;
                } else { j++; }
            }
        }
    }
    
    while (1)
    {
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &flag, &status);
        if (flag)
        {
            if (status.MPI_TAG == TAG_ED) {
                
                MPI_Recv(&(ED_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, MPI_ANY_SOURCE, TAG_ED, comm, &status);
                
                tetrad_Index = (int) ED_Forces[0][3 * max_Atoms + 1];
                index = displacement[tetrad_Index];
                for (k = 0; k < 3 * io.tetrad[tetrad_Index].num_Atoms_In_Tetrad; k++) {
                    total_ED_Forces[index++] += ED_Forces[0][k];
                    langevin_Forces[index++] += ED_Forces[1][k];
                }
                
            } else if (status.MPI_TAG == TAG_NB) {
                
                MPI_Recv(&(NB_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, MPI_ANY_SOURCE, TAG_NB, comm, &status);
                
                pair_Index   = (int) NB_Forces[0][(int) 3 * max_Atoms + 1];
                
                tetrad_Index = pair_List[pair_Index][0];
                index = displacement[tetrad_Index];
                for (k = 0; k < 3 * io.tetrad[tetrad_Index].num_Atoms_In_Tetrad; k++) {
                    total_NB_Forces[index++] += NB_Forces[0][k];
                }
                
                tetrad_Index = pair_List[pair_Index][1];
                index = displacement[tetrad_Index];
                for (k = 0; k < 3 * io.tetrad[tetrad_Index].num_Atoms_In_Tetrad; k++) {
                    total_NB_Forces[index++] += NB_Forces[1][k];
                }
                
            }
            
            if (i < io.prm.num_Tetrads) { // receive and send ED forces parameters
                MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, TAG_ED, comm);
                
            } else {
                while (j < io.prm.num_Tetrads*(io.prm.num_Tetrads-1)/2) {
                    if (pair_List[j][0] + pair_List[j][1] != -2) {
                        indexes[0] = j; indexes[1] = pair_List[j][0]; indexes[2] = pair_List[j][1];
                        MPI_Send(indexes, 2, MPI_INT, status.MPI_SOURCE, TAG_NB, comm);
                        j++; break;
                    } else { j++; }
                }
            }
            i++;
            
            if (i >= num_PL + io.prm.num_Tetrads + size - 1 && j >= io.prm.num_Tetrads*(io.prm.num_Tetrads-1)/2) {
                break;
            }
        }
    }
    cout << "Forces calculation completed.\n" << endl;
}


/*
 * Master total forces management
 * Function:  Receive all ED & NB forces from workers, 
 *            update the total forces on DNA segement.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::calculate_Total_Forces(void) {
    
    // Calculate the average ED forces
    for (int i = 0; i < 3 * io.crd.total_Atoms; i++) {
        total_ED_Forces[i] /= 4;
        if (edmd.circular == false) {
            if (i < 3) total_ED_Forces[i] /= (i+1);
            if (i > io.crd.total_Atoms-4) total_ED_Forces[i] /= (io.crd.total_Atoms-i);
        }
        
        //cout << total_ED_Forces[i] << "\t";
        //if ((i+1)%10 == 0) cout << endl;
    }
    //cout << endl;
    
    // Add flv & NB_Forces together -  Langevin dynamics
    for (int i = 0; i < 3 * io.crd.total_Atoms; i++) {
        total_NB_Forces[i] += langevin_Forces[i];
   
        //cout << total_ED_Forces[i] << "\t";
        //if ((i+1)%10 == 0) cout << endl;
    }
    //cout << endl;
    
}


/*
 * Master updates velocities
 * Function:  Update Velocities & Berendsen temperature control
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::velocities(void) {
    
    cout << "Velocities calculation starting..." << endl;
    
    edmd.update_Velocities(total_Velocities, total_ED_Forces, total_NB_Forces, masses, io.crd.total_Atoms);
    
    for (int i = 0;  i < 3 * io.crd.total_Atoms; i++) {
        //cout << total_Velocities[i] << "\t";
        //if ((i+1)%10 == 0) cout << endl;
    }
    //cout << endl;
    
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
void Master_Management::coordinates(void) {
    
    cout << "Coordinates calculation starting..." << endl;
    
    edmd.update_Coordinates(total_Coordinates, total_Velocities, io.crd.total_Atoms);

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
    
    io.write_Results(output_File, total_Velocities, total_Coordinates);
    
    for (int stop = -1, i = 1; i < size; i++) {
        MPI_Send(&stop, 1, MPI_INT, i, TAG_DEATH, comm);
    }
    
    cout << "Simulation ended.\n" << endl;
    
}