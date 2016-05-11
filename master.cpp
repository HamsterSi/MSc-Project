
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
    
    int i, j;
    comm = MPI_COMM_WORLD;
    
    prm_File = "./data//GC90c12.prm";
    crd_File = "./data//GC90_6c.crd";
    
    // Read prm file (Initialise tetrads)
    io.read_Prm(prm_File);
    
    // Read crd file ("true" is for circular DNA and "flase" is for linar DNA)
    io.read_Crd(crd_File, true);
    
    // Read in initial coordinates
    io.read_Initial_Crds();
    cout << "Data read completed." << endl;
    
    // Allocate memory for arrays
    for (max_Atoms = 0, i = 0; i < io.prm.num_Tetrads; i++) {
        max_Atoms = max_Atoms > io.tetrad[i].num_Atoms_In_Tetrad ? max_Atoms : io.tetrad[i].num_Atoms_In_Tetrad;
    }
    ED_Forces = new float[3 * max_Atoms];
    NB_Forces = new float*[2];
    NB_Forces[0] = new float [3 * max_Atoms]; NB_Forces[1] = new float [3 * max_Atoms];
    
    masses = new float[3 * io.crd.total_Atoms];
    displacement = new float[io.crd.num_BP];
    total_ED_Forces   = new float[3 * io.crd.total_Atoms];
    total_NB_Forces   = new float[3 * io.crd.total_Atoms];
    total_Velocities  = new float[3 * io.crd.total_Atoms];
    total_Coordinates = new float[3 * io.crd.total_Atoms];
    
    // Initialise the whole masses array
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms_In_Tetrad; j++) {
            masses[j] = io.tetrad[i].masses[j];
        }
    }
    
    // Generate diplacements of BP
    for (displacement[0] = 0.0, i = 0; i < io.crd.num_BP; i++) {
        displacement[i+1] = displacement[i] + 3 * io.crd.num_Atoms_In_BP[i];
    }
}


/*
 * The destructor of Master_Management class
 * Function:  Write out final result and deallocate arrays
 *
 * Parameter: None
 *
 * Return:    None
 */
Master_Management::~Master_Management(void) {
    
    io.write_Results(output_File, total_Velocities, total_Coordinates);
    
    delete []masses;
    delete []displacement;
    delete []ED_Forces;
    delete []total_ED_Forces;
    delete []NB_Forces[0];
    delete []NB_Forces[1];
    delete []NB_Forces;
    delete []total_NB_Forces;
    delete []total_Velocities;
    delete []total_Coordinates;
    delete []noise_Factor;
    delete []langevin_Forces;
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
void Master_Management::master_ED_Forces(void) {
    
    int i, j, source, dest;
    MPI_Datatype MPI_Tetrad;
    MPI_Status *status;
    
    // Statrt ED forces calculation
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        
        MPI_Library::create_MPI_Tetrad(MPI_Tetrad, io.tetrad[i].num_Atoms_In_Tetrad, io.tetrad[i].num_Evecs);
        MPI_Send(&io.tetrad[i], 1, MPI_Tetrad, dest, TAG_ED, comm);
        
        MPI_Recv(ED_Forces, 3 * io.tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, source, TAG_ED, comm, status);
        
        // Add ED forces to according places int the total_ED_Forces;
        start_Index = displacement[i];
        end_Index   = displacement[i+4] - 1;
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms_In_Tetrad; j++) {
            total_ED_Forces[j] = ED_Forces[start_Index++];
        }
        
        MPI_Library::free_MPI_Tetrad(MPI_Tetrad);
    }
    
    // Calculate the average ED forces
    for (i = 0; i < io.crd.total_Atoms; i++) {
        total_ED_Forces[i] /= 4;
        if (edmd.circular != true) {
            if (i < 3) total_ED_Forces[i] /= (i+1);
            if (i > io.crd.total_Atoms-4) total_ED_Forces[i] /= (io.crd.total_Atoms-i);
        }
    }
}


/*
 * Master NB forces management
 * Function:  Distrubute NB force calculation among workers,
 *            Main job is to send tetrad parameters to workers to compute NB forces.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::master_NB_Forces(void) {
    
    int i, j, source, dest;
    float Energies[2];
    MPI_Datatype MPI_Tetrad;
    MPI_Status *status;
    
    // Generate pair lists of tetrads
    int pair_List[(io.prm.num_Tetrads*(io.prm.num_Tetrads-1)/2)][2];
    //edmd.generate_Pair_Lists(pair_List, io.prm.num_Tetrads, io.tetrad);
    
    // Distrubut NB forces calculation
    for (i = 0; i < io.prm.num_Tetrads*(io.prm.num_Tetrads-1)/2; i++) {
        if (pair_List[i][0] + pair_List[i][1] != -2 ) {
            
            int index_I = pair_List[i][0], index_J = pair_List[i][1];
            Tetrad tetrad[2] = {io.tetrad[index_I], io.tetrad[index_J]};
            MPI_Send(tetrad, 2, MPI_Tetrad, dest, TAG_NB, comm);
            
            MPI_Recv(NB_Forces, 3 * max_Atoms, MPI_FLOAT, 0, TAG_NB, comm, status);
            MPI_Recv(Energies, 2, MPI_FLOAT, 0, TAG_NB, comm, status);
            
            // Add NB forces to according places int the total_ED_Forces;
            start_Index = displacement[pair_List[i][0]];
            end_Index   = displacement[pair_List[i][0]+4] - 1;
            for (j = 0; j < 3 * io.tetrad[pair_List[i][0]].num_Atoms_In_Tetrad; j++) {
                total_NB_Forces[j] = NB_Forces[0][start_Index++];
            }
            
            start_Index = displacement[pair_List[i][1]];
            end_Index   = displacement[pair_List[i][1]+4] - 1;
            for (j = 0; j < 3 * io.tetrad[pair_List[i][1]].num_Atoms_In_Tetrad; j++) {
                total_NB_Forces[j] = NB_Forces[1][start_Index++];
            }
        }
    }
}


/*
 * Master calculates LV forces
 * Function:  Calculate the LV random forces.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::master_LV_Forces(void) {
    int i, j;
    
    // Intialise random number generator
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.generate_Stochastic_Term(i);
    }
    
    // Calculate the randon forces
    {
        noise_Factor    = new float[3 * io.crd.total_Atoms];
        langevin_Forces = new float[3 * io.crd.total_Atoms];
        for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
            noise_Factor[i] = sqrt(2.0 * edmd.gamma * edmd.scaled * masses[i] / edmd.dt);
        }
        //langevin_Forces = random_Forces(io.crd.total_Atoms, noise_Factor);
    }
    
    // Add flv & NB_Forces together -  Langevin dynamics
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        total_NB_Forces[i] += langevin_Forces[i];
    }
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
void Master_Management::master_Total_Forces(void) {
    int i, j;
    
}


/*
 * Master updates velocities
 * Function:  Update Velocities & Berendsen temperature control
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::master_Velocities(void) {
    edmd.update_Velocities(total_Velocities, total_ED_Forces, total_NB_Forces, masses, io.prm.num_Tetrads);
}


/*
 * Master updates coordinates
 * Function:  Update Coordinates
 *
 * Parameter: None
 *
 * Return:    None
 */
void Master_Management::master_Coordinates(void) {
    edmd.update_Coordinates(total_Coordinates, total_Velocities, io.prm.num_Tetrads);
}