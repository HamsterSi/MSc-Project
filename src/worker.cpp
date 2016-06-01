
#include "worker.hpp"

/*
 * The constructor of Master_Management class
 * Function:  Construct Worker Management class
 *
 * Parameter: None
 *
 * Return:    None
 */
Worker_Management::Worker_Management(void) {
    
    comm = MPI_COMM_WORLD;
    
    // Get the rank of worker process
    MPI_Comm_rank(comm, &rank);
    
}


/*
 * The destructor of Master_Management class
 * Function:  Deallocate memory etc.
 *
 * Parameter: None
 *
 * Return:    None
 */
Worker_Management::~Worker_Management(void) {
    
    // Deallocate memory spaces of tetrads
    for (int i = 0; i < num_Tetrads; i++) {
        
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
 * Receive parameters
 * Function:  Workers receive the number of tetrads, number of atoms in every tetrad
 *            and number of evecs from the master process
 *
 * Parameter: None
 *
 * Return:    None
 */
void Worker_Management::parameters_Receiving(void) {
    
    int i, j, signal = 1, parameters[2];
    
    // Receive parameters from the master process
    MPI_Recv(parameters, 2, MPI_INT, 0, TAG_DATA, comm, &status);
    num_Tetrads = parameters[0]; // The number of tetrads
    max_Atoms   = parameters[1]; // The maximum number of atoms in tetrads
    
    int num_Atoms_N_Evecs[2 * num_Tetrads];
    
    // Receive the number of atoms & evecs in tetrads
    MPI_Recv(num_Atoms_N_Evecs, 2 * num_Tetrads, MPI_INT, 0, TAG_DATA, comm, &status);
    
    // Allocate memory for tetrads
    tetrad = new Tetrad[num_Tetrads];
    for (i = 0; i < num_Tetrads; i++) {
        tetrad[i].num_Atoms_In_Tetrad = num_Atoms_N_Evecs[2*i];
        tetrad[i].num_Evecs = num_Atoms_N_Evecs[2*i+1];
        
        tetrad[i].avg_Structure = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
        tetrad[i].masses        = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
        tetrad[i].abq           = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
        tetrad[i].eigenvalues   = new float [tetrad[i].num_Evecs];
        tetrad[i].eigenvectors  = new float*[tetrad[i].num_Evecs];
        for (j = 0; j < tetrad[i].num_Evecs; j++) {
            tetrad[i].eigenvectors[j] = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
        }
        tetrad[i].coordinates   = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
        tetrad[i].velocities    = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
        
        tetrad[i].ED_Forces     = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
        tetrad[i].random_Forces = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
        tetrad[i].NB_Forces     = new float[3 * tetrad[i].num_Atoms_In_Tetrad];
    }
    
    // Send feedback to the master process that has received all parameters.
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_DATA, comm);
    
}


/*
 * Receive tetrads
 * Function:  Workers receive tetrads from the master process
 *
 * Parameter: None
 *
 * Return:    None
 */
void Worker_Management::tetrads_Receiving(void) {
    
    int i, signal = 1;
    MPI_Datatype MPI_Tetrad;
    
    // Receive all tetrads from the master process
    for (i = 0; i < num_Tetrads; i++) {
        
        MPI_Recv(tetrad[i].avg_Structure, 3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+rank+i+1, comm, &status);
        
        MPI_Recv(tetrad[i].masses,        3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+rank+i+2, comm, &status);
        
        MPI_Recv(tetrad[i].abq,           3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+rank+i+3, comm, &status);
        
        MPI_Recv(tetrad[i].eigenvalues,   tetrad[i].num_Evecs,             MPI_FLOAT, 0, TAG_TETRAD+rank+i+4, comm, &status);
        
        for (int j = 0; j < tetrad[i].num_Evecs; j++) {
            MPI_Recv(tetrad[i].eigenvectors[j],  3 * tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+rank+i+5+j, comm, &status);
        }
        
        MPI_Recv(tetrad[i].velocities,    3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+rank+i+6, comm, &status);
        
        MPI_Recv(tetrad[i].coordinates,   3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+rank+i+7, comm, &status);

        /*
        MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &tetrad[i]);
        MPI_Recv(&tetrad[i], 1, MPI_Tetrad, 0, TAG_TETRAD+i, comm, &status);
        MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);*/
    }
    /*
    for (int i = 0; i < 3 * tetrad[89].num_Atoms_In_Tetrad; i++) {
        cout << tetrad[89].coordinates[i] << " "; if ((i+1)%10 == 0) cout << endl;
    }
    cout << endl;*/
    
    // Send feedback to the master process that has received all tetrads
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_TETRAD, comm);
}


/*
 * Function:  Receive tetrad parameters, Compute ED forces of tetrad &
 *            Send the calculated ED forces back to the master.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Worker_Management::ED_Calculation(void) {
    
    int index;
    float ED_Forces[2][3 * max_Atoms + 2];
    
    // Receive the tetrad index from the master process
    MPI_Recv(&index, 1, MPI_INT, 0, TAG_ED, comm, &status);
    
    // Calculate ED forces (ED energy) & random forces
    edmd.calculate_ED_Forces(&tetrad[index], ED_Forces[0], edmd.scaled, 3 * max_Atoms);
    edmd.calculate_Random_Forces(&tetrad[index], ED_Forces[1]);
    
    // Needs to send tetrad index back so the master process can operates on according tetrad
    ED_Forces[0][3 * max_Atoms + 1] = index;
    
    // Send the calculated ED forces, ED energy, random forces & index back
    MPI_Send(&(ED_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, 0, TAG_ED, comm);
    
    /*
    if (index == 89 ) {
        for (int i = 0; i < 3 * tetrad[index].num_Atoms_In_Tetrad; i++) {
            cout << ED_Forces[0][i] << "\t"; if ((i+1)%10 == 0) cout << endl;
        } cout << endl; }*/
    
    //cout << "Rank " << setw(3) << rank << " computed ED forces on Tetrad " << setw(3) << index << endl;
    
}


/*
 * Function:  Receive tetrad parameters, Compute NB forces of tetrad &
 *            Send the calculated NB forces back to the master.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Worker_Management::NB_Calculation(void) {
    
    int indexes[2];
    float NB_Forces[2][3 * max_Atoms + 2];
    
    // Receive the tetrad indexes for NB forces calculation
    MPI_Recv(indexes, 2, MPI_INT, 0, TAG_NB, comm, &status);
    
    // Calculate NB forces, NB energy, etc.
    edmd.calculate_NB_Forces(&tetrad[indexes[0]], &tetrad[indexes[1]], NB_Forces[0], NB_Forces[1], 3 * max_Atoms);
    
    // Need to send both tetrad indexes back to the master process
    NB_Forces[0][3 * max_Atoms + 1] = indexes[0];
    NB_Forces[1][3 * max_Atoms + 1] = indexes[1];
    
    // Send NB forces, energies & indexes back
    MPI_Send(&(NB_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, 0, TAG_NB, comm);
    
    /*
    if (indexes[0] == 83) {
        for (int i = 0 ; i < 2; i++) {
            for (int j = 0 ; j < 3 * max_Atoms+2; j++) {
                cout << NB_Forces[i][j] << " "; if ((i+1)%10 == 0) cout << endl;
            } cout << endl;
        } cout << endl;
    }*/
    
    //cout << "Rank " << setw(3) << rank << " computed NB forces on Tetrad " << setw(3) << indexes[0] << " and " << setw(3) << indexes[1] << endl;
}