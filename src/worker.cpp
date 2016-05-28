
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
    /*
    delete []ED_Forces[0];
    delete []ED_Forces[1];
    delete []ED_Forces;
    delete []NB_Forces[0];
    delete []NB_Forces[1];
    delete []NB_Forces;
    delete []num_Atoms_N_Evecs;
     */
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
    int signal = 1, parameters[2];
    
    MPI_Recv(parameters, 2, MPI_INT, 0, TAG_DATA, comm, &status);
    
    num_Tetrads = parameters[0];
    max_Atoms   = parameters[1];
    num_Atoms_N_Evecs = new int[2 * num_Tetrads];
    
    ED_Forces    = new float*[2];
    ED_Forces[0] = new float [3 * max_Atoms + 2];
    ED_Forces[1] = new float [3 * max_Atoms + 2];
    NB_Forces    = new float*[2];
    NB_Forces[0] = new float [3 * max_Atoms + 2];
    NB_Forces[1] = new float [3 * max_Atoms + 2];
    
    MPI_Recv(num_Atoms_N_Evecs, 2 * num_Tetrads, MPI_INT, 0, TAG_DATA, comm, &status);
    
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
    
    tetrad = new Tetrad[num_Tetrads];
    for (i = 0; i < num_Tetrads; i++) {  // Allocate memory
        tetrad[i].num_Atoms_In_Tetrad = num_Atoms_N_Evecs[2*i];
        tetrad[i].num_Evecs = num_Atoms_N_Evecs[2*i+1];
    
        tetrad[i].avg_Structure = new float[3 * num_Atoms_N_Evecs[2*i]];
        tetrad[i].masses        = new float[3 * num_Atoms_N_Evecs[2*i]];
        tetrad[i].abq           = new float[3 * num_Atoms_N_Evecs[2*i]];
        tetrad[i].eigenvalues   = new float [num_Atoms_N_Evecs[2*i+1]];
        tetrad[i].eigenvectors  = new float*[num_Atoms_N_Evecs[2*i+1]];
        for (int j = 0; j < num_Atoms_N_Evecs[2*i+1]; j++) {
            tetrad[i].eigenvectors[j] = new float[3 * num_Atoms_N_Evecs[2*i]];
        }
        
        tetrad[i].coordinates   = new float[3 * num_Atoms_N_Evecs[2*i]];
        tetrad[i].velocities    = new float[3 * num_Atoms_N_Evecs[2*i]];
    }
    
    for (i = 0; i < num_Tetrads; i++) {
        
        MPI_Recv(tetrad[i].avg_Structure, 3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+1, comm, &status);
        MPI_Recv(tetrad[i].masses,        3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+2, comm, &status);
        MPI_Recv(tetrad[i].abq,           3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+3, comm, &status);
        MPI_Recv(tetrad[i].eigenvalues,   tetrad[i].num_Evecs,             MPI_FLOAT, 0, TAG_TETRAD+i+4, comm, &status);
        MPI_Recv(&(tetrad[i].eigenvectors[0][0]),  tetrad[i].num_Evecs*3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+5, comm, &status);
        MPI_Recv(tetrad[i].coordinates,   3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+6, comm, &status);
        MPI_Recv(tetrad[i].velocities,    3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+7, comm, &status);
    
        /*
        MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &tetrad[i]);
        MPI_Recv(&tetrad[i], 1, MPI_Tetrad, 0, TAG_TETRAD+i, comm, &status);
        MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);*/
    }
    
    /*
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cout << rank << ": ";
    for (int j = 0; j < num_Tetrads; j++) {
        cout << rank << ": ";
        cout << "index: " << j << endl;
        for (int i = 0; i < tetrad[j].num_Atoms_In_Tetrad; i++) {
            cout << tetrad[j].coordinates[i] << " ";
            if ((i+1)%10 == 0) cout << endl;
        }
        cout << endl << endl;
    }
    cout << endl << endl;*/
    
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_DATA, comm);
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
    
    MPI_Recv(&index, 1, MPI_INT, 0, TAG_ED, comm, &status);
    
    edmd.calculate_ED_Forces(&tetrad[index], ED_Forces[0], edmd.scaled, 3*max_Atoms);
    edmd.calculate_Random_Forces(&tetrad[index], ED_Forces[1]);
    
    ED_Forces[0][3 * max_Atoms + 1] = index;
    
    /*
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cout << "Rank " << rank << ": " << endl;
    for (int i = 0; i < 3 * tetrad[index].num_Atoms_In_Tetrad; i++) {
        cout << setw(11) << ED_Forces[1][i] << " ";
        if ((i+1)%10 == 0) cout << endl;
    }
    cout << endl << endl;*/
    
    cout << "Rank " << setw(3) << rank << " computed ED forces on Tetrad " << setw(3) << index << endl;
    
    MPI_Send(&(ED_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, 0, TAG_ED, comm);
    
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
    
    int i, j, indexes[3];
    
    MPI_Recv(indexes, 3, MPI_INT, 0, TAG_NB, comm, &status);
    
    edmd.calculate_NB_Forces(&tetrad[indexes[1]], &tetrad[indexes[2]], NB_Forces, 3*max_Atoms, 3*max_Atoms);
    
    NB_Forces[0][3*max_Atoms+1] = indexes[0];
    
    /*
    if (indexes[0] == 83) {
        for (int i = 0 ; i < 2; i++) {
            for (int j = 0 ; j < 3*max_Atoms+2; j++) {
                cout << NB_Forces[i][j] << " ";
                if ((i+1)%10 == 0) cout << endl;
            }
        }
        cout << endl;
    }*/
    
    //cout << "Rank " << setw(3) << rank << " computed NB forces on Tetrad " << setw(3) << indexes[0] << " and " << setw(3) << indexes[1] << endl;
    
    MPI_Send(&(NB_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, 0, TAG_NB, comm);
    
}