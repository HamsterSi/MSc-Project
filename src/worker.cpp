
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
    delete []num_Atoms_N_Evecs;
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
    int signal = -1;
    
    MPI_Recv(&parameters, 2, MPI_INT, 0, TAG_DATA, comm, &status);
    
    num_Atoms_N_Evecs = new int[2 * parameters[0]];
    
    MPI_Recv(num_Atoms_N_Evecs, 2 * parameters[0], MPI_INT, 0, TAG_DATA, comm, &status);
    
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
    int signal = -1;
    MPI_Datatype MPI_Tetrad;
    
    tetrad = new Tetrad[parameters[0]];
    for (int i = 0; i < parameters[0]; i++) {  // Allocate memory
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
    
    for (int i = 0; i < parameters[0]; i++) {
        MPI_Recv(tetrad[i].avg_Structure, 3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+1, comm, &status);
        
        MPI_Recv(tetrad[i].masses,        3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+2, comm, &status);
        
        MPI_Recv(tetrad[i].abq,           3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+3, comm, &status);
        
        MPI_Recv(tetrad[i].eigenvalues,   tetrad[i].num_Evecs,             MPI_FLOAT, 0, TAG_TETRAD+i+4, comm, &status);
        
        MPI_Recv(&(tetrad[i].eigenvectors[0][0]),  tetrad[i].num_Evecs*3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+5, comm, &status);
        
        MPI_Recv(tetrad[i].coordinates,   3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+6, comm, &status);
        
        MPI_Recv(tetrad[i].velocities,    3*tetrad[i].num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_TETRAD+i+7, comm, &status);
        
        //MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &tetrad[i]);
        //MPI_Recv(&tetrad[i], 1, MPI_Tetrad, 0, TAG_TETRAD+i, comm, &status);
        //MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);
    }
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
    
    float* ED_Forces = new float[3 * parameters[1] + 1];
    
    edmd.calculate_ED_Forces(tetrad[index], ED_Forces, edmd.scaled, 3*parameters[1]);
    
    MPI_Send(ED_Forces, 3 * parameters[1] + 1, MPI_FLOAT, 0, TAG_ED, comm);
    
    delete []ED_Forces;
    
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
    
    int index[2];
    
    MPI_Recv(index, 2, MPI_INT, 0, TAG_NB, comm, &status);
    
    Tetrad te[2] = {tetrad[index[0]], tetrad[index[1]]};
    
    float** NB_Forces = new float*[2];
    NB_Forces[0] = new float [3 * parameters[1] + 2];
    NB_Forces[1] = new float [3 * parameters[1] + 2];
    
    edmd.calculate_NB_Forces(te, NB_Forces, 3*parameters[1], 3*parameters[1]);
    NB_Forces[0][3*parameters[1]+1] = index[0];
    NB_Forces[1][3*parameters[1]+1] = index[1];
     
    MPI_Send(NB_Forces, 2 * 3 * parameters[1] + 4, MPI_FLOAT, 0, TAG_NB, comm);
    
    delete[] NB_Forces[0];
    delete[] NB_Forces[1];
    delete[] NB_Forces;

}