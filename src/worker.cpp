
#include "worker.hpp"


/*
 * Function:  The constructor of Worker class.
 *
 * Parameter: None
 *
 * Return:    None
 */
Worker::Worker(void) {
    
    comm    =   MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank); // Get the rank of worker process
    
}




/*
 * Function:  The destructor of Worker class.  Deallocate memory.
 *
 * Parameter: None
 *
 * Return:    None
 */
Worker::~Worker(void) {
    
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
 * Function:  Workers receive the number of tetrads, number of atoms in every tetrad
 *            and number of evecs from master
 *
 * Parameter: None
 *
 * Return:    None
 */
void Worker::recv_Parameters(void) {
    
    int i, j, signal = 1, parameters[3];
    
    // Receive parameters from the master process
    MPI_Recv(parameters, 3, MPI_INT, 0, TAG_DATA, comm, &status);
    num_Tetrads = parameters[0]; // The number of tetrads
    max_Atoms   = parameters[1]; // The maximum number of atoms in tetrads
    iteration   = parameters[2];
    
    if (iteration == 0) {
        
        int * num_Atoms_n_Evecs = new int[2 * num_Tetrads];;
        
        // Receive the number of atoms & evecs in tetrads
        MPI_Recv(num_Atoms_n_Evecs, 2 * num_Tetrads, MPI_INT, 0, TAG_DATA, comm, &status);
        
        // Allocate memory for all tetrads
        tetrad = new Tetrad[num_Tetrads];
        
        for (i = 0; i < num_Tetrads; i++) {
            tetrad[i].num_Atoms = num_Atoms_n_Evecs[2*i];
            tetrad[i].num_Evecs = num_Atoms_n_Evecs[2*i+1];
            
            tetrad[i].avg_Structure = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].masses        = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].abq           = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].eigenvalues   = new float [tetrad[i].num_Evecs];
            tetrad[i].eigenvectors  = new float*[tetrad[i].num_Evecs];
            for (j = 0; j < tetrad[i].num_Evecs; j++) {
                tetrad[i].eigenvectors[j] = new float[3 * tetrad[i].num_Atoms];
            }
            tetrad[i].coordinates   = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].velocities    = new float[3 * tetrad[i].num_Atoms];
            
            tetrad[i].ED_Forces     = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].random_Forces = new float[3 * tetrad[i].num_Atoms];
            tetrad[i].NB_Forces     = new float[3 * tetrad[i].num_Atoms];
            
            tetrad[i].energies[0] = tetrad[i].energies[1] = 0.0;
            tetrad[i].energies[2] = tetrad[i].temperature = 0.0;
        }
        
        delete []num_Atoms_n_Evecs;
    }
    
    // Send feedback to master that has received all parameters.
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_DATA, comm);
    
}




/*
 * Function:  Workers receive tetrads from master
 *
 * Parameter: None
 *
 * Return:    None
 */
void Worker::recv_Tetrads(void) {
    
    int i, signal = 1;
    MPI_Datatype MPI_Tetrad;
    
    // Receive tetrads from master
    for (i = 0; i < num_Tetrads; i++) {
        
        // If it is the 1st iteration, then it should receive all parameters of tetrads
        if (iteration == 0) {
            MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &tetrad[i]);
            MPI_Recv(&tetrad[i], 1, MPI_Tetrad, 0, TAG_TETRAD+rank+i, comm, &status);
            MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);
        
        // 2nd, 3rd, 4th... iterations only needs to receive velocities & coordinates
        } else {
            MPI_Recv(tetrad[i].velocities,  3 * tetrad[i].num_Atoms, MPI_FLOAT, 0, TAG_TETRAD+rank+i+1, comm, &status);
            MPI_Recv(tetrad[i].coordinates, 3 * tetrad[i].num_Atoms, MPI_FLOAT, 0, TAG_TETRAD+rank+i+2, comm, &status);
        }
    }
    
    // Send feedback to master that has received all tetrads
    MPI_Send(&signal, 1, MPI_INT, 0, TAG_TETRAD, comm);
    
}




/*
 * Function:  Compute ED forces of tetrads.
 *            Send calculated ED/random forces back to master.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Worker::ED_Calculation(void) {
    
    int i, index;
    float ED_Forces[2][3 * max_Atoms + 2];
    
    // Receive the tetrad index from the master process
    MPI_Recv(&index, 1, MPI_INT, 0, TAG_ED, comm, &status);
    
    // Calculate ED forces (ED energy) & random forces
    edmd.calculate_ED_Forces(&tetrad[index], edmd.scaled);
    edmd.calculate_Random_Forces(&tetrad[index], rank);
    
    // Assign ED forces & random Forces to the 2D array for sending once
    for (i = 0; i < 3 * tetrad[index].num_Atoms; i++) {
        ED_Forces[0][i] = tetrad[index].ED_Forces[i];
        ED_Forces[1][i] = tetrad[index].random_Forces[i];
    }
    
    // Need to send the ED Energy back
    ED_Forces[0][3 * max_Atoms] = tetrad[index].energies[0];
    
    // Need to send tetrad index back
    ED_Forces[0][3 * max_Atoms + 1] = index;
    
    // Send the calculated ED forces, ED energy, random forces & index back
    MPI_Send(&(ED_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, 0, TAG_ED, comm);
    
}




/*
 * Function:  Compute NB forces of tetrads.
 *            Send calculated NB forces & energies back to master.
 *
 * Parameter: None
 *
 * Return:    None
 */
void Worker::NB_Calculation(void) {
    
    int i, indexes[2];
    float NB_Forces[2][3 * max_Atoms + 2];
    
    // Receive tetrad indexes for NB forces calculation
    MPI_Recv(indexes, 2, MPI_INT, 0, TAG_NB, comm, &status);
    
    // Calculate NB forces, NB energy & Electrostatic Energy
    edmd.calculate_NB_Forces(&tetrad[indexes[0]], &tetrad[indexes[1]]);
    
    // Assign NB forces to the 2D array to send back once
    for (i = 0; i < 3 * max_Atoms; i++) {
        NB_Forces[0][i] = tetrad[indexes[0]].NB_Forces[i];
        NB_Forces[1][i] = tetrad[indexes[1]].NB_Forces[i];
    }
    
    // Need to send NB Energy & Electrostatic Energy back (bacause the energies in both tetrads are the same when calculation, so only need to send one set)
    NB_Forces[0][3 * max_Atoms] = tetrad[indexes[0]].energies[1];
    NB_Forces[1][3 * max_Atoms] = tetrad[indexes[0]].energies[2];
    
    // Need to send both tetrad indexes back
    NB_Forces[0][3 * max_Atoms + 1] = indexes[0];
    NB_Forces[1][3 * max_Atoms + 1] = indexes[1];
    
    // Send NB forces, energies & indexes back to master
    MPI_Send(&(NB_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_FLOAT, 0, TAG_NB, comm);
    
}





/*
 * Function:  Receive the terminate signal from master & terminate work
 *
 * Parameter: None
 *
 * Return:    None
 */
int Worker::terminate(void) {
    
    int signal;
    
    // Receive the terminate signal from master
    MPI_Recv(&signal, 1, MPI_INT, 0, TAG_DEATH, MPI_COMM_WORLD, &status);

    return signal;
}



