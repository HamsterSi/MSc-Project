
#include "worker.hpp"

/*
 * Function:  Receive tetrad parameters, Compute ED forces of tetrad &
 *            Send the calculated ED forces back to the master.
 *
 * Parameter: None
 *
 * Return:    None
 */
int Worker_Management::worker_ED_Forces(void) {
    
    Tetrad tetrad;
    EDMD edmd;
    float *ED_Forces, scaled, ED_Energy;
    
    MPI_Datatype MPI_Tetrad;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status *status;
    
    MPI_Library::create_MPI_Tetrad(MPI_Tetrad, tetrad.num_Atoms_In_Tetrad, tetrad.num_Evecs);
    MPI_Recv(&tetrad, 1, MPI_Tetrad, 0, TAG_ED, comm, status);
    
    ED_Forces = new float[3 * tetrad.num_Atoms_In_Tetrad];
    edmd.calculate_ED_Forces(tetrad, ED_Forces, edmd.scaled, &ED_Energy);
    
    MPI_Send(ED_Forces, 3 * tetrad.num_Atoms_In_Tetrad, MPI_FLOAT, 0, TAG_ED, comm);
    MPI_Send(&ED_Energy, 1, MPI_FLOAT, 0, TAG_ED, comm);
    
    delete []ED_Forces;
    MPI_Library::free_MPI_Tetrad(MPI_Tetrad);
    
    return 1;
}


/*
 * Function:  Receive tetrad parameters, Compute NB forces of tetrad &
 *            Send the calculated NB forces back to the master.
 *
 * Parameter: None
 *
 * Return:    None
 */
int Worker_Management::worker_NB_Forces(void) {
    
    Tetrad tetrad[2];
    EDMD edmd;
    int source, dest;
    float **NB_Forces, Energies[2]; //NB_Energy, Electrostatic_Energy;
    
    MPI_Datatype MPI_Tetrad;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status *status;
    
    MPI_Recv(tetrad, 2, MPI_Tetrad, 0, TAG_NB, comm, status);
    
    NB_Forces = new float*[2];
    int num_Atoms = max(tetrad[0].num_Atoms_In_Tetrad, tetrad[1].num_Atoms_In_Tetrad);
    NB_Forces[0] = new float [3 * num_Atoms];
    NB_Forces[1] = new float [3 * num_Atoms];
    
    //edmd.calculate_NB_Forces(tetrad, NB_Forces, &NB_Energy, &Electrostatic_Energy);
    edmd.calculate_NB_Forces(tetrad, NB_Forces, &Energies[0], &Energies[1]);
    
    MPI_Send(NB_Forces, 3 * num_Atoms, MPI_FLOAT, 0, TAG_NB, comm);
    MPI_Send(Energies, 2, MPI_FLOAT, 0, TAG_NB, comm);
    
    delete[] NB_Forces[0];
    delete[] NB_Forces[1];
    delete[] NB_Forces;
    
    return 1;
}