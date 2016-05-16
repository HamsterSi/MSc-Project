
#include "worker.hpp"

Worker_Management::Worker_Management(void) {
    comm = MPI_COMM_WORLD;
}

/*
 *
 */
void Worker_Management::data_Receiving(void) {
    
    MPI_Recv(data, 3, MPI_INT, 0, TAG_DATA, comm, &status);
}

/*
 *
 */
void Worker_Management::tetrad_Receiving(void) {

    MPI_Datatype MPI_Tetrad;
    
    tetrad = new Tetrad[data[0]];
    
    MPI_Library::create_MPI_Tetrad(MPI_Tetrad, data[1], data[2]);

    MPI_Recv(tetrad, data[0], MPI_Tetrad, 0, TAG_TETRAD, comm, &status);
    
    MPI_Library::free_MPI_Tetrad(MPI_Tetrad);
    
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
    
    float* ED_Forces = new float[3 * data[1] + 1];
    
    edmd.calculate_ED_Forces(tetrad[index], ED_Forces, edmd.scaled, 3*data[1]);
    
    MPI_Send(ED_Forces, 3 * data[1] + 1, MPI_FLOAT, 0, TAG_ED, comm);
    
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
    NB_Forces[0] = new float [3 * data[1] + 2];
    NB_Forces[1] = new float [3 * data[1] + 2];
    
    edmd.calculate_NB_Forces(te, NB_Forces, 3*data[1], 3*data[1]);
    NB_Forces[0][3*data[1]+1] = index[0];
    NB_Forces[1][3*data[1]+1] = index[1];
     
    MPI_Send(NB_Forces, 2 * 3 * data[1] + 4, MPI_FLOAT, 0, TAG_NB, comm);
    
    delete[] NB_Forces[0];
    delete[] NB_Forces[1];
    delete[] NB_Forces;

}