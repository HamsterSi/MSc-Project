/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                 MSc in High Performance Computing, EPCC                      *
 *                      The University of Edinburgh                             *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  master.cpp
 * Brief: Implementation of the Master class functions 
 */

#include "master.hpp"


Master::Master(void) {
    
    num_Pairs = 0;
    
    comm      = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size); // Get size of MPI processes

}



Master::~Master(void) {
    
    // Deallocate memory of arrays
    array.deallocate_2D_Array(pair_Lists);
    delete [] velocities;
    delete [] coordinates;
    
    // Free the MPI_Datatyep
    MPI_Library::free_MPI_Force(&MPI_Force);
    
}



void Master::initialise(void) {
    
    io.read_Cofig(&edmd);
    io.read_Prm();
    io.read_Crd();
    
    // Initialise coordinates & velocities of tetrads from the crd file
    io.initialise_Tetrad_Crds();
    
    // Inintialise the frequencies of iterations & outputs
    io.ntwt -= io.ntwt % io.ntsync; if (io.ntwt == 0) io.ntwt = 1;
    io.ntpr -= io.ntpr % io.ntsync; if (io.ntpr == 0) io.ntpr = 1;
    
    // The total number of NB pairs & allocate memory for pair lists
    num_Pairs   = io.prm.num_Tetrads * (io.prm.num_Tetrads - 1) / 2;
    pair_Lists  = array.allocate_2D_Array(num_Pairs, 2);
    
    // Create new MPI_Datatype to receive ED/NB forces of all tetrads from rank 1
    MPI_Library::create_MPI_Force(&MPI_Force, io.prm.num_Tetrads, io.tetrad);
    
    // Allocate memory for velocities & coordinates
    velocities  = new double [3 * io.crd.total_Atoms];
    coordinates = new double [3 * io.crd.total_Atoms];
    
    // Print information of the simulation
    cout << endl << "Initialising simulation..." << endl << endl;
    cout << "The number of MPI Processes : " << size << endl << endl;
    
    cout << "Siulation parameters:" << endl;
    cout << ">>> Total number of iterations  : " << io.nsteps << endl;
    cout << ">>> Frequency of synchronization: " << io.ntsync << endl;
    cout << ">>> Frequency of writing energy & temperature, trajectory: " << io.ntwt << endl;
    cout << ">>> Frequency of updating th coordinates file: " << io.ntpr << endl << endl;
    
    cout << "DNA Information:" << endl;
    cout << ">>> The number of DNA Base Pairs  : " << io.crd.num_BP      << endl;
    cout << ">>> The number of DNA Tetrads     : " << io.prm.num_Tetrads << endl;
    cout << ">>> Total number of atoms in DNA  : " << io.crd.total_Atoms << endl << endl;
    
    cout << "Inputs:" << endl;
    cout << ">>> The tetrad parameter file path: " << io.prm_File << endl;
    cout << ">>> The coordinate file path      : " << io.crd_File << endl << endl;
    
    cout << "Outputs:" << endl;
    cout << ">>> Energy & temperature file path: " << io.energy_File  << endl;
    cout << ">>> The trajectory file path      : " << io.trj_File     << endl;
    cout << ">>> The new coordinates file path : " << io.new_Crd_File << endl << endl;
    
}



void Master::send_Parameters(void) {
    
    int i, * tetrad_Para = new int[2 * io.prm.num_Tetrads];
    double edmd_Para[9] = { edmd.dt, edmd.gamma, edmd.tautp,
        edmd.temperature, edmd.scaled, edmd.mole_Cutoff,
        edmd.atom_Cutoff, edmd.mole_Least, io.prm.num_Tetrads };
    MPI_Status status;
    
    // Assign the number of atoms & evecs of tetrads into the sending array
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        tetrad_Para[2 * i] = io.tetrad[i].num_Atoms;
        tetrad_Para[2*i+1] = io.tetrad[i].num_Evecs;
    }
    
    // Broadcast the simulation parameters and tetrad parameters to all workers
    cout << "Master sending parameters to workers..." << endl;
    MPI_Bcast(edmd_Para  , 9, MPI_DOUBLE, 0, comm);
    MPI_Bcast(tetrad_Para, 2 * io.prm.num_Tetrads, MPI_INT, 0, comm);
    
    delete [] tetrad_Para;
    
}



void Master::send_Tetrads(void) {
    
    MPI_Datatype MPI_Tetrad;
    
    cout << "Master sending tetrads to workers...\n" << endl;
    
    // Create MPI_Datatype "MPI_Tetrad" for all tetrads
    MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, io.prm.num_Tetrads, io.tetrad);
    
    // Broadcast tetrads to all workers
    MPI_Bcast(io.tetrad, 1, MPI_Tetrad, 0, comm);
    
    // Free the MPI_Datatype "MPI_Tetrad"
    MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);
    
    cout << "Simulation starting..." << endl;
}



void Master::cal_Centre_of_Mass(double** com) {
    
    // Calculate the centre of mass (actually, centre of geom)
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        
        com[i][0] = com[i][1] = com[i][2] = 0.0;
        
        for (int j = 0; j < 3 * io.tetrad[i].num_Atoms; j += 3) {
            com[i][0] += io.tetrad[i].coordinates[ j ];
            com[i][1] += io.tetrad[i].coordinates[j+1];
            com[i][2] += io.tetrad[i].coordinates[j+2];
        }
        
        com[i][0] /= io.tetrad[i].num_Atoms;
        com[i][1] /= io.tetrad[i].num_Atoms;
        com[i][2] /= io.tetrad[i].num_Atoms;
        
    }
    
}



void Master::generate_Pair_Lists(void) {
    
    int i, j;
    double r, ** com = array.allocate_2D_Array(io.prm.num_Tetrads, 3);
    
    cal_Centre_of_Mass(com); // Calculate the centre of mass
    
    // Loop to generate pairlists
    for (num_Pairs = 0, i = 0; i < io.prm.num_Tetrads; i++) {
        for (j = i + 1; j < io.prm.num_Tetrads; j++) {
            
            r = (com[i][0] - com[j][0]) * (com[i][0] - com[j][0]) + (com[i][1] - com[j][1]) * (com[i][1] - com[j][1]) + (com[i][2] - com[j][2]) * (com[i][2] - com[j][2]);
            
            // If r exceeds mole_Cutoff then no interaction between these two mols
            if ((r < (edmd.mole_Cutoff * edmd.mole_Cutoff)) && (abs(i - j) > edmd.mole_Least) &&
                (abs(i - j) < (io.prm.num_Tetrads - edmd.mole_Least))) {
                
                if ((abs(i - j) % 2 == 1 && (i - j) < 0) || (abs(i - j) % 2 == 0 && (i - j) > 0)) {
                    pair_Lists[num_Pairs][0] = j; pair_Lists[num_Pairs][1] = i;
                } else {
                    pair_Lists[num_Pairs][0] = i; pair_Lists[num_Pairs][1] = j;
                }
            
                num_Pairs++;
            }
        }
    }
    
    array.deallocate_2D_Array(com);
    
}



void Master::initialise_Forces_n_Energies(void) {

    // Set all elements of the force arrays to zero
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        
        for (int j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            io.tetrad[i].ED_Forces[j]     = 0.0;
            io.tetrad[i].random_Forces[j] = 0.0;
            io.tetrad[i].NB_Forces[j]     = 0.0;
        }
        
        io.tetrad[i].ED_Forces[3 * io.tetrad[i].num_Atoms] =
        io.tetrad[i].NB_Forces[3 * io.tetrad[i].num_Atoms] =
        io.tetrad[i].NB_Forces[3 * io.tetrad[i].num_Atoms + 1] = 0.0;
    }
    
}




void Master::send_n_Recv(int* i, int* j, int dest, int index[], MPI_Request* send_Request, MPI_Request* recv_Request) {
    
    // Send new tetrad index for ED calculation & receive the calculation finish signal,
    // Calculate the random forces
    if (*i < io.prm.num_Tetrads) {

        index[0] = * i;
        index[1] = 0.0;
        index[2] = TAG_ED;
        
        MPI_Isend(index, 3, MPI_INT, dest, TAG_INDEX, comm, send_Request);
        
        MPI_Irecv(index, 3, MPI_INT, dest, TAG_INDEX, comm, recv_Request);
        
    // i >= number of tetrads, then send two tetrad indexes for NB calculation,
    // and receive the calculation finish signal
    } else if (*j < num_Pairs) {
        
        index[0] = pair_Lists[*j][0];
        index[1] = pair_Lists[*j][1];
        index[2] = TAG_NB;
        
        MPI_Isend(index, 3, MPI_INT, dest, TAG_INDEX, comm, send_Request);
        
        MPI_Irecv(index, 3, MPI_INT, dest, TAG_INDEX, comm, recv_Request);
        
        (*j)++;
    }
    
}



void Master::calculate_Forces(void) {

    int i, j, rank, index[size-2][3];
    MPI_Request send_Request[size-2], recv_Request[size-2];
    MPI_Status send_Status, recv_Status;

    // Initialise the ED/NB forces & energies of tetrads to 0,
    // and send signal to rank 1 to initialise forces as well.
    initialise_Forces_n_Energies();
    index[0][0] = index[0][1] = index[0][2] = TAG_CLEAN;
    MPI_Send(index[0], 3, MPI_INT, 1, TAG_INDEX, comm);
    
    // Send new tasks to worker from rank 2 ~ & receive their finish singal.
    for (i = 0, j = 0; i < size - 2; i++) {
        send_n_Recv(&i, &j, i + 2, index[i], &send_Request[i], &recv_Request[i]);
    }
    
    // Loop to receive finish singal from workers & generate new tasks and send to idle workers
    for (; i < num_Pairs + io.prm.num_Tetrads + size - 2 && j <= num_Pairs; i++) {
        
        // Wait any worker to finish its calculation
        MPI_Waitany(size - 2, recv_Request, &rank, &recv_Status);
        MPI_Wait(&(send_Request[rank]), &send_Status);
        
        // Send new tetrad index(es) to idle workers
        send_n_Recv(&i, &j, rank + 2, index[rank], &send_Request[rank], &recv_Request[rank]);
        
    }
    
    // Send signal to require all ED/NB forces & receive all forces from rank 1
    index[0][0] = index[0][1] = index[0][2] = TAG_FORCE;
    MPI_Send(index[0], 3, MPI_INT, 1, TAG_INDEX, comm);
    MPI_Recv(io.tetrad, 1, MPI_Force, 1, TAG_FORCE, comm, &recv_Status);
    
    double temp[6][3];
    for (int i = 0; i < 6; i++) {  temp[i][0] = temp[i][1] = temp[i][2] = 0.0;  }
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        for (int j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            temp[5][0] += io.tetrad[i].ED_Forces[j];
            temp[5][1] += io.tetrad[i].NB_Forces[j];
            temp[5][2] += io.tetrad[i].coordinates[j];
        }
    }
    //cout << "M new forces: " << temp[5][0] << " " << temp[5][1] << " " << temp[5][2] << endl;

}




void Master::update_Velocities(void) {
    
    // Calculate velocities of every tetrad
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Velocities(&io.tetrad[i]);
    }
    
}



void Master::update_Coordinates(void) {
    
    // Calculate coordinates of all tetrads
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Coordinates(&io.tetrad[i]);
    }

}



void Master::data_Processing(void) {
    
    int i, j, index;
    
    // Initialise velocities & coordinates array
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        velocities[i] = coordinates[i] = 0.0;
    }
    
    // Gather all velocities & coordinates into a single array
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = io.displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; index++, j++) {
            velocities [index] += io.tetrad[i].velocities [j];
            coordinates[index] += io.tetrad[i].coordinates[j];
        }
    }
    
    // Process the beginning & end tetrads
    index = io.displs[io.crd.num_BP - 3];
    for (i = 0; index < io.displs[io.crd.num_BP]; index++, i++) {
        velocities [index] += velocities [i]; velocities [i] = velocities [index];
        coordinates[index] += coordinates[i]; coordinates[i] = coordinates[index];
    }

    // Divide velocities & coordinates by 4
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        velocities[i] *= 0.25; coordinates[i] *= 0.25;
    }

    // Restore the velocities & coordinates back to tetrads
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = io.displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; index++, j++) {
            io.tetrad[i].velocities [j] = velocities [index];
            io.tetrad[i].coordinates[j] = coordinates[index];
        }
    }
    
}



void Master::write_Energy(int istep) {
    
    double energies[4] = {0.0};
    
    // Gather energies & temperature of tetrads together
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        energies[0] += io.tetrad[i].ED_Forces[3 * io.tetrad[i].num_Atoms];
        energies[1] += io.tetrad[i].NB_Forces[3 * io.tetrad[i].num_Atoms];
        energies[2] += io.tetrad[i].NB_Forces[3 * io.tetrad[i].num_Atoms + 1];
        energies[3] += io.tetrad[i].temperature;
    }
    
    // Calculate the average temperature of tetrads
    energies[3] /= io.prm.num_Tetrads;
    
    // Wrtie out energies
    io.write_Energies(istep, energies);
    
}



void Master::write_Trajectories(int istep) {

    // io.displs[io.crd.num_BP - 3] is the index that the afterwards 3
    // tetrads the same as the first three tetrads
    io.write_Trajectory(istep, io.displs[io.crd.num_BP - 3], coordinates);
    
}



void Master::write_Info(int istep) {
    
    write_Energy(istep + io.ntsync);
    
    write_Trajectories(istep);
    
}



void Master::write_Crds(void) {
    
    io.update_Crd_File(velocities, coordinates);
    
}



void Master::finalise(void) {

    int i, index[3];
    
    // Send signal to terminate all workers
    for (index[2] = TAG_SIGNAL, i = 1; i < size; i++) {
        MPI_Send(index, 3, MPI_INT, i, TAG_INDEX, comm);
    }
    
    cout << "\nSimulation ecompleted..." << endl << endl;
    
}

