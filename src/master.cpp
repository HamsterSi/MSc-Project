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
    
    max_Atoms = 0;
    num_Pairs = 0;
    
    comm      = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size); // Get size of MPI processes

}



Master::~Master(void) {
    
    io.array.deallocate_2D_Array(pair_Lists);
    delete [] velocities;
    delete [] coordinates;
    
}



void Master::initialise(void) {
    
    io.read_Cofig(&edmd);

    io.read_Prm();
    
    io.read_Crd();
    
    // Generate the displacement array for DNA base pairs
    io.generate_Displacements();
    
    // Initialise coordinates & velocities of tetrads from the crd file
    io.initialise_Tetrad_Crds();
    
    // Inintialise the frequencies of iterations & outputs
    io.ncycs = io.nsteps/io.ntsync;
    io.ntpr -= io.ntpr % io.ntsync; if (io.ntpr == 0) io.ntpr = 1;
    io.ntwt -= io.ntwt % io.ntsync; if (io.ntwt == 0) io.ntwt = 1;
    
    // Pick out the maximum number of atoms in tetrads
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        max_Atoms = max_Atoms > io.tetrad[i].num_Atoms ? max_Atoms : io.tetrad[i].num_Atoms;
    }
    
    // Allocate memory for pair lists, velocities & coordinates
    num_Pairs = io.prm.num_Tetrads * (io.prm.num_Tetrads - 1) / 2;
    pair_Lists = io.array.allocate_2D_Array(num_Pairs, 2);
    velocities  = new double [3 * io.crd.total_Atoms];
    coordinates = new double [3 * io.crd.total_Atoms];
    
    // Print information of the simulation
    cout << endl << "Initialising simulation..." << endl;
    cout << ">>> MPI Processes: " << size << endl;
    cout << ">>> DNA Shape: ";
    if (io.circular == true) cout << "Circular" << endl;
    else cout << "Linear" << endl;
    cout << ">>> Reading prm & crd file...\nData reading completed.\n" << endl;
    cout << "The number of DNA Base Pairs: " << io.crd.num_BP << endl;
    cout << "The number of DNA Tetrads   : " << io.prm.num_Tetrads << endl;
    cout << "Total number of atoms in DNA: " << io.crd.total_Atoms << endl;
    
}



void Master::send_Parameters(void) {
    
    int i, signal, parameters[2] = {io.prm.num_Tetrads, max_Atoms};
    int * tetrad_Para = new int[2 * io.prm.num_Tetrads];
    double edmd_Para[8] = {edmd.dt, edmd.gamma, edmd.tautp, edmd.temperature,
           edmd.scaled, edmd.mole_Cutoff, edmd.atom_Cutoff, edmd.mole_Least};
    
    // Assign the number of atoms & evecs of tetrads into the sending array
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        tetrad_Para[2 * i] = io.tetrad[i].num_Atoms;
        tetrad_Para[2*i+1] = io.tetrad[i].num_Evecs;
    }
    
    // Send all the parameters to workers
    for (i = 0; i < size - 1; i++) {
        MPI_Send(parameters, 2, MPI_INT   , i + 1, TAG_DATA, comm);
        MPI_Send(edmd_Para , 8, MPI_DOUBLE, i + 1, TAG_DATA, comm);
        MPI_Send(tetrad_Para, 2 * io.prm.num_Tetrads, MPI_INT, i + 1, TAG_DATA, comm);
    }
    
    // Feedback that worker processes have received all parameters
    for (i = 0; i < size - 1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_DATA, comm, &status);
    }
    
    delete [] tetrad_Para;
    
}



void Master::send_Tetrads(void) {
    
    int i, j, signal;
    MPI_Datatype MPI_Tetrad;
    
    // Send tetrads to all worker processes
    for (i = 1; i < size; i++) {
        for (j = 0; j < io.prm.num_Tetrads; j++) {
            
            // Create MPI_Datatype for every Tetrad instance & send them to workers
            MPI_Library::create_MPI_Tetrad(&MPI_Tetrad, &io.tetrad[j]);
            MPI_Send(&io.tetrad[j], 1, MPI_Tetrad, i, TAG_TETRAD+j, comm);
            MPI_Library::free_MPI_Tetrad(&MPI_Tetrad);
            
        }
    }
    
    // Feedback that all worker processes have received tetrads
    for (i = 0; i < size - 1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_TETRAD, comm, &status);
    }
    
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
    double r, ** com = io.array.allocate_2D_Array(io.prm.num_Tetrads, 3);
    
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
    
    io.array.deallocate_2D_Array(com);
    
}



void Master::send_Tetrad_Index(int* i, int* j, int dest, double** buffer) {
    
    if (*i < io.prm.num_Tetrads) { // Send tetrad index for ED calculation
        
        // Assign data to the send buffer for ED calculation
        buffer[0][3 * max_Atoms + 1] = (double) (*i);
        io.array.assignment(3 * io.tetrad[*i].num_Atoms, io.tetrad[*i].coordinates, buffer[0]);
        
        // Send data to available workers
        MPI_Send(&(buffer[0][0]), 3 * max_Atoms + 2, MPI_DOUBLE, dest, TAG_ED, comm);
        
        // Calculate random forces
        edmd.calculate_Random_Forces(&io.tetrad[*i]);
        
    } else if (*j < num_Pairs) { // i >= num_Tetrads, send tetrad indexes for NB calculation
        int index = pair_Lists[*j][0];
        buffer[0][3 * max_Atoms + 1] = (double) index;
        io.array.assignment(3 * io.tetrad[index].num_Atoms, io.tetrad[index].coordinates, buffer[0]);
        
        index = pair_Lists[*j][1];
        buffer[1][3 * max_Atoms + 1] = (double) index;
        io.array.assignment(3 * io.tetrad[index].num_Atoms, io.tetrad[index].coordinates, buffer[1]);
        
        // Send data to available workers
        MPI_Send(&(buffer[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, dest, TAG_NB, comm);
        (*j)++;
    }
    
}



void Master::recv_ED_Forces(double** buffer) {
    
    int i, index = (int) buffer[0][3 * max_Atoms + 1];
    
    io.tetrad[index].ED_Energy = buffer[0][3 * max_Atoms];
    
    for (i = 0; i < 3 * io.tetrad[index].num_Atoms; i++) {
        io.tetrad[index].ED_Forces[i]   = buffer[0][i];
        io.tetrad[index].coordinates[i] = buffer[1][i];
    }
    
}



void Master::recv_NB_Forces(double** buffer, int it) {
    
    int i, index = (int) buffer[it][3 * max_Atoms + 1];
    
    // Assign the NB & EL energies to tetrad
    io.tetrad[index].NB_Energy += buffer[0][3 * max_Atoms];
    io.tetrad[index].EL_Energy += buffer[1][3 * max_Atoms];
    
    // Add the NB forces into tetrad
    for (i = 0; i < 3 * io.tetrad[index].num_Atoms; i++) {
        io.tetrad[index].NB_Forces[i] += buffer[it][i];
    }
    
}



void Master::clip_NB_Forces(void) {
    
    double max_Forces = 1.0;

    for (int i  = 0; i < io.prm.num_Tetrads; i++) {
        for (int j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            io.tetrad[i].NB_Forces[j] = min( max_Forces, io.tetrad[i].NB_Forces[j]);
            io.tetrad[i].NB_Forces[j] = max(-max_Forces, io.tetrad[i].NB_Forces[j]);
        }
    }
    
}



void Master::cal_Forces(void) {
    
    int i, j;
    double ** buffer = io.array.allocate_2D_Array(2, 3 * max_Atoms + 2);

    // Initialise the forces & energies in the tetrads
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            io.tetrad[i].ED_Forces[j] = io.tetrad[i].random_Forces[j] =
            io.tetrad[i].NB_Forces[j] = 0.0;
        }
        io.tetrad[i].ED_Energy = io.tetrad[i].NB_Energy =
        io.tetrad[i].EL_Energy = 0.0;
    }
    
    // Send tetrad indexes & coordinates for ED/NB forces calculation at the beginning
    for (i = 0, j = 0; i < size - 1; i++) {
        send_Tetrad_Index(&i, &j, i + 1, buffer);
    }
    
    // Receive ED or NB forces & energies from workers & send new tetrad indexes & coordinates
    for (; i < num_Pairs + io.prm.num_Tetrads + size - 1 && j <= num_Pairs; i++) {

        MPI_Recv(&(buffer[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status); // Receive ED/NB forces from workers
        
        if (status.MPI_TAG == TAG_ED) {
            recv_ED_Forces(buffer);    // TAG_ED, store ED forces & energies into tetrad
            
        } else if (status.MPI_TAG == TAG_NB) {
            recv_NB_Forces(buffer, 0); // TAG_NB, store NB forces & energies into tetrads
            recv_NB_Forces(buffer, 1);
        }
        
        // Send tetrad index & coordinates if there are some more
        send_Tetrad_Index(&i, &j, status.MPI_SOURCE, buffer);
        
    }
    
    // Clip NB forces into range (-1.0, 1.0)
    clip_NB_Forces();
    
    io.array.deallocate_2D_Array(buffer);
    
}



void Master::cal_Velocities(void) {
    
    // Calculate velocities of every tetrad
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Velocities(&io.tetrad[i]);
    }
    
}



void Master::cal_Coordinate(void) {
    
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
        energies[0] += io.tetrad[i].ED_Energy;
        energies[1] += io.tetrad[i].NB_Energy;
        energies[2] += io.tetrad[i].EL_Energy;
        energies[3] += io.tetrad[i].temperature;
    }
    
    // Calculate the average temperature of tetrads
    energies[3] /= io.prm.num_Tetrads;
    
    // Wrtie out energies
    io.write_Energies(istep, energies);
    
}



void Master::write_Forces(void) {
    
    int i, j, index;
    double * ED_Forces     = new double [3 * io.crd.total_Atoms];
    double * random_Forces = new double [3 * io.crd.total_Atoms];
    double * NB_Forces     = new double [3 * io.crd.total_Atoms];
    
    // Initialise arrays
    for (i = 0; i < 3 * io.crd.total_Atoms; i++) {
        ED_Forces[i] = random_Forces[i] = NB_Forces[i] = 0.0;
    }
    
    // Gather all forcees together into three arrays
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (index = io.displs[i], j = 0; j < 3 * io.tetrad[i].num_Atoms; index++, j++) {
            ED_Forces[index]     += io.tetrad[i].ED_Forces[j];
            random_Forces[index] += io.tetrad[i].random_Forces[j];
            NB_Forces[index]     += io.tetrad[i].NB_Forces[j];
        }
    }
    
    // Write out all forces
    io.write_Forces(ED_Forces, random_Forces, NB_Forces);
    
    delete []ED_Forces;
    delete []random_Forces;
    delete []NB_Forces;
    
}



void Master::write_Trajectories(int istep) {
    
    int index = io.displs[io.crd.num_BP - 3];

    io.write_Trajectory(istep, index, coordinates);
    
}



void Master::write_Crds(void) {
    
    io.update_Crd(velocities, coordinates);
    
}



void Master::finalise(void) {
    
    cout << "Writing Energies out at     " << io.energy_File << endl;
    cout << "Writing Forces   out at     " << io.forces_File << endl;
    cout << "Writing Trajectories out at " << io.trj_File << endl;
    cout << "Writing New Crd  out at     " << io.new_Crd_File << endl << endl;

    // Send signal to stop all worker processes
    for (int signal = -1, i = 1; i < size; i++) {
        MPI_Send(&signal, 1, MPI_INT, i, TAG_DEATH, comm);
    }
    
    cout << "Simulation ended.\n" << endl;
    
}

