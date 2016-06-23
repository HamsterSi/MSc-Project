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
    effective_Pairs = 0;
    
    comm      = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size); // Get size of MPI processes

}



Master::~Master(void) {
    
    delete [] pair_Lists;
    delete [] velocities;
    delete [] coordinates;
    
}



void Master::initialise(void) {
    
    io.read_Cofig(&edmd);

    io.read_Prm();
    
    io.read_Crd();
    
    // Initialise the displacement array & generate diplacements of BP
    io.generate_Displacements();
    
    // Initialise coordinates (velocities) of tetrads from crd file
    io.initialise_Tetrad_Crds();
    
    // Inintialise the frequencies of iterations and writing out files
    io.ncycs = io.nsteps/io.ntsync;
    io.ntpr -= io.ntpr % io.ntsync; if (io.ntpr == 0) io.ntpr = 1;
    io.ntwt -= io.ntwt % io.ntsync; if (io.ntwt == 0) io.ntwt = 1;
    
    // Pick out the maximum number of atoms in tetrads
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        max_Atoms = max_Atoms > io.tetrad[i].num_Atoms ? max_Atoms : io.tetrad[i].num_Atoms;
    }
    
    // Allocate memory for storing velocities & coordinates of the whole DNA
    velocities  = new double [3 * io.crd.total_Atoms];
    coordinates = new double [3 * io.crd.total_Atoms];
    
    cout << endl << "Simulation starting..." << endl;
    cout << ">>> MPI Processes: " << size << endl;
    cout << ">>> DNA Shape: ";
    if (io.circular == true) cout << "Circular" << endl;
    else cout << "Linear" << endl;
    
    cout << "Reading prm & crd file...\nData reading completed." << endl << endl;
    cout << "The number of DNA Base Pairs: " << io.crd.num_BP << endl;
    cout << "The number of DNA Tetrads   : " << io.prm.num_Tetrads << endl;
    cout << "Total number of atoms in DNA: " << 3 * io.crd.total_Atoms << endl;
    
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

    delete [] tetrad_Para;
    
    // Feedback that worker processes have received all parameters
    for (i = 0; i < size - 1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_DATA, comm, &status);
    }
    
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



void Master::send_Vels_n_Crds(void) {
    
    int i, j, signal;
    
    // Send velocities & coordinates of all tetrad to all worker processes
    for (i = 1; i < size; i++) {
        for (j = 0; j < io.prm.num_Tetrads; j++) {
            MPI_Send(io.tetrad[j].velocities,  3 * io.tetrad[j].num_Atoms, MPI_DOUBLE, i, TAG_CRDS+j+1, comm);
            MPI_Send(io.tetrad[j].coordinates, 3 * io.tetrad[j].num_Atoms, MPI_DOUBLE, i, TAG_CRDS+j+2, comm);
        }
    }
    
    // Feedback that all worker processes have received all velocities & coordinates
    for (i = 0; i < size - 1; i++) {
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_CRDS, comm, &status);
    }
    
}



void Master::generate_Pair_Lists(void) {
    
    num_Pairs = io.prm.num_Tetrads * (io.prm.num_Tetrads - 1) / 2;
    pair_Lists = new int [2 * num_Pairs];
    
    int i, j, k, pairs;
    double r, ** com = new double * [io.prm.num_Tetrads];
    for (i = 0; i < io.prm.num_Tetrads; i++) { com[i] = new double[3]; }
    
    // The centre of mass (actually, centre of geom)
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        
        com[i][0] = com[i][1] = com[i][2] = 0.0;
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms; ) {
            com[i][0] += io.tetrad[i].coordinates[ j ];
            com[i][1] += io.tetrad[i].coordinates[j+1];
            com[i][2] += io.tetrad[i].coordinates[j+2];
            j += 3;
        }
        com[i][0] /= io.tetrad[i].num_Atoms;
        com[i][1] /= io.tetrad[i].num_Atoms;
        com[i][2] /= io.tetrad[i].num_Atoms;
    }
    
    // Loop to generate pairlists
    for (pairs = 0, effective_Pairs = 0, i = 0; i < io.prm.num_Tetrads; i++) {
        for (j = i + 1; j < io.prm.num_Tetrads; pairs++, j++) {
            
            // If r exceeds mole_Cutoff then no interaction between these two mols
            for (r = 0.0, k = 0; k < 3; k++) {
                r += (com[i][k] - com[j][k]) * (com[i][k] - com[j][k]);
            }
            
            if ((r < (edmd.mole_Cutoff * edmd.mole_Cutoff)) && (abs(i - j) > edmd.mole_Least) &&
                (abs(i - j) < (io.prm.num_Tetrads - edmd.mole_Least))) {
                
                pair_Lists[2 * pairs] = i; pair_Lists[2*pairs+1] = j;
                effective_Pairs++;
                
            } else { pair_Lists[2 * pairs] = -1; pair_Lists[2*pairs+1] = -1; }
            
        }
    }
    
    for (i = 0; i < io.prm.num_Tetrads; i++) { delete [] com[i]; }
    delete [] com;
    
}



void Master::send_Tetrad_Index(int* i, int* j, int source) {
    
    // Send tetrad index for ED calculation & Calculate random forces
    if ((*i) < io.prm.num_Tetrads) {
        MPI_Send(i, 1, MPI_INT, source, TAG_ED, comm);
        edmd.calculate_Random_Forces(&io.tetrad[(*i)]);
        
    // i >= num_Tetrads, send tetrad indexes for NB calculation
    } else {
        for (; (*j) < num_Pairs; (*j)++) {
            if (pair_Lists[2 * (*j)] != -1) {
                int indexes[2] = { pair_Lists[2 * (*j)],  pair_Lists[2 * (*j) + 1] };
                MPI_Send(indexes, 2, MPI_INT, source, TAG_NB, comm); (*j)++; break;
            }
        }
    }
    
}



void Master::cal_Forces(void) {
    
    int i, j, k, flag, index;
    double max_Forces = 1.0, temp_Forces[2][3 * max_Atoms + 2];
    
    // Initialise the forces & energies in the tetrads
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            io.tetrad[i].ED_Forces[j]     = 0.0;
            io.tetrad[i].random_Forces[j] = 0.0;
            io.tetrad[i].NB_Forces[j]     = 0.0;
        }
        io.tetrad[i].ED_Energy = 0.0;
        io.tetrad[i].NB_Energy = 0.0;
        io.tetrad[i].EL_Energy = 0.0;
    }
    
    // Send tetrad indexes for ED/NB forces calculation at the beginning
    for (i = 0, j = 0; i < size - 1; i++) { send_Tetrad_Index(&i, &j, i + 1); }
    
    // When there are still forces waiting for calculating,
    // receive ED/NB forces from workers & send new tetrad indexes to workers
    for (; i < effective_Pairs + io.prm.num_Tetrads + size - 1 && j <= num_Pairs; i++) {

        // Receive ED/NB forces from workers
        MPI_Recv(&(temp_Forces[0][0]), 2 * (3 * max_Atoms + 2), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
        
        // If the MPI Tag indicates it's ED forces then store ED forces & random forces
        if (status.MPI_TAG == TAG_ED) {
            index = (int) temp_Forces[0][3 * max_Atoms + 1];
            io.tetrad[index].ED_Energy = temp_Forces[0][3 * max_Atoms];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].ED_Forces[k]   = temp_Forces[0][k];
                io.tetrad[index].coordinates[k] = temp_Forces[1][k];
            }
            
        // The MPI Tag shows it's NB forces, stroe NB forces in tetrads
        } else if (status.MPI_TAG == TAG_NB) {
            index = (int) temp_Forces[0][3 * max_Atoms + 1];
            io.tetrad[index].NB_Energy += temp_Forces[0][3 * max_Atoms];
            io.tetrad[index].EL_Energy += temp_Forces[1][3 * max_Atoms];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].NB_Forces[k] += temp_Forces[0][k];
            }
            
            index = (int) temp_Forces[1][3 * max_Atoms + 1];
            io.tetrad[index].NB_Energy += temp_Forces[0][3 * max_Atoms];
            io.tetrad[index].EL_Energy += temp_Forces[1][3 * max_Atoms];
            for (k = 0; k < 3 * io.tetrad[index].num_Atoms; k++) {
                io.tetrad[index].NB_Forces[k] += temp_Forces[1][k];
            }
            
        }
        
        // If there are some more need to be calculated, send indexes.
        send_Tetrad_Index(&i, &j, status.MPI_SOURCE);
        
    }
    
    // Clip NB forces & add random forces into the NB forces
    for (i  = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            io.tetrad[i].NB_Forces[j] = min( max_Forces, io.tetrad[i].NB_Forces[j]);
            io.tetrad[i].NB_Forces[j] = max(-max_Forces, io.tetrad[i].NB_Forces[j]);
        }
    }
    
}



void Master::cal_Velocities(void) {
    
    // Calculate velocities of every tetrad
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Velocities(&io.tetrad[i], i);
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
        velocities[index] += velocities[i]; coordinates[index] += coordinates[i];
        velocities[i] = velocities[index];  coordinates[i] = coordinates[index];
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
    
    int i, index = io.displs[io.crd.num_BP - 3];

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

