/********************************************************************************
 *                                                                              *
 *          Porting the Essential Dynamics/Molecular Dynamics method            *
 *             for large-scale nucleic acid simulations to ARCHER               *
 *                                                                              *
 *                               Zhuowei Si                                     *
 *              EPCC supervisors: Elena Breitmoser, Iain Bethune                *
 *     External supervisor: Charlie Laughton (The University of Nottingham)     *
 *                                                                              *
 *                  MSc in High Performance Computing, EPCC                     *
 *                       The University of Edinburgh                            *
 *                                                                              *
 *******************************************************************************/

/**
 * File:  master.cpp
 * Brief: Implementation of the Master class functions 
 */

#include "master.hpp"


Master::Master(void) {
    
    num_Pairs = 0;
    max_Atoms = 0;
    
    comm      = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size); // Get size of MPI processes

}



Master::~Master(void) {
    
    // Deallocate memory of arrays
    array.deallocate_2D_Double_Array(pair_Lists);
    array.deallocate_2D_Int_Array(ED_Index);
    array.deallocate_2D_Int_Array(NB_Index);
    array.deallocate_2D_Double_Array(NB_Forces);
    delete [] velocities;
    delete [] coordinates;
    
    // Free the MPI_Datatype
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        mpi.free_MPI_ED_Forces(&(MPI_ED_Forces[i]));
    }
    mpi.free_MPI_Crds(&MPI_Crds);
    
}



void Master::initialise(void) {
    
    // Read the comfiguration file, the tetrad parameter file, the coordinate
    // file & Initialise the coordinates & velocities of tetrads
    io.read_Cofig(&edmd);
    io.read_Prm();
    io.read_Crd();
    io.initialise_Tetrad_Crds();
    
    // Inintialise the output frequencies
    io.ntwt -= io.ntwt % io.ntsync; if (io.ntwt == 0) io.ntwt = 1;
    io.ntpr -= io.ntpr % io.ntsync; if (io.ntpr == 0) io.ntpr = 1;
    
    // Create MPI_Datatype for message passing
    MPI_ED_Forces = new MPI_Datatype [io.prm.num_Tetrads]; // For every tetrad
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        if(max_Atoms < io.tetrad[i].num_Atoms) max_Atoms = io.tetrad[i].num_Atoms;
        mpi.create_MPI_ED_Forces(&(MPI_ED_Forces[i]), &(io.tetrad[i]));
    }
    mpi.create_MPI_Crds(&MPI_Crds, io.prm.num_Tetrads, io.tetrad);// For all tetrads
    
    // Allocate memory for arrays
    num_Pairs  = io.prm.num_Tetrads * (io.prm.num_Tetrads - 1) / 2;
    pair_Lists = array.allocate_2D_Double_Array(num_Pairs, 2);
    ED_Index = array.allocate_2D_Int_Array(size - 1, 2);
    NB_Index = array.allocate_2D_Int_Array(size - 1, 2);
    NB_Forces  = array.allocate_2D_Double_Array(io.prm.num_Tetrads, 3 * max_Atoms + 2);
    velocities  = new double [3 * io.crd.total_Atoms];
    coordinates = new double [3 * io.crd.total_Atoms];
    
    // Print information of the EDMD simulation
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
    double edmd_Para[11] = { edmd.dt, edmd.gamma, edmd.tautp, edmd.temperature,
        edmd.scaled, edmd.mole_Cutoff, edmd.atom_Cutoff, edmd.mole_Least,
        (double)io.prm.num_Tetrads, (double)num_Pairs, (double)max_Atoms };
    
    // Assign the number of atoms & evecs of tetrads into the sending array
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        tetrad_Para[2 * i] = io.tetrad[i].num_Atoms;
        tetrad_Para[2*i+1] = io.tetrad[i].num_Evecs;
    }
    
    // Broadcast the simulation and tetrad parameters
    cout << "Master sending parameters to workers..." << endl;
    MPI_Bcast(edmd_Para, 11, MPI_DOUBLE, 0, comm);
    MPI_Bcast(tetrad_Para, 2 * io.prm.num_Tetrads, MPI_INT, 0, comm);
    
    delete [] tetrad_Para;
    
}



void Master::send_Tetrads(void) {
    
    MPI_Datatype MPI_Tetrad;
    
    cout << "Master sending tetrads to workers...\n" << endl;
    
    // Create MPI_Datatype "MPI_Tetrad" for all tetrads
    mpi.create_MPI_Tetrad(&MPI_Tetrad, io.prm.num_Tetrads, io.tetrad);
    
    // Broadcast tetrads to all workers
    MPI_Bcast(io.tetrad, 1, MPI_Tetrad, 0, comm);
    
    // Free the MPI_Datatype "MPI_Tetrad"
    mpi.free_MPI_Tetrad(&MPI_Tetrad);
    
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
    double r, ** com = array.allocate_2D_Double_Array(io.prm.num_Tetrads, 3);
    
    cal_Centre_of_Mass(com);
    
    // Loop to generate pair lists
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
    
    array.deallocate_2D_Double_Array(com);
    
}



void Master::generate_Indexes(void) {
    
    int i, loop;
    
    // Divide NB force calculation into simuilar chunk (balanced workload)
    // For every part of the NB force caulation, it has the start point &
    // how many NB forces to be calculated
    
    // Distribute the workload to workers in balance
    for (i = 0; i < size - 1; i++) {
        ED_Index[i][1] = io.prm.num_Tetrads / (size - 1);
        NB_Index[i][1] = num_Pairs / (size - 1);
    }
    
    // If can not be divided exactly, then the remaining works are assigned
    // to parts of the workers.
    loop = io.prm.num_Tetrads - (ED_Index[0][1] * (size - 1));
    for (i = 0; i < loop; i++) { ED_Index[i][1] += 1; }
    loop = num_Pairs - (NB_Index[0][1] * (size - 1));
    for (i = 0; i < loop; i++) { NB_Index[i][1] += 1; }
    
    // Set the start index of the workload
    ED_Index[0][0] = NB_Index[0][0] = 0;
    for (i = 1; i < size - 1; i++) {
        ED_Index[i][0] = ED_Index[i - 1][0] + ED_Index[i - 1][1];
        NB_Index[i][0] = NB_Index[i - 1][0] + NB_Index[i - 1][1];
    }
    
}


void Master::send_Workload_Indexes(void) {
    
    MPI_Request send_Request[3][size - 1];
    MPI_Status send_Status[3][size - 1];
    
    // Send the pair lists & the workload displacements to workers
    for (int i = 0; i < size - 1; i++) {
        MPI_Isend(&(pair_Lists[0][0]), 2 * num_Pairs, MPI_DOUBLE, i + 1,
                  TAG_PAIRS,     comm, &(send_Request[0][i]));
        MPI_Isend(&(NB_Index[0][0]), 2 * (size - 1),  MPI_INT,    i + 1,
                  TAG_PAIRS + 1, comm, &(send_Request[1][i]));
        MPI_Isend(&(ED_Index[0][0]), 2 * (size - 1),  MPI_INT,    i + 1,
                  TAG_PAIRS + 2, comm, &(send_Request[2][i]));
    }
    MPI_Waitall(size - 1, send_Request[0], send_Status[0]);
    MPI_Waitall(size - 1, send_Request[1], send_Status[1]);
    MPI_Waitall(size - 1, send_Request[2], send_Status[2]);
    
}




void Master::calculate_Forces(void) {
    
    int i, j, rank, signal = TAG_FORCE;
    MPI_Request send_Request[size - 1], recv_Request[io.prm.num_Tetrads];
    MPI_Status send_Status[size - 1], recv_Status[io.prm.num_Tetrads];
    
    // Send a signle to indicate workers to prepare the force calculation
    // Broadcast the cooridnates
    for (i = 0; i < size - 1; i++) {
        MPI_Isend(&signal, 1, MPI_INT, i + 1, TAG_FORCE, comm, &(send_Request[i]));
    }
    MPI_Bcast(io.tetrad, 1, MPI_Crds, 0, comm);
    MPI_Waitall(size - 1, send_Request, send_Status);
    
    // Receive all the ED forces from workers
    for (i = 0; i < io.prm.num_Tetrads; i++) {
        MPI_Irecv(&(io.tetrad[i]), 1, MPI_ED_Forces[i], MPI_ANY_SOURCE, TAG_ED + i, comm, &(recv_Request[i]));
    }
    
    // Reduce & sum up the NB forces & process and assign the NB forces to tetrads
    MPI_Reduce(MPI_IN_PLACE, &(NB_Forces[0][0]), io.prm.num_Tetrads * (3 * max_Atoms + 2), MPI_DOUBLE, MPI_SUM, 0, comm);
    process_NB_Forces();
    MPI_Waitall(io.prm.num_Tetrads, recv_Request, recv_Status);
    
}




void Master::process_NB_Forces(void) {
    
    int i, j;
    for (i  = 0; i < io.prm.num_Tetrads; i++) {
        for (j = 0; j < 3 * io.tetrad[i].num_Atoms; j++) {
            
            // Clip the NB forces between -1.0 and 1.0
            if (NB_Forces[i][j] < -1.0)  NB_Forces[i][j]  =  -1.0;
            else if (NB_Forces[i][j] > 1.0) NB_Forces[i][j] = 1.0;
            
            // Assign NB forces to tetrads
            io.tetrad[i].NB_Forces[j] = NB_Forces[i][j];
            
            // Clear the NB force array for next use
            NB_Forces[i][j] = 0.0;
            
        }
        
        // NB energy & Electrostatic Energy
        io.tetrad[i].NB_Forces[j]     = NB_Forces[i][j];
        io.tetrad[i].NB_Forces[j + 1] = NB_Forces[i][j + 1];
        
        NB_Forces[i][j] = NB_Forces[i][j + 1] = 0.0;
    }
    
}




void Master::update_Velocity(void) {
    
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Velocities(&io.tetrad[i]);
    }
    
}



void Master::update_Coordinate(void) {
    
    for (int i = 0; i < io.prm.num_Tetrads; i++) {
        edmd.update_Coordinates(&io.tetrad[i]);
    }

}



void Master::merge_Vels_n_Crds(void) {
    
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
    
    // Process the first & last 3 tetrads
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



void Master::write_Info(int istep) {
    
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
    io.write_Energies(istep + io.ntsync, energies);
    
    // Write trajectory
    // io.displs[io.crd.num_BP - 3] is the index that the afterwards 3
    // tetrads the same as the first three tetrads
    io.write_Trajectory(istep, io.displs[io.crd.num_BP - 3], coordinates);
    
}



void Master::write_Crds(void) {
    
    io.update_Crd_File(velocities, coordinates);
    
}



void Master::finalise(void) {
    
    MPI_Request send_Request[size - 1];
    MPI_Status send_Status[size - 1];
    
    // Send the terminate signal to all workers
    for (int signal = TAG_END, i = 0; i < size - 1; i++) {
        MPI_Isend(&signal, 1, MPI_INT, i + 1, TAG_END, comm, &(send_Request[i]));
    }
    MPI_Waitall(size - 1, send_Request, send_Status);
    
}

