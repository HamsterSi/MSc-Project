//
//  projectTest.cpp
//  
//
//  Created by Zhuowei Si on 04/04/2016.
//
//

#include <iostream>
#include "mpi.h"

#include "./mpilibrary.hpp"


/* Initialise the MPI library */
void MPI_Library::initialize(int* argc, char **argv[]){
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &my_Size);
    
    if (my_Rank == 0) {
        if (my_Size < 2) {
            cout << "Not more worker process available" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);  
        }
        active = (char*) malloc(my_Size);
        active_Workers = 0;
        waiting_Start  = 0;
        for (int i = 0; i < my_Size-1; i++) active[i] = 0;
        if (DEBUG) cout << "[Master] Initialised Master" << endl;
        return 0;
    } else {
        int worker_Status;
        MPI_Recv(&worker_Status, 1, MPI_INT, 0, CONTROL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return handleRecievedCommand();
    }
}

int MPI_Library::get_Size(int* size){
    
    MPI_Comm_size(MPI_COMM_WORLD, size);
    return *size;
}

int MPI_Library::get_Rank(int* rank){
    
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    return *rank;
}

/* Finalize the MPI library */
void MPI_Library::finalise(void) {
    if (my_Rank == 0) {
        if (active != NULL) free(active);
        for (int i = 0; i < my_Size-1; i++) {
            if (DEBUG) cout << "[Master] Shutting down process" << i << endl;
            int worker_Status = STOP;
            MPI_Send(&worker_Status, 1, MPI_INT, i+1, CONTROL_TAG, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}


int MPI_Library::master_Poll(void) {
    if (my_Rank == 0) {
        MPI_Status status;
        int worker_Status;
        MPI_Recv(&worker_Status, 1, MPI_INT, MPI_ANY_SOURCE, CONTROL_TAG, MPI_COMM_WORLD, &status);
        
        switch (worker_Status) {
            case SLEEPING:
                if (DEBUG) cout << "[Master] Received sleep command from " << status.MPI_SOURCE << endl;
                active_Workers--;
                active[status.MPI_SOURCE-1]=0;
                break;
                
            case COMPLETE:
                if (DEBUG) cout << "[Master] Received shutdown command" << endl;
                return 0;
                break;
                
            case START:
                waiting_Start++;
                int returnRank = startAwaitingProcessesIfNeeded(processesAwaitingStart, status.MPI_SOURCE);
                
                int i, rank, parent = status.MPI_SOURCE;
                int waiting_Id = waiting_Start;
                if (waiting_Start) {
                    for (i = 0; i < my_Size-1; i++) {
                        if (!active[i]) {
                            active_Workers++;
                            active[i] = 1;
                            worker_Status = WAKE;
                            waiting_Id = waiting_Start ? parent : -1;
                            if (DEBUG) cout << "[Master] Starting process" << i+1 << endl;
                            MPI_Send(&worker_Status, 1, MPI_INT, i+1, CONTROL_TAG, MPI_COMM_WORLD);
                        }
                    }
                }
                
                MPI_Send(&returnRank, 1, MPI_INT, status.MPI_SOURCE, PID_TAG, MPI_COMM_WORLD);
                break;
                
            default:
                break;
        }
        return 1;
        
    } else {
        cout << "Worker process called master poll" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 0;
    }
}

static int startAwaitingProcessesIfNeeded(int awaitingId, int parent) {
    int awaitingProcessMPIRank=-1;
    if(processesAwaitingStart) {
        int i;
        for(i=0;i<numProcs-1;i++){
            if(!active[i]) {
                active_Workers++;
                active[i]=1;
                struct Control_Package out_command = createCommandPackage(WAKE);
                out_command.data = awaitingId == processesAwaitingStart ? parent : -1;
                if (DEBUG) printf("[Master] Starting process %d\n", i+1);
                MPI_Send(&out_command, 1, COMMAND_TYPE, i+1, CONTROL_TAG, MPI_COMM_WORLD);
                if (awaitingId == processesAwaitingStart) awaitingProcessMPIRank = i+1;	// Will return this rank to the caller
                processesAwaitingStart--;
                if (processesAwaitingStart == 0) break;
            }
            if(i==numProcs-2){
                //If I reach this point, I must have looped through the whole array and found no available processes
                if(QuitOnNoProcs) {
                    errorMessage("No more processes available");
                }
                
                if(IgnoreOnNoProcs) {
                    fprintf(stderr,"[ProcessPool] Warning. No processes available. Ignoring launch request.\n");
                    processesAwaitingStart--;
                }
                // otherwise, do nothing; a process may become available on the next iteration of the loop
            }
        }
    }
    return awaitingProcessMPIRank;
}

int startWorkerProcess() {
    if (myRank == 0) {
        processesAwaitingStart++;
        return startAwaitingProcessesIfNeeded(processesAwaitingStart, 0);
    } else {
        int workerRank;
        struct Control_Package out_command = createCommandPackage(START);
        MPI_Send(&out_command, 1, COMMAND_TYPE, 0, CONTROL_TAG, MPI_COMM_WORLD);
        // Receive the rank that this worker has been placed on - if you change the default option from aborting when
        // there are not enough MPI processes then this may be -1
        MPI_Recv(&workerRank, 1, MPI_INT, 0, PID_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return workerRank;
    }
}

void shutdownPool() {
    if (myRank != 0) {
        if (DEBUG) printf("[Worker] Commanding a pool shutdown\n");
        struct Control_Package out_command = createCommandPackage(COMPLETE);
        MPI_Send(&out_command, 1, COMMAND_TYPE, 0, CONTROL_TAG, MPI_COMM_WORLD);
    }
}

int workerSleep() {
    if (myRank != 0) {
        if (in_command.command==WAKE) {
            // The command was to wake up, it has done the work and now it needs to switch to sleeping mode
            struct Control_Package out_command = createCommandPackage(SLEEPING);
            MPI_Send(&out_command, 1, COMMAND_TYPE, 0, CONTROL_TAG, MPI_COMM_WORLD);
            if (pollRecvCommandRequest != MPI_REQUEST_NULL) MPI_Wait(&pollRecvCommandRequest, MPI_STATUS_IGNORE);
        }
        return handleRecievedCommand();
    } else {
        errorMessage("Master process called worker poll");
        return 0;
    }
}


int getCommandData() {
    return in_command.data;
}

int getactive_Workers()
{
    return active_Workers;
}

int shouldWorkerStop() {
    if (pollRecvCommandRequest != MPI_REQUEST_NULL) {
        int flag;
        MPI_Test(&pollRecvCommandRequest, &flag, MPI_STATUS_IGNORE);
        if(flag && in_command.command==STOP){
            // If there's a message waiting, and it's a stop call then return 1 to denote stop
            return 1;
        }
    }
    return 0;
}

static int handleRecievedCommand() {
    // We have just (most likely) received a command, therefore decide what to do
    if (in_command.command==WAKE) {
        // If we are told to wake then post a recv for the next command and return true to continues
        MPI_Irecv(&in_command, 1, COMMAND_TYPE, 0, CONTROL_TAG, MPI_COMM_WORLD, &pollRecvCommandRequest);
        if (DEBUG) printf("[Worker] Process %d woken to work\n", myRank);
        return 1;
    } else if (in_command.command == STOP) {
        // Stopping so return zero to denote stop
        if (DEBUG) printf("[Worker] Process %d commanded to stop\n", myRank);
        return 0;
    } else {
        errorMessage("Unexpected control command");
        return 0;
    }
}

static void errorMessage(char * message) {
    fprintf(stderr,"%4d: [ProcessPool] %s\n", myRank, message);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

static struct Control_Package createCommandPackage(enum process_Status desiredCommand) {
    struct Control_Package package;
    package.command = desiredCommand;
    return package;
}
