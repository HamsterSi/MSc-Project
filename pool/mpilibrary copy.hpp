//
//  projectTest.hpp
//  
//
//  Created by Zhuowei Si on 04/04/2016.
//
//

#ifndef projectTest_hpp
#define projectTest_hpp

// MPI P2P tag to use for command communications, it is important not to reuse this
#define CONTROL_TAG 16384
#define PID_TAG 16383

// Pool options
#define QuitOnNoProcs 1
#define IgnoreOnNoProcs 0
#define DEBUG 0

// The process status
#define STOP        0
#define SLEEPING    1
#define WAKE        2
#define START       3
#define COMPLETE    4

class MPI_Library{
  
private:
    
    // Internal pool global state
    int my_Rank;
    int my_Size;
    char* active = NULL;
    int active_Workers;
    int waiting_Start;
    
    
    
    static int myRank;
    static int numProcs;
    static char* active=NULL;
    static int processesAwaitingStart;
    static struct Control_Package in_command;
    static MPI_Request pollRecvCommandRequest = MPI_REQUEST_NULL;
    
    // Internal pool functions
    static void errorMessage(char*);
    static int startAwaitingProcessesIfNeeded(int, int);
    static int handleRecievedCommand();
    static void initialiseType();
    static struct Control_Package createCommandPackage(enum process_Status);
    
public:
    // Setup MPI library.
    static void initialize(int* argc, char **argv[]);
    
    static int get_Size(int* size);
    
    static int get_Rank(int* rank);
    
    static void finalize(void);
    
    // Initialises the process pool
    int processPoolInit();
    
    // Finalises the process pool
    void processPoolFinalise();
    
    // Called by the master in loop, blocks until state change, 1=continue and 0=stop
    int masterPoll();
    
    // Called by worker once task completed and will sleep until it receives instructions, 1=continue with new task and 0=stop
    int workerSleep();
    
    // Determines whether the current worker should stop or not (i.e. whether the pool is shutting down)
    int shouldWorkerStop();
    
    // Called by the master or a worker to start a new worker process
    int startWorkerProcess();
    
    // Called by a worker to shut the pool down
    void shutdownPool();
    
    // Retrieves the optional data associated with the command, provides an example of how this can be done
    int getCommandData();
    
    // Returns the number of active workers
    int getactive_Workers();

    
};


#endif /* projectTest_hpp */
