//
//  projectTest.hpp
//  
//
//  Created by Zhuowei Si on 04/04/2016.
//
//

#ifndef projectTest_hpp
#define projectTest_hpp

#include <iostream>
#include <cstddef>
#include "mpi.h"

#include "tetrad.hpp"

#define TAG_ED 8888
#define TAG_NB 9999s

using namespace std;

class MPI_Library{
    
public:
    // Setup MPI library.
    static void initialize(int* argc, char **argv[]);
    
    static int get_Size(int* size);
    
    static int get_Rank(int* rank);
    
    static void create_MPI_Tetrad(MPI_Datatype MPI_Tetrad, int num_Atoms_In_Tetrad, int num_Evecs);
    
    static void free_MPI_Tetrad(MPI_Datatype MPI_Tetrad);
    
    static void finalize(void);
    
};


#endif /* projectTest_hpp */
