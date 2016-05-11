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
#include "mpi.h"

using namespace std;

class MPI_Library{
    
public:
    // Setup MPI library.
    static void initialize(int* argc, char **argv[]);
    
    static int get_Size(int* size);
    
    static int get_Rank(int* rank);
    
    static void finalize(void);
    
};


#endif /* projectTest_hpp */
