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

#define TAG_DATA   5555
#define TAG_TETRAD 6666
#define TAG_ED     7777
#define TAG_NB     8888
#define TAG_DEATH  9999

using namespace std;

class MPI_Library{
    
public:
    
    static void create_MPI_Tetrad(MPI_Datatype MPI_Tetrad, int num_Atoms_In_Tetrad, int num_Evecs);
    
    static void free_MPI_Tetrad(MPI_Datatype MPI_Tetrad);
    
};


#endif /* projectTest_hpp */
