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

#define TAG_TETRAD 5
#define TAG_DATA   1
#define TAG_ED     2
#define TAG_NB     3
#define TAG_DEATH  4

using namespace std;

class MPI_Library{
    
public:
    
    static void create_MPI_Tetrad(MPI_Datatype* MPI_Tetrad, Tetrad* tetrad);
    
    static void free_MPI_Tetrad(MPI_Datatype* MPI_Tetrad);
    
};


#endif /* projectTest_hpp */

