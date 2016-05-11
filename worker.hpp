
#ifndef worker_hpp
#define worker_hpp

#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "mpi.h"

#include "pool/pool.h"
#include "mpilibrary.hpp"
#include "edmd.hpp"
#include "tetrad.hpp"
#include "io.hpp"

using namespace std;

class Worker_Management {
    
public:
    
    int worker_ED_Forces(void);
    
    int worker_NB_Forces(void);
    
};


#endif /* worker_hpp */
