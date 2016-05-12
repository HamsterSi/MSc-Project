
#ifndef simulation_hpp
#define simulation_hpp

#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "mpi.h"

#include "pool/pool.h"
#include "mpilibrary.hpp"
#include "master.hpp"
#include "worker.hpp"

using namespace std;

void master_Code(void);

void worker_Code(void);

#endif /* simulation_hpp */