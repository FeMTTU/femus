#ifndef _femus_init_
#define _femus_init_

#include "FEMTTUConfig.h"


#ifdef HAVE_MPI
#include "mpi.h"  //For MPI_COMM_WORLD
#endif

class FemusInit { 

public:

FemusInit(int & argc, char** & argv/*, MPI_Comm comm_world_in = MPI_COMM_WORLD*/);

~FemusInit();



};



#endif