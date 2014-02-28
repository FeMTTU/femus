#ifndef _femus_init_
#define _femus_init_

#include "FemusExtLib_conf.hpp"

// #ifndef FEMUS_HAVE_LASPACK  //use this class only with petsc-mpi //NO, use it always, and put the switches inside

#ifdef HAVE_MPI
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "mpi.h"  //For MPI_COMM_WORLD
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif

class FemusInit { 

public:

FemusInit(int & argc, char** & argv/*, MPI_Comm comm_world_in = MPI_COMM_WORLD*/);

~FemusInit();



};

// #endif


#endif