#ifndef _femttu_init_
#define _femttu_init_

#include "FEMTTUConfig.h"

//use this class only with petsc-mpi
#ifdef HAVE_MPI
#include <mpi.h>  //For MPI_COMM_WORLD
#endif

// PETSC ----------------------
#ifdef HAVE_PETSC
// # include "petsc_macroM.hpp"
// EXTERN_C_FOR_PETSC_BEGIN
# include <petsc.h>
# include <petscerror.h>
// EXTERN_C_FOR_PETSC_END
#endif

// ========================================
/// Class MGFemusInit
// ========================================

class FemTTUInit { 
public:
/// Constructor
 FemTTUInit(int & argc,char** & argv,MPI_Comm comm_world_in=MPI_COMM_WORLD);
/// Destructor
 ~FemTTUInit();
};

#endif // end _femus_init ----------------------