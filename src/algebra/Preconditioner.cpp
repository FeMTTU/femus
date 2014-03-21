/*=========================================================================

 Program: FEMUS
 Module: Preconditioner
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Preconditioner.hpp"
#include "PetscPreconditioner.hpp"
#include <cstdlib>  
#include "FEMTTUConfig.h"

//------------------------------------------------------------------
// Preconditioner members

Preconditioner *
Preconditioner::build(const SolverPackage solver_package) {
  // Build the appropriate solver
  switch (solver_package) {

#if HAVE_PETSC == 1
  case PETSC_SOLVERS: {
    return new PetscPreconditioner();
  }
#endif
#if HAVE_TRILINOS == 1
      case TRILINOS_SOLVERS:
        {
  	AutoPtr<Preconditioner > ap(new AztecPreconditioner);
  	return ap;
        }
#endif

  default:
    std::cerr << "ERROR:  Unrecognized solver package: "
              << solver_package
              << std::endl;
    abort();
  }

  return NULL;
}



//------------------------------------------------------------------
// Explicit instantiations
// template class Preconditioner<Number>;



