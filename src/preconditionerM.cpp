#include "preconditionerM.hpp"
#include "petsc_preconditionerM.hpp"

//here, only the Petsc Preconditioner is settled, why not also the Laspack one?


#include <cstdlib>  //for abort

#include "FemusExtLib_conf.hpp"

//------------------------------------------------------------------
// Preconditioner members

PreconditionerM *
PreconditionerM::build(const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {

/*
#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<PreconditionerM > ap(new LaspackPreconditionerM);
	return ap;
      }
#endif
*/

#ifdef FEMUS_HAVE_PETSC
    case PETSC_SOLVERS:
      {
	return new PetscPreconditionerM();
      }
#endif

/*
#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERS:
      {
	AutoPtr<PreconditionerM > ap(new AztecPreconditionerM);
	return ap;
      }
#endif
*/
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



