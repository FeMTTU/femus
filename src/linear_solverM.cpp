#include "linear_solverM.hpp"

#include "FemusExtLib_conf.hpp"

#include "PrecondtypeEnum.hpp"

// C++ includes
#include <memory>
#include <cstdlib>

#include "petsc_linear_solverM.hpp"
#include "Preconditioner.hpp"

//------------------------------------------------------------------
// LinearSolver members

// =============================================================
std::auto_ptr<LinearSolverM > LinearSolverM::build(const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)  {
#ifdef FEMUS_HAVE_PETSC
    case PETSC_SOLVERS:      {
	std::auto_ptr<LinearSolverM > ap(new PetscLinearSolverM);
	return ap;
      }
#endif
#ifdef HAVE_TRILINOS
    case TRILINOS_SOLVERSM:  {
	std::auto_ptr<LinearSolverM > ap(new AztecLinearSolverM);
	return ap;
      }
#endif

    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package<< std::endl;
      abort();
    }
    
  std::auto_ptr<LinearSolverM > ap(NULL);
  return ap;    
}

// ============================================================
PreconditionerType LinearSolverM::preconditioner_type () const
{
  if(_preconditioner)    return _preconditioner->type();
  return _preconditioner_type;
}

// ===========================================================
void LinearSolverM::set_preconditioner_type (const PreconditionerType pct)
{
  if(_preconditioner)    _preconditioner->set_type(pct);
  else    _preconditioner_type = pct;
}

// =============================================================
void LinearSolverM::attach_preconditioner(Preconditioner * preconditioner)
{
  if(this->_is_initialized)  {
    std::cerr<<"Preconditioner must be attached before the solver is initialized!"<<std::endl;
    abort();
  } 
  _preconditioner_type = SHELL_PRECOND;
  _preconditioner = preconditioner;
}

//------------------------------------------------------------------
// Explicit instantiations
// template class LinearSolver<Number>;



