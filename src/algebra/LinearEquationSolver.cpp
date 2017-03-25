/*=========================================================================

  Program: FEMUS
  Module: LinearEquationSolver
  Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/


//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "LinearEquationSolver.hpp"
#include "FemusConfig.hpp"
#include "PrecondtypeEnum.hpp"
#include "petscksp.h"
#include "petscvec.h"

// Local Includes
#include "AsmPetscLinearEquationSolver.hpp"
#include "GmresPetscLinearEquationSolver.hpp"
#include "FieldSplitPetscLinearEquationSolver.hpp"
#include "Preconditioner.hpp"

namespace femus {

  // =============================================================
  std::auto_ptr<LinearEquationSolver> LinearEquationSolver::build(const unsigned &igrid, Solution *other_solution,const MgSmoother & smoother_type, const SolverPackage solver_package) {
    // Build the appropriate solver

    switch (solver_package)  {
#ifdef HAVE_PETSC
    case PETSC_SOLVERS:  {
      switch(smoother_type){
      case ASM_SMOOTHER:{
	std::auto_ptr<LinearEquationSolver> ap(new AsmPetscLinearEquationSolver(igrid, other_solution));
	return ap;
      }
      case GMRES_SMOOTHER:{
      	std::auto_ptr<LinearEquationSolver> ap(new GmresPetscLinearEquationSolver(igrid, other_solution));
	return ap;
      }
      case FIELDSPLIT_SMOOTHER:{
        std::auto_ptr<LinearEquationSolver> ap(new FieldSplitPetscLinearEquationSolver(igrid, other_solution));
        return ap;
      }
      }
    }
#endif
#ifdef HAVE_TRILINOS
    case TRILINOS_SOLVERS:  {
      std::auto_ptr<LinearEquationSolver> ap(new AztecLinearEquationSolver);
      return ap;
    }
#endif
    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package<< std::endl;
      abort();
    }

    std::auto_ptr<LinearEquationSolver> ap(NULL);
    return ap;
  }

  // ============================================================
  PreconditionerType LinearEquationSolver::preconditioner_type () const {
    if (_preconditioner)    return _preconditioner->type();
    return _preconditioner_type;
  }

  // ===========================================================
  void LinearEquationSolver::set_preconditioner_type (const PreconditionerType pct) {
    if (_preconditioner)    _preconditioner->set_type(pct);
    else    _preconditioner_type = pct;
  }

  // =============================================================
  void LinearEquationSolver::attach_preconditioner(Preconditioner * preconditioner) {
    if (this->_is_initialized)  {
      std::cerr<<"Preconditioner must be attached before the solver is initialized!"<<std::endl;
      abort();
    }
    _preconditioner_type = SHELL_PRECOND;
    _preconditioner = preconditioner;
  }

} //end namespace femus



