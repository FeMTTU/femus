/*=========================================================================

 Program: FEMUS
 Module: LinearEquationSolver
 Authors: Eugenio Aulisa, Simone Bnà
 
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
#include "FEMTTUConfig.h"
#include "PrecondtypeEnum.hpp"
#include "petscksp.h"
#include "petscvec.h" 

// C++ includes
#include <memory>
// Local Includes
#include "PetscLinearEquationSolver.hpp"
#include "Preconditioner.hpp"


namespace femus {



using namespace std;

// =============================================================
std::auto_ptr<LinearEquationSolver> LinearEquationSolver::build(const unsigned &igrid, mesh *other_mesh, const SolverPackage solver_package) {
  // Build the appropriate solver
  switch (solver_package)  {
#ifdef HAVE_PETSC
  case PETSC_SOLVERS:      {
    std::auto_ptr<LinearEquationSolver> ap(new PetscLinearEquationSolver(igrid, other_mesh));
    return ap;
  }
#endif
#ifdef HAVE_TRILINOS
  case TRILINOS_SOLVERS:  {
    std::auto_ptr<LinearEquationSolver> ap(new AztecLinearEquationSolver);
    return ap;


} //end namespace femus


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


// =============================================================
unsigned LinearEquation::GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol, 
				  const unsigned &idof_gmt) const {
  
   //return KKIndex[kkindex_sol]+idof_gmt;
     
   unsigned soltype =  _SolType[index_sol]; 
   unsigned isubdom = (soltype<3)?_msh->npart[idof_gmt]:(_msh->epart[idof_gmt % _msh->GetElementNumber()]);
   unsigned idof_metis = _msh->GetMetisDof(idof_gmt,soltype);   
   return KKoffset[kkindex_sol][isubdom] + idof_metis - _msh->MetisOffset[soltype][isubdom];
}





} //end namespace femus



