#include "LinearSolver.hpp"
#include "FEMTTUConfig.h"
#include "PrecondtypeEnum.hpp"
#include "petscksp.h"
#include "petscvec.h" 

// C++ includes
#include <memory>
// Local Includes
#include "PetscLinearSolver.hpp"
#include "Preconditioner.hpp"

using namespace std;

vector <unsigned > LinearSolverM::indexa;
vector <unsigned > LinearSolverM::indexb;
vector <unsigned > LinearSolverM::indexc;
vector <unsigned > LinearSolverM::indexd;

vector <PetscInt> LinearSolverM::indexai;
vector <PetscInt> LinearSolverM::indexbi;
vector <PetscInt> LinearSolverM::indexci;
vector <PetscInt> LinearSolverM::indexdi;

//------------------------------------------------------------------
// LinearSolver members

// =============================================================
// std::auto_ptr<LinearSolverM > LinearSolverM::build(vector < vector < double> > &vt,
//     const double Lref, const char infile[], const SolverPackage solver_package) {
//   // Build the appropriate solver
//   switch (solver_package)  {
// #if HAVE_PETSC == 1
//   case PETSC_SOLVERS:      {
//     std::auto_ptr<LinearSolverM > ap(new PetscLinearSolver(infile,vt,Lref));
//     return ap;
//   }
// #endif
// #if HAVE_TRILINOS == 1
//   case TRILINOS_SOLVERS:  {
//     std::auto_ptr<LinearSolverM > ap(new AztecLinearSolverM);
//     return ap;
//   }
// #endif
// 
//   default:
//     std::cerr << "ERROR:  Unrecognized solver package: "
//               << solver_package<< std::endl;
//     abort();
//   }
// 
//   std::auto_ptr<LinearSolverM > ap(NULL);
//   return ap;
// }
// 
// // =============================================================
// std::auto_ptr<LinearSolverM > LinearSolverM::build(const unsigned &igrid,elem *elc, const SolverPackage solver_package) {
//   Build the appropriate solver
//   switch (solver_package)  {
// #if HAVE_PETSC == 1
//   case PETSC_SOLVERS:      {
//     std::auto_ptr<LinearSolverM > ap(new PetscLinearSolver(igrid,elc));
//     return ap;
//   }
// #endif
// #if HAVE_TRILINOS == 1
//   case TRILINOS_SOLVERS:  {
//     std::auto_ptr<LinearSolverM > ap(new AztecLinearSolverM);
//     return ap;
//   }
// #endif
// 
//   default:
//     std::cerr << "ERROR:  Unrecognized solver package: "
//               << solver_package<< std::endl;
//     abort();
//   }
// 
//   std::auto_ptr<LinearSolverM > ap(NULL);
//   return ap;
// }


// =============================================================
std::auto_ptr<LinearSolverM > LinearSolverM::build(const unsigned &igrid, mesh *other_mesh, const SolverPackage solver_package) {
  // Build the appropriate solver
  switch (solver_package)  {
#if HAVE_PETSC == 1
  case PETSC_SOLVERS:      {
    std::auto_ptr<LinearSolverM > ap(new PetscLinearSolver(igrid, other_mesh));
    return ap;
  }
#endif
#if HAVE_TRILINOS == 1
  case TRILINOS_SOLVERS:  {
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


void LinearSolverM::set_num_elem_vanka_block(const unsigned num_elem_vanka_block) {
  _num_elem_vanka_block = num_elem_vanka_block;
};


// ============================================================
PreconditionerType LinearSolverM::preconditioner_type () const {
  if (_preconditioner)    return _preconditioner->type();
  return _preconditioner_type;
}

// ===========================================================
void LinearSolverM::set_preconditioner_type (const PreconditionerType pct) {
  if (_preconditioner)    _preconditioner->set_type(pct);
  else    _preconditioner_type = pct;
}

// =============================================================
void LinearSolverM::attach_preconditioner(Preconditioner * preconditioner) {
  if (this->_is_initialized)  {
    std::cerr<<"Preconditioner must be attached before the solver is initialized!"<<std::endl;
    abort();
  }
  _preconditioner_type = SHELL_PRECOND;
  _preconditioner = preconditioner;
}

//------------------------------------------------------------------
// Explicit instantiations
// template class LinearSolver<Number>;



