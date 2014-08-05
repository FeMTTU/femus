/*=========================================================================

 Program: FEMuS
 Module: LinearEquationSolver
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMuS
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __linear_equation_solver_h__
#define __linear_equation_solver_h__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FEMTTUConfig.h"
#include "SolverPackageEnum.hpp"
#include "PrecondtypeEnum.hpp"
#include "SolvertypeEnum.hpp"
#include "LinearEquation.hpp"
#include "MgSmootherEnum.hpp"
#include "vector"
#include <memory>
#include <cstdio>


namespace femus {


using std::vector;

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class SparseMatrix;
class NumericVector;
class Preconditioner;

/** 
* This class provides a uniform interface for linear solvers.  This base
* class is overloaded to provide linear solvers from different packages
* like PETSC or LASPACK.
*/ 

class LinearEquationSolver : public LinearEquation {
  // ================================
  // DATA
  // ==============================
protected:
  /// Enum stating which type of iterative solver to use.
  SolverType _solver_type;

  /// Enum statitng with type of preconditioner to use.
  PreconditionerType _preconditioner_type;

  /// Holds the Preconditioner object to be used for the linear solves.
  Preconditioner *_preconditioner;

  /// Flag indicating if the data structures have been initialized.
  bool _is_initialized;

  /// The number of elements in a Vanka Block
  //unsigned _num_elem_vanka_block;

public:

  /// Boolean flag to indicate whether we want to use an identical preconditioner to the previous solve.
  bool same_preconditioner;

  // =================================
  // CONSTR/DESTR
  // ================================
  /**  Constructor. Initializes Solver data structure */
  LinearEquationSolver (const unsigned &igrid, mesh* other_msh);
  /// Destructor.
  virtual ~LinearEquationSolver();

  /// Release all memory and clear data structures.
  virtual void clear() {}

  /// Initialize data structures if not done so already.
  //virtual void init() = 0;

  /** Builds a \p LinearEquationSolver using the linear solver in \p solver_package */
  static std::auto_ptr<LinearEquationSolver> build(const unsigned &igrid, mesh* other_msh,const MgSmoother & smoother_type, const SolverPackage solver_package =LSOLVER);
  
  // =================================
  // SETTING FUNCTIONS
  // ================================
  ///Set the tolerance for the solver
  virtual void set_tolerances(const double &rtol, const double &atol,
                              const double &divtol, const unsigned &maxits) = 0;

  /// Sets the type of solver to use.
  void set_solver_type (const SolverType st)  {
    _solver_type = st;
  }

  /// Sets the type of preconditioner to use.
  void set_preconditioner_type (const PreconditionerType pct);

  /// Attaches a Preconditioner object to be used
  void attach_preconditioner(Preconditioner * preconditioner);

  // =================================
  // RETURN FUNCTIONS
  // ================================

  /// @returns true if the data structures are
  bool initialized () const {
    return _is_initialized;
  }

  /// Returns the type of solver to use.
  SolverType solver_type () const {
    return _solver_type;
  }

  /// Returns the type of preconditioner to use.
  PreconditionerType preconditioner_type () const;

  // =================================
  // OUR virtual FUNCTIONS
  // ================================

  /// Set the number of elements of the Vanka Block
  virtual void SetElementBlockNumber(const unsigned & block_elemet_number){
    std::cout<<"Warning SetElementBlockNumber(const unsigned &) is not available for this smoother\n"; 
  };
  
  virtual void SetElementBlockNumber(const char all[], const unsigned & overlap=1){
    std::cout<<"Warning SetElementBlockNumber(const char [], const unsigned & ) is not available for this smoother\n"; 
  };
  
  
  virtual void SetNumberOfSchurVariables(const unsigned short & NSchurVar){
    std::cout<<"Warning SetNumberOfSchurVariables(const unsigned short &) is not available for this smoother\n"; 
  };
  
  virtual void SetDirichletBCsHandling(const unsigned int &DirichletBCsHandlingMode){
    std::cout<<"Warning SetDirichletBCsHandling(const unsigned int &) is not available for this smoother\n"; 
  };
  
  /// Call the Vanka(Schur) smoother-solver using the PetscLibrary.
  virtual std::pair< int, double> solve(const vector <unsigned> &VankaIndex, const bool &ksp_clean) = 0;
  /// Call the Gmres smoother-solver
  // virtual std::pair< int, double> solve( const bool &clean = true) = 0;
  
//   /// This function calls the solver
//   virtual std::pair< int, double> solve (SparseMatrix&,  // System Matrix
//                                          NumericVector&, // Solution vector
//                                          NumericVector&, // RHS vector
//                                          const double,      // Stopping tolerance
//                                          const  int) = 0; // N. Iterations
// 
//   ///   *This function calls the solver with preconditioner
//   virtual std::pair< int, double> solve (SparseMatrix&,  // System Matrix
//                                          SparseMatrix&,  // Preconditioning Matrix
//                                          NumericVector&, // Solution vector
//                                          NumericVector&, // RHS vector
//                                          const double,      // Stopping tolerance
//                                          const  int) = 0; // N. Iteration

  /// Prints a useful message about why the latest linear solve con(di)verged.
 // virtual void print_converged_reason() = 0;
};

// -------------------- inline functions ---------

// =============================================
inline LinearEquationSolver::LinearEquationSolver(const unsigned &igrid, mesh* other_msh) :
  LinearEquation(other_msh),
  _solver_type(GMRES),
  _preconditioner(NULL),
  _is_initialized(false),
  same_preconditioner(false) {

  if(igrid==0){
    _preconditioner_type=LU_PRECOND;
  }
  else{
    _preconditioner_type=ILU_PRECOND;
  }
}

inline LinearEquationSolver::~LinearEquationSolver() {
  this->clear();
}


} //end namespace femus



#endif 
