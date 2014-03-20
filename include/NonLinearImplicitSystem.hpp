/*=========================================================================

 Program: FEMUS
 Module: NonLinearImplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __nonlinear_implicit_system_h_
#define __nonlinear_implicit_system_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ImplicitSystem.hpp"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------


class NonLinearImplicitSystem : public ImplicitSystem {

public:

/** Constructor.  Optionally initializes required data structures. */
  NonLinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number);
  
  virtual ~NonLinearImplicitSystem();
  
  /** Solves the system. */
  virtual void solve ();
  
  /** Clear all the data structures associated with the system. */
  virtual void clear();

  /** Init the system PDE structures */
  virtual void init();
   
   /**
   * The \p NonlinearSolver defines the default interface used to
   * solve the nonlinear_implicit system.  This class handles all the
   * details of interfacing with various nonlinear algebra packages
   * like PETSc or LASPACK. Up to now also for the nonlinear case we use linear_solvers, in future we will add the nonlinear solver
   */
   vector<LinearSolver*> _LinSolver;
   
  /** Returns  the number of iterations taken for the most recent nonlinear solve. */
  unsigned int n_nonlinear_iterations() const { return _n_nonlinear_iterations; }

  /** Returns the final residual for the nonlinear system solve. */
  double final_nonlinear_residual() const { return _final_nonlinear_residual; }

protected:
  
   /**
   * The number of nonlinear iterations required to solve the nonlinear
   * system R(x)=0.
   */
  unsigned int _n_nonlinear_iterations;

  /**
   * The final residual for the nonlinear system R(x)
   */
  double _final_nonlinear_residual;

private:

  void CreateSystemPDEStructure();
  
};

#endif