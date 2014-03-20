/*=========================================================================

 Program: FEMUS
 Module: LinearImplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __linear_implicit_system_h_
#define __linear_implicit_system_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ImplicitSystem.hpp"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

class LinearImplicitSystem : public ImplicitSystem {

public:

/** Constructor.  Optionally initializes required data structures. */
  LinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number);

  virtual ~LinearImplicitSystem();
  
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
 

protected:
  
  /** The number of linear iterations required to solve the linear system Ax=b. */
  unsigned int _n_linear_iterations;

  /** The final residual for the linear system Ax=b. */
  double _final_linear_residual;

  void Prolongator(const unsigned &gridf);
  
  void Restrictor(const unsigned &gridf, const unsigned &gridn, 
					    const unsigned &non_linear_iteration, const unsigned &linear_iteration, const bool &full_cycle);
  
  int BuildProlongatorMatrix(unsigned gridf, const char pdename[]);
  
  
private:

  void CreateSystemPDEStructure();

};

#endif