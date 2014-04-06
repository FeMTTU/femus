/*=========================================================================

 Program: FEMUS
 Module: NewmarkTransientSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __newmark_transient_system_h__
#define __newmark_transient_system_h__

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------
#include "TransientSystem.hpp"
#include "string"
#include "vector"

/**
 * This class provides a specific system class for the time integration of system PDE 
 * using the Newmark algorithm.
 */

// ------------------------------------------------------------
// TransientSystem class definition
template <class Base>
class NewmarkTransientSystem : public TransientSystem<Base>
{
public:

  /** Constructor.  Initializes required data structures.  */
  NewmarkTransientSystem (MultiLevelProblem& ml_probl,
		   const std::string& name,
		   const unsigned int number);

  /** Destructor. */
  virtual ~NewmarkTransientSystem ();
  
 
  /** Update the acceleration using the Newmark algorithm */
  void UpdateAcceleration(const std::vector<std::string>& vel_vars, const std::vector<std::string>& acc_vars);
  
  /** Set the Newmark parameters */
  void SetNewmarkParameters(const double gamma, const double delta);
  
private:

  double _gamma;
  double _delta;
  double _a5;
  double _a1;
  double _a2;  
  
};

// -----------------------------------------------------------
// Useful typedefs
 typedef NewmarkTransientSystem<LinearImplicitSystem> NewmarkTransientImplicitSystem;
 typedef NewmarkTransientSystem<LinearImplicitSystem> NewmarkTransientLinearImplicitSystem;
 typedef NewmarkTransientSystem<NonLinearImplicitSystem> NewmarkTransientNonlinearImplicitSystem;
 typedef NewmarkTransientSystem<MonolithicFSINonLinearImplicitSystem> NewmarkTransientMonolithicFSINonlinearImplicitSystem;
 typedef NewmarkTransientSystem<ExplicitSystem> NewmarkTransientExplicitSystem;
 typedef NewmarkTransientSystem<System> NewmarkTransientBaseSystem;


#endif