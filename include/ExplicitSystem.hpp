/*=========================================================================

 Program: FEMUS
 Module: ExplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __explicit_system_h_
#define __explicit_system_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "System.hpp"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------


class ExplicitSystem : public System {

public:
  
  /** Constructor.  Optionally initializes required data structures. */
  ExplicitSystem (NonLinearMultiLevelProblem& es, const std::string& name, const unsigned int number);

  virtual ~ExplicitSystem();
  
  /** Solves the system. */
  virtual void solve () {};
  
   /** Clear all the data structures associated with the system. */
  virtual void clear();

  /** Init the system PDE structures */
  virtual void init();
  
protected:

private:


};

#endif