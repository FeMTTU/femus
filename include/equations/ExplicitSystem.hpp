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


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------


class ExplicitSystem : public System {

public:
  
  /** Constructor.  Optionally initializes required data structures. */
  ExplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number,const MgSmoother & smoother_type);

  virtual ~ExplicitSystem();
  
  /** The type of the parent. */
  typedef System Parent;
  
  /** Solves the system. */
  virtual void solve () {};
  
   /** Clear all the data structures associated with the system. */
  virtual void clear();

  /** Init the system PDE structures */
  virtual void init();
  
  /**
   * @returns \p "Explicit".  Helps in identifying
   * the system type in an equation system file.
  */
  virtual std::string system_type () const { return "Explicit"; }
  
protected:

private:


};


} //end namespace femus



#endif