/*=========================================================================

 Program: FEMUS
 Module: ImplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __implicit_system_h_
#define __implicit_system_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ExplicitSystem.hpp"


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------


class ImplicitSystem : public ExplicitSystem { 

public:

  /** Constructor.  Optionally initializes required data structures. */
  ImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number,const MgSmoother & smoother_type);

  virtual ~ImplicitSystem();
  
  /** The type of the parent. */
  typedef ExplicitSystem Parent;
  
  /** Solves the system. */
  virtual void solve () {};
  
  /** Clear all the data structures associated with the system. */
  virtual void clear();

  /** Init the system PDE structures */
  virtual void init();
  
  /**
   * @returns \p "Implicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "Implicit"; }

  // the sparse matrix must be putted here A, now is in linsysPDE
  
  /** Set a parameter option for the SparseMatrix A */
  virtual void SetMatrixOption(MatOption op, bool flag) {};
  
protected:
  
 
private:

 
};


} //end namespace femus



#endif