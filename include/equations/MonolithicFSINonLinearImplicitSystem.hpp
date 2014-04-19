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

#ifndef __MonolithicFSI_nonlinear_implicit_system_h_
#define __MonolithicFSI_nonlinear_implicit_system_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "NonLinearImplicitSystem.hpp"


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------


class MonolithicFSINonLinearImplicitSystem : public NonLinearImplicitSystem {

public:

/** Constructor.  Optionally initializes required data structures. */
  MonolithicFSINonLinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number);
  
  virtual ~MonolithicFSINonLinearImplicitSystem();
  
  /** The type of the parent. */
  typedef NonLinearImplicitSystem Parent;
  
  /** Clear all the data structures associated with the system. */
  virtual void clear();

  /** Init the system PDE structures */
  virtual void init();
   
  /**
   * @returns \p "MonolithicFSINonlinearImplicit".  Helps in identifying
   * the system type in an equation system file.
  */
  virtual std::string system_type () const { return "MonolithicFSINonlinearImplicit"; }
      
protected:
 
 /** Create the Restrictor matrix for the Multigrid solver */ 
 void Restrictor(const unsigned &gridf, const unsigned &gridn, const unsigned &non_linear_iteration, 
                 const unsigned &linear_iteration, const bool &full_cycle);
 
 /** Create the Prolongator and Restrictor Operators for the Algebraic Multigrid Solver */ 
 void BuildProlongatorMatrix(unsigned gridf);
  
private:
 
};



} //end namespace femus



#endif