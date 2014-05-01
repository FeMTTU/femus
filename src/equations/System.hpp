/*=========================================================================

 Program: FEMUS
 Module: System
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __System_hpp__
#define __System_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelProblem.hpp"


namespace femus {




//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class System;
class MultiLevelProblem;
class String;


class System {

public:

  /** Constructor.  Optionally initializes required data structures. */
  System (MultiLevelProblem& ml_prob, const std::string& name, const unsigned int number, const MgSmoother & smoother_type);

  /** destructor */
  virtual ~System();
  
  unsigned int number() const;
  
  const std::string & name() const;
  
  /**
   * @returns the type of system, helpful in identifying
   * which system type to use when reading equation system
   * data from file.  Should be overloaded in derived classes.
  */
  virtual std::string system_type () const { return "Basic"; }
  
  /** Associate the solution variables to the system PDE */
  void AddSolutionToSytemPDE(const char solname[]);
  
  /** Register a user function to use in assembling the system matrix and RHS. */
  void AttachAssembleFunction (void fptr(MultiLevelProblem &ml_prob, unsigned level, 
				      const unsigned &gridn, const bool &assembe_matrix));

  /** Solves the system.  Should be overloaded in derived systems. */
  virtual void solve () {};
  
   /** Clear all the data structures associated with the system. */
  virtual void clear();
  
  /** Init the system PDE structures */
  virtual void init();
  
  /** Get the index of the Solution "solname" for this system */
  unsigned GetSolPdeIndex(const char solname[]);
  

protected:
  
  /** Constant reference to the \p EquationSystems object used for the simulation. */
  MultiLevelProblem& _equation_systems;
  
  vector<mesh*> _msh;
  
  vector<Solution*> _solution;
  
  MultiLevelSolution* _ml_sol;
    
  vector <unsigned> _SolSystemPdeIndex;
 
  unsigned _gridn;
  
  unsigned _gridr;
  
  /** Function that assembles the system. */
  void (* _assemble_system_function) (MultiLevelProblem &ml_prob, unsigned level, 
				      const unsigned &gridn, const bool &assembe_matrix);
  
  /** The number associated with this system */
  const unsigned int _sys_number;
  
  /** A name associated with this system. */
  const std::string _sys_name;
  
};

// System inline methods
inline
const std::string & System::name() const
{
  return _sys_name;
}

inline
unsigned int System::number() const
{
  return _sys_number;
}


} //end namespace femus



#endif