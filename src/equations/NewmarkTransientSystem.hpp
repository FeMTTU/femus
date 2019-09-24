/*=========================================================================

 Program: FEMuS
 Module: NewmarkTransientSystem
 Authors: Simone Bnà

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_NewmarkTransientSystem_hpp__
#define __femus_equations_NewmarkTransientSystem_hpp__

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------
#include "TransientSystem.hpp"
#include <string>
#include <vector>
#include "assert.h"

namespace femus {
    
/**
 * This class provides a specific system class for the time integration of system PDE
 * using the Newmark algorithm.
 */

// ------------------------------------------------------------
// TransientSystem class definition
template <class Base>
class NewmarkTransientSystem : public TransientSystem<Base> {

public:

    /** Constructor.  Initializes required data structures. */
    NewmarkTransientSystem (MultiLevelProblem& ml_probl,
                            const std::string& name,
                            const unsigned int number,
                            const LinearEquationSolverType & smoother_type);

    /** Destructor. */
    virtual ~NewmarkTransientSystem ();

    /** Update the acceleration using the Newmark algorithm */
    void UpdateAcceleration(const std::vector<std::string>& vel_vars, const std::vector<std::string>& acc_vars);

    /** Set the Newmark parameters */
    void SetNewmarkParameters(const double gamma, const double delta);

private:

    // member data
    double _gamma;
    double _delta;
    double _a5;
    double _a1;
    double _a2;

};


template <class Base>
NewmarkTransientSystem<Base>::NewmarkTransientSystem(
            MultiLevelProblem& ml_probl,
            const std::string& name,
            const unsigned int number, 
            const LinearEquationSolverType & smoother_type):
            TransientSystem<Base>(ml_probl, name, number, smoother_type),
		   _gamma(0.5),
		   _delta(0.5)
{
  _a5 = -1.*(1. - _gamma)/_gamma;
  _a1 = 1./(_gamma*this->_dt);
  _a2 = -1./(_gamma*this->_dt);
}

/** Destructor. */
template <class Base>
NewmarkTransientSystem<Base>::~NewmarkTransientSystem()
{
  
}

template <class Base>
void NewmarkTransientSystem<Base>::UpdateAcceleration(const std::vector<std::string>& vel_vars, const std::vector<std::string>& acc_vars)
{
//   const double gamma = 0.5;
//   const double a5    = -1.*(1. - gamma)/gamma;
//   const double a1    = 1./(gamma*this->_dt);
//   const double a2    = -1./(gamma*this->_dt);
  
//  unsigned dim = this->_msh[0]->GetDimension();
  // vel and acc must have the same number of variables
  assert(vel_vars.size() == acc_vars.size());
  
  const unsigned dim = vel_vars.size();
 
  unsigned axyz[3];
  unsigned vxyz[3];
//   const char accname[3][3] = {"AX","AY","AZ"};
//   const char velname[3][2] = {"U","V","W"};
   
  for(unsigned i=0; i<dim; i++) {
     axyz[i] = this->_equation_systems.GetIndex(acc_vars.at(i));
     vxyz[i] = this->_equation_systems.GetIndex(vel_vars.at(i));
  }
   
  for (int ig=0;ig< this->_gridn;ig++) {
    for(unsigned i=0; i<dim; i++) {
      this->_solution[ig]->_Sol[axyz[i]]->scale(_a5);
      this->_solution[ig]->_Sol[axyz[i]]->add(_a1,*(this->_solution[ig]->_Sol[vxyz[i]]));
      this->_solution[ig]->_Sol[axyz[i]]->add(_a2,*(this->_solution[ig]->_SolOld[vxyz[i]]));
    }
  }
}

template <class Base>
void NewmarkTransientSystem<Base>::SetNewmarkParameters(const double gamma, const double delta)
{
  _gamma = gamma;
  _delta = delta;
  
  _a5 = -1.*(1. - _gamma)/_gamma;
  _a1 = 1./(_gamma*this->_dt);
  _a2 = -1./(_gamma*this->_dt);
  
}

// -----------------------------------------------------------
// Useful typedefs
typedef NewmarkTransientSystem<LinearImplicitSystem> NewmarkTransientImplicitSystem;
typedef NewmarkTransientSystem<LinearImplicitSystem> NewmarkTransientLinearImplicitSystem;
typedef NewmarkTransientSystem<NonLinearImplicitSystem> NewmarkTransientNonlinearImplicitSystem;
typedef NewmarkTransientSystem<MonolithicFSINonLinearImplicitSystem> NewmarkTransientMonolithicFSINonlinearImplicitSystem;
typedef NewmarkTransientSystem<ExplicitSystem> NewmarkTransientExplicitSystem;
typedef NewmarkTransientSystem<System> NewmarkTransientBaseSystem;

} //end namespace femus



#endif
