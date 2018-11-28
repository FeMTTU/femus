/*=========================================================================

 Program: FEMuS
 Module: ImplicitRungeKuttaSystem
 Authors: Simone Bnà

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_ImplicitRungeKuttaSystem_hpp__
#define __femus_equations_ImplicitRungeKuttaSystem_hpp__

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
class ImplicitRungeKuttaSystem : public TransientSystem<Base> { 

public:

    /** Constructor.  Initializes required data structures. */
    ImplicitRungeKuttaSystem (MultiLevelProblem& ml_probl,
                            const std::string& name,
                            const unsigned int number,
                            const MgSmoother & smoother_type);

    /** Destructor. */
    virtual ~ImplicitRungeKuttaSystem ();

    void AddSolutionToSystemPDE(const char solname[]);
    
    inline void SetRKStage(const unsigned & RK){
      _RK = RK;
    }
    
    inline unsigned GetRKStage(){
      return _RK;
    }
    
private:
   unsigned _RK;
};


template <class Base>
ImplicitRungeKuttaSystem<Base>::ImplicitRungeKuttaSystem(
            MultiLevelProblem& ml_probl,
            const std::string& name,
            const unsigned int number, 
            const MgSmoother & smoother_type):
            TransientSystem<Base>(ml_probl, name, number, smoother_type),
		    _RK(1.0)
{
 
}

/** Destructor. */
template <class Base>
ImplicitRungeKuttaSystem<Base>::~ImplicitRungeKuttaSystem()
{
  
}

template <class Base>
void ImplicitRungeKuttaSystem<Base>::AddSolutionToSystemPDE(const char solname[]){
    
  for(unsigned i = 0; i < _RK; i++){
    std::ostringstream solnameki;
    solnameki << solname << "k" << i+1;
    this->Base::AddSolutionToSystemPDE(solnameki.str().c_str());
  }
  
}

// -----------------------------------------------------------
// Useful typedefs
typedef ImplicitRungeKuttaSystem<LinearImplicitSystem> ImplicitRungeKuttaLinearImplicitSystem;
typedef ImplicitRungeKuttaSystem<NonLinearImplicitSystem> ImplicitRungeKuttaNonlinearImplicitSystem;
typedef ImplicitRungeKuttaSystem<MonolithicFSINonLinearImplicitSystem> ImplicitRungeKuttaMonolithicFSINonlinearImplicitSystem;

} //end namespace femus



#endif






