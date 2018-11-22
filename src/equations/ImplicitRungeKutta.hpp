/*=========================================================================

 Program: FEMuS
 Module: ImplicitRungeKutta
 Authors: Kara Erdi

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_ImplicitRungeKutta_hpp__
#define __femus_equations_ImplicitRungeKutta_hpp__

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------
#include "TransientSystem.hpp"
#include <string>
#include <vector>


namespace femus {

/**
 * This class provides a specific system class for the time integration of system PDE
 * using the Newmark algorithm.
 */

// ------------------------------------------------------------
// TransientSystem class definition
template <class Base>
class ImplicitRungeKutta : public TransientSystem<Base> {

public:

    /** Constructor.  Initializes required data structures. */
    ImplicitRungeKutta (MultiLevelProblem& ml_probl,
                            const std::string& name,
                            const unsigned int number);

    /** Destructor. */
    virtual ~ImplicitRungeKutta ();

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

// -----------------------------------------------------------
// Useful typedefs
typedef ImplicitRungeKutta<LinearImplicitSystem> ImplicitRungeKuttaImplicitSystem;
typedef ImplicitRungeKutta<LinearImplicitSystem> ImplicitRungeKuttaLinearImplicitSystem;
typedef ImplicitRungeKutta<NonLinearImplicitSystem> ImplicitRungeKuttaNonlinearImplicitSystem;
typedef ImplicitRungeKutta<MonolithicFSINonLinearImplicitSystem> ImplicitRungeKuttaMonolithicFSINonlinearImplicitSystem;
typedef ImplicitRungeKutta<ExplicitSystem> ImplicitRungeKuttaExplicitSystem;
typedef ImplicitRungeKutta<System> ImplicitRungeKuttaBaseSystem;



} //end namespace femus



#endif
