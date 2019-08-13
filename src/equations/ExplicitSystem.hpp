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

#ifndef __femus_equations_ExplicitSystem_hpp__
#define __femus_equations_ExplicitSystem_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "System.hpp"

namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------


/**
 * The explicit system class
 */

class ExplicitSystem : public System {

public:

    /** Constructor.  Optionally initializes required data structures. */
    ExplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number,const LinearEquationSolverType & smoother_type);

    /** Destructor */
    virtual ~ExplicitSystem();

    /** The type of the parent. */
    typedef System Parent;

    /** Init the system PDE structures */
    virtual void init();

    /** @deprecated Init the system PDE structures */
    virtual void init_two(){};

    /**
     * @returns \p "Explicit".  Helps in identifying
     * the system type in an equation system file.
    */
    virtual std::string system_type () const {
        return "Explicit";
    }

protected:

private:


};


} //end namespace femus



#endif
