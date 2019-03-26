/*=========================================================================

 Program: FEMuS
 Module: NonLinearImplicitSystemWithPrimalDualActiveSetMethod
 Authors: Giorgio Bornia, Saikanth Ratnavale

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_NonLinearImplicitSystemWithPrimalDualActiveSetMethod_hpp__
#define __femus_equations_NonLinearImplicitSystemWithPrimalDualActiveSetMethod_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "NonLinearImplicitSystem.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

/**
 * The non linear implicit system abstract class
 */

class NonLinearImplicitSystemWithPrimalDualActiveSetMethod : public NonLinearImplicitSystem {

public:

    /** Constructor.  Optionally initializes required data structures. */
    NonLinearImplicitSystemWithPrimalDualActiveSetMethod (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number, const LinearEquationSolverType & smoother_type );

    /**
     * @returns \p "NonlinearImplicit".  Helps in identifying
     * the system type in an equation system file.
    */
    virtual std::string system_type () const {
        return "NonlinearImplicitWithPrimalDualActiveSetMethod";
    }
    
    /** Set the active set flag name */
    void SetActiveSetFlagName(const std::string & name_in ) {
        _active_flag_name = name_in;
    };
    
    /** Solves the system. */
    virtual void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);

protected:

        
    std::string _active_flag_name;

};



} //end namespace femus



#endif
