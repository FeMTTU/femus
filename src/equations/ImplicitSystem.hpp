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

#ifndef __femus_equations_ImplicitSystem_hpp__
#define __femus_equations_ImplicitSystem_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ExplicitSystem.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

/**
* The linear implicit system class
*/

class ImplicitSystem : public ExplicitSystem {

public:

    /** Constructor.  Optionally initializes required data structures. */
    ImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number,const LinearEquationSolverType & smoother_type);

    /** Destructor */
    virtual ~ImplicitSystem();

    /** The type of the parent. */
    typedef ExplicitSystem Parent;

    /** Init the system PDE structures */
    virtual void init();

    /** @deprecated Init the system PDE structures */
    virtual void init_two() {};

    /**
     * @returns \p "Implicit".  Helps in identifying
     * the system type in an equation system file.
     */
    virtual std::string system_type () const {
        return "Implicit";
    }

    // the sparse matrix must be putted here A, now is in linsysPDE

    /** Set a parameter option for the SparseMatrix A */
    virtual void SetMatrixOption(MatOption op, bool flag) {};

protected:


private:


};


} //end namespace femus



#endif
