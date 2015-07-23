/*=========================================================================

 Program: FEMuS
 Module: NonLinearImplicitSystem
 Authors: Simone Bnà

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_NonLinearImplicitSystem_hpp__
#define __femus_equations_NonLinearImplicitSystem_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "LinearImplicitSystem.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

/**
 * The non linear implicit system abstract class
 */

class NonLinearImplicitSystem : public LinearImplicitSystem {

public:

    /** Constructor.  Optionally initializes required data structures. */
    NonLinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number, const MgSmoother & smoother_type );

    /** destructor */
    virtual ~NonLinearImplicitSystem();

    /** The type of the parent. */
    typedef LinearImplicitSystem Parent;

    /** Solves the system. */
    virtual void solve ();
    virtual void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);
    /** Clear all the data structures associated with the system. */
    virtual void clear();

    /** Init the system PDE structures */
    virtual void init();

    /**
     * @returns \p "NonlinearImplicit".  Helps in identifying
     * the system type in an equation system file.
    */
    virtual std::string system_type () const {
        return "NonlinearImplicit";
    }

    /** Returns the final residual for the nonlinear system solve. */
    double final_nonlinear_residual() const {
        return _final_nonlinear_residual;
    }

    /** Set the max number of non-linear iterations for the nonlinear system solve. */
    void SetMaxNumberOfNonLinearIterations(unsigned int max_nonlin_it) {
        _n_max_nonlinear_iterations = max_nonlin_it;
    };

    /** Set the max nonlinear convergence tolerance */
    void SetNonLinearConvergenceTolerance(double nonlin_convergence_tolerance) {
        _max_nonlinear_convergence_tolerance = nonlin_convergence_tolerance;
    };

    /** Checks for the non the linear convergence */
    bool IsNonLinearConverged(const unsigned gridn);

protected:

    /** The final residual for the nonlinear system R(x) */
    double _final_nonlinear_residual;

    /** The max number of non-linear iterations */
    unsigned int _n_max_nonlinear_iterations;

    /** The max non linear tolerance **/
    double _max_nonlinear_convergence_tolerance;


private:

    /** To be Added */
    void CreateSystemPDEStructure();

};



} //end namespace femus



#endif
