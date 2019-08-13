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

protected:
    /** Debug function typedef */
    typedef void (*DebugFunc) (const MultiLevelProblem& ml_prob);
    
public:

    /** Constructor.  Optionally initializes required data structures. */
    NonLinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number, const LinearEquationSolverType & smoother_type );

    /** destructor */
    virtual ~NonLinearImplicitSystem();

    /** The type of the parent. */
    typedef LinearImplicitSystem Parent;

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

    /** Returns the final residual for the nonlinear system solve. */
    const unsigned GetNonlinearIt() const { return _nonliniteration; }
    
    /** Set the max number of non-linear iterations for the nonlinear system solve. */
    void SetDebugFunction(DebugFunc debug_func_in) { _debug_function = debug_func_in; 
                                                     _debug_function_is_initialized = true; 
    }
    
    /** Flag to print fields to file after each nonlinear iteration */
    void SetDebugNonlinear(const bool my_value);
    
    /** Set the max number of non-linear iterations for the nonlinear system solve. */
    void SetMaxNumberOfNonLinearIterations(unsigned int max_nonlin_it) {
        _n_max_nonlinear_iterations = max_nonlin_it;
    };

    /** Set the max nonlinear convergence tolerance */
    void SetNonLinearConvergenceTolerance(double nonlin_convergence_tolerance) {
        _max_nonlinear_convergence_tolerance = nonlin_convergence_tolerance;
    };

    /** Checks for the non the linear convergence */
    bool HasNonLinearConverged(const unsigned gridn, double &nonLinearEps);

    void SetMaxNumberOfResidualUpdatesForNonlinearIteration( const unsigned & maxNumberOfIterations){
      _n_max_linear_iterations = 1;
      _maxNumberOfResidualUpdateIterations = maxNumberOfIterations;
    }
    void SetResidualUpdateConvergenceTolerance(const double & tolerance){
      _linearAbsoluteConvergenceTolerance = tolerance;
    }
    
    void compute_convergence_rate() const;
    
    /** Solves the system. */
    virtual void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);
    
protected:

    /** The final residual for the nonlinear system R(x) */
    double _final_nonlinear_residual;

    /** The max number of non-linear iterations */
    unsigned int _n_max_nonlinear_iterations;

    /** The max non linear tolerance **/
    double _max_nonlinear_convergence_tolerance;

    unsigned _maxNumberOfResidualUpdateIterations;
    
    /** Vector of all nonlinear iterations for convergence rate */
    std::vector< NumericVector* >  _eps_fine;
    
    /** Flag for printing fields at each nonlinear iteration */
    bool _debug_nonlinear;
    
    /** Debug function pointer */
    DebugFunc _debug_function;
    
    /**  */
    bool _debug_function_is_initialized;

    /** Current nonlinear iteration index */
    unsigned _last_nonliniteration;
    
    /** Current nonlinear iteration index */
    unsigned _nonliniteration;

private:

    /** To be Added */
    void CreateSystemPDEStructure();

};



} //end namespace femus



#endif
