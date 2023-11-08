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

//==== Constructors / Destructor - BEGIN ========
public:

    /** Constructor.  Optionally initializes required data structures. */
    NonLinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number, const LinearEquationSolverType & smoother_type );

    /** destructor */
    virtual ~NonLinearImplicitSystem();
//==== Constructors / Destructor - END ========
    
//==== Basic - BEGIN ========
public:

    /** The type of the parent. */
    typedef LinearImplicitSystem Parent;

    /**
     * @returns \p "NonlinearImplicit".  Helps in identifying
     * the system type in an equation system file.
    */
    virtual std::string system_type () const {
        return "NonlinearImplicit";
    }
    
    /** Init the system PDE structures */
    virtual void init();

//==== Basic - END ========

    
//==== Matrix, Assemble - BEGIN ========
public:

    /** Only call assemble function */
    virtual void assemble_call_before_boundary_conditions(const unsigned int n_times);
    

//==== Matrix, Assemble - END ========


    
//==== Solver - BEGIN ========
public:
    
    /** Solves the system. */
    virtual void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);
    
    
//==== Solver - END ========    


//==== Solver, Convergence rate - BEGIN ========
public:

    void compute_convergence_rate() const;
    
protected:

    /** Vector of all nonlinear iterations for convergence rate */
    std::vector< NumericVector* >  _eps_fine;
    
//==== Solver, Convergence rate - END ========
    
public:

    /** Returns the final residual for the nonlinear system solve. */
    double final_nonlinear_residual() const {
        return _final_nonlinear_residual;
    }

    /** Returns the final residual for the nonlinear system solve. */
    const unsigned GetNonlinearIt() const { return _nonliniteration; }
    
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
    
    void print_iteration_to_file(const unsigned nonLinearIterator) const;
    
protected:

   clock_t total_mg_time_begin() const;
   
   double total_mg_time_end(const clock_t start_mg_time) const;
  
   clock_t nonlinear_time_begin() const;
   
   void nonlinear_time_end(const clock_t start_nl_time) const;
        
   void compute_assembly_vs_net_solver_times(const double totalSolverTime, const double totalAssemblyTime);

    /** The final residual for the nonlinear system R(x) */
    double _final_nonlinear_residual;

    /** The max number of non-linear iterations */
    unsigned int _n_max_nonlinear_iterations;

    /** The max non linear tolerance **/
    double _max_nonlinear_convergence_tolerance;

    unsigned _maxNumberOfResidualUpdateIterations;
    
    /** Flag for printing fields at each nonlinear iteration */
    bool _debug_nonlinear;
    
    /** Current nonlinear iteration index */
    unsigned _last_nonliniteration;
    
    /** Current nonlinear iteration index */
    unsigned _nonliniteration;

};



} //end namespace femus



#endif
