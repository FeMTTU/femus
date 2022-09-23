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
     * Helps in identifying
     * the system type in an equation system file.
    */
    virtual std::string system_type () const {
        return "NonlinearImplicitWithPrimalDualActiveSetMethod";
    }
    
    /** Set the active set flag name */
    void SetActiveSetFlagName(const std::vector<std::string> & name_in ) {
        _active_flag_name = name_in;
    }
    
    /** Set the active set flag name */
    std::vector<std::string> GetActiveSetFlagName() const {
        return _active_flag_name;
    }
    
    /** Solves the system. */
    virtual void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);
    
    
    void nonlinear_solve_single_level(const MgSmootherType& mgSmootherType, double & totalAssembyTime, const unsigned int grid0, const unsigned int igridn);

protected:
  
    std::vector< std::string > _active_flag_name;

    
// debug function sub-class - BEGIN ======   
protected:
    /** Debug function typedef */
    typedef void (*DebugFunc) (const MultiLevelProblem& ml_prob);  ///@deprecated
    
    /** Debug function typedef */
    typedef void (*DebugFuncLevel) (const MultiLevelProblem& ml_prob, 
                                    const unsigned level,
                                    const unsigned iteration,
                                    const std::vector<std::string> state_vars,  
                                    const std::vector<std::string> ctrl_vars  
    );
    
    /** @deprecated Debug function pointer */
    DebugFunc _debug_function; 

public:
    
    void set_state_vars(std::vector<std::string> state_vars_in ) {  _state_vars = state_vars_in;  } 

    void set_ctrl_vars(std::vector<std::string> ctrl_vars_in ) {  _ctrl_vars = ctrl_vars_in;  } 
    
    std::vector<std::string> get_state_vars() const {  return _state_vars;  } 

    std::vector<std::string> get_ctrl_vars() const { return _ctrl_vars;  } 
    
protected:
    /** Debug function pointer */
    DebugFuncLevel _debug_function_level;
    
    /**  */
    bool _debug_function_is_initialized;

    
    std::vector<std::string> _state_vars;  
    
    std::vector<std::string> _ctrl_vars;  
    
public:

    /** @deprecated Set the max number of non-linear iterations for the nonlinear system solve. */
    void SetDebugFunction(DebugFunc debug_func_in) { _debug_function = debug_func_in; 
                                                     _debug_function_is_initialized = true; 
    }
    
    /** Set the max number of non-linear iterations for the nonlinear system solve. */
    void SetDebugFunctionLevel(DebugFuncLevel debug_func_in) { _debug_function_level = debug_func_in; 
                                                     _debug_function_is_initialized = true; 
    }
    
    void print_iteration_and_do_additional_computations_with_given_function(const unsigned nonLinearIterator) const;
    
    void print_iteration_and_do_additional_computations_with_given_function_level(const unsigned nonLinearIterator, const unsigned level, 
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars  ) const;
    
    void do_additional_computations_with_given_function() const;
   
    void do_additional_computations_with_given_function_level(const unsigned level,  const unsigned nonLinearIterator, 
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars  ) const;
   
// debug function sub-class - END ======   
    
    
    
};



} //end namespace femus



#endif
