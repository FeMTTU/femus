/*=========================================================================

Program: FEMuS
Module: ParsedFunction
Authors: Simone Bnà

Copyright (c) FEMuS
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "FemusConfig.hpp"

#ifdef HAVE_FPARSER

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ParsedFunction.hpp"

//C++ include
#include <iostream>
#include <cstdlib>

namespace femus {
  
 ParsedFunction::ParsedFunction() {
   _independent_variables = "x"; 
   _expression = "0.";
   _pfunc.Parse(_expression,_independent_variables);
   _pfunc.Optimize();
 }
  
  
 ParsedFunction::ParsedFunction(const std::string expression , const std::string independent_variables) : FunctionBase() {
  
    _independent_variables = independent_variables;
    _expression = expression;
  
    _pfunc.AddConstant("pi", std::acos(-1.));
    _pfunc.AddConstant("e", std::exp(1.));

    // Parse (and optimize if possible) the subexpression.
    // Add some basic constants, to Real precision.
    int res = _pfunc.Parse(_expression, _independent_variables);
    if(res >= 0) {
      std::cout << std::string(res+7, ' ') << "^\n"
                << _pfunc.ErrorMsg() << "\n\n";
      exit(1);
    }
    
    _pfunc.Optimize();

 }
 
 void ParsedFunction::SetExpression(const std::string expression) {
    _expression = expression; 
 }
 
  void ParsedFunction::SetIndependentVariables(const std::string independent_variables) {
    _independent_variables = independent_variables;
 }
 
 void ParsedFunction::Parse() {
   
    _pfunc.AddConstant("pi", std::acos(-1.));
    _pfunc.AddConstant("e", std::exp(1.));

    // Parse (and optimize if possible) the subexpression.
    // Add some basic constants, to Real precision.
    int res = _pfunc.Parse(_expression, _independent_variables);
    if(res >= 0) {
      std::cout << std::string(res+7, ' ') << "^\n"
                << _pfunc.ErrorMsg() << "\n\n";
      exit(1);
    }
    
    _pfunc.Optimize();
   
 }
 
 
}

#endif