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

#ifndef __ParsedFunction_gghpp__
#define __ParsedFunction_gghpp__

#include "FEMTTUConfig.h"

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FunctionBase.hpp"

#ifdef HAVE_FPARSER

#include "fparser.hh"

namespace femus {

class ParsedFunction : public FunctionBase {

public:

    /** Constructor */
    ParsedFunction();

    ParsedFunction(std::string expression, std::string independent_variables);

    /** Destructor */
    ~ParsedFunction() {};

    /** set the epression of the mathematical function */
    void SetExpression(const std::string expression);

    /** Set the independent_variables of the mathematical function */
    void SetIndependentVariables(const std::string independent_variables);

    /** Parse the expression */
    void Parse();

    virtual double operator() (double* x)
    {
        return _pfunc.Eval(x);
    }

private:

    FunctionParserBase<double> _pfunc;
    std::string _independent_variables;
    std::string _expression;


};

}

#else

namespace femus {

class ParsedFunction : public FunctionBase

{

public:

    /** constructors */
    ParsedFunction() : _dummy(0) {};

    
    ParsedFunction(std::string expression, std::string independent_variables) {};

    /** Destructor */
    ~ParsedFunction() {};

    /** () operator */
    virtual double operator() (double* x)
    {
        return _dummy;
    }

private:
    double _dummy;

};

}

#endif

#endif
