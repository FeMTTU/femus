/*=========================================================================

Program: FEMuS
Module: FunctionBase
Authors: Simone Bnà

Copyright (c) FEMuS
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_FunctionBase_hpp__
#define __femus_algebra_FunctionBase_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------


namespace femus {

class FunctionBase {

public:

    /** Constructor */
    FunctionBase();

    /** Destructor */
    ~FunctionBase();
    
    virtual double operator() (double* x) = 0;
    
};

}

#endif