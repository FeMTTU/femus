/*=========================================================================

 Program: FEMUS
 Module: FemusInit
 Authors: Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_utils_FemusInit_hpp__
#define __femus_utils_FemusInit_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

//use this class only with petsc-mpi
#ifdef HAVE_MPI
#include <mpi.h>  //For MPI_COMM_WORLD
#endif

// PETSC ----------------------
#ifdef HAVE_PETSC
// EXTERN_C_FOR_PETSC_BEGIN
# include <petsc.h>
# include <petscerror.h>
// EXTERN_C_FOR_PETSC_END
#endif


// ========================================
/// Class FemusInit
// ========================================

#include "adept.h"

#include "uq.hpp"

namespace femus {



class FemusInit {
  
public:
  
// === Constructors / Destructor  - BEGIN =================
    /// Constructor
    FemusInit(int & argc,char** & argv, MPI_Comm comm_world_in = MPI_COMM_WORLD);

    /// Destructor
    ~FemusInit();
// === Constructors / Destructor  - END =================


// === Automatic differentiation (needed for automatic Jacobian computations) - BEGIN =================
    static adept::Stack _adeptStack;
// === Automatic differentiation (needed for automatic Jacobian computations) - END =================
    
    
// === UQ Quadratures - BEGIN =================
    static uq _uqHermite; 
    static uq _uqLegendre; 
// === UQ Quadratures - END =================
     
};


} //end namespace femus



#endif // end _femus_init ----------------------
