/*=========================================================================

 Program: FEMUS
 Module: SlepcSVD
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_SlepcSVD_hpp__
#define __femus_algebra_SlepcSVD_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

#if defined(HAVE_PETSC) && defined(HAVE_SLEPC)

// C++ includes
#include <map>
#include <vector>
#include <cstdio>
// Local includes
#include "NumericVector.hpp"
#include "PetscMacro.hpp"
#include "ParalleltypeEnum.hpp"
#include "Casts.hpp"

/// Petsc include files.
EXTERN_C_FOR_PETSC_BEGIN
//#include <petscvec.h>
#include "slepc.h"
//#include "slepcsvd.h"
EXTERN_C_FOR_PETSC_END


namespace femus {

class SlepcSVD {  
  
public:
  double get_2norm_condition_number();
  SlepcSVD();
  ~SlepcSVD();

private:
  SVD m_svd;
  
};


} //end namespace femus

#endif
#endif