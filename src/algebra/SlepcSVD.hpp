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

// Local includes
#include "PetscMacro.hpp"
#include "SparseMatrix.hpp"

/// Petsc include files.
EXTERN_C_FOR_PETSC_BEGIN
#include <petscmat.h>
#include "slepcsvd.h"
EXTERN_C_FOR_PETSC_END


namespace femus {

class SlepcSVD {  
  
public:
  
  SlepcSVD();
  
  ~SlepcSVD();
  
  void set_operator(Mat A);
  
  void set_operator(SparseMatrix* Amat);
  
  void init();

  double compute_2norm_condition_number();  

private:
  
  SVD m_svd;
  
};


} //end namespace femus

#endif
#endif