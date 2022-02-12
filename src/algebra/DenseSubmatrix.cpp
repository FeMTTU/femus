/*=========================================================================

 Program: FEMUS


} //end namespace femus


 Module: DenseSubMatrix
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "DenseSubmatrix.hpp"


namespace femus {




// =================================================
void DenseSubMatrix::left_multiply (const DenseMatrixBase& M2) {
  // (*this) <- M2 * M3 Where: (*this) = (m x n), M2= (m x p), M3= (p x n)
  // M3 is a simply a copy of *this
  DenseSubMatrix M3(*this);
  // Call the multiply function in the base class
  this->multiply(*this, M2, M3);
}

// =============================================================
void DenseSubMatrix::right_multiply (const DenseMatrixBase& M3) {
  // (*this) <- M2 * M3 Where: (*this) = (m x n), M2=(m x p), M3=(p x n)
  // M2 is simply a copy of *this
  DenseSubMatrix M2(*this);
  // Call the multiply function in the base class
  this->multiply(*this, M2, M3);
}


} //end namespace femus


