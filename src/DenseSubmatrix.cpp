// $Id: dense_submatrix.C 3391 2009-05-26 03:50:35Z benkirk $

// Local Includes
#include "DenseSubmatrix.hpp"

// // ------------------------------------------------------------
// // Dense Matrix member functions

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

