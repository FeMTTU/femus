/*=========================================================================

 Program: FEMUS
 Module: DenseSubMatrix
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_DenseSubmatrix_hpp__
#define __femus_algebra_DenseSubmatrix_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

// Local Includes
#include "DenseMatrixBase.hpp"
#include "DenseMatrix.hpp"
#include "DenseSubvector.hpp"

#include <cassert>

namespace femus {



/**
 * Defines a dense submatrix for use in Finite Element-type computations.
 * Useful for storing element stiffness matrices before summation
 * into a global matrix, particularly when you have systems of equations.
 */

// ------------------------------------------------------------
// DenseSubMatrix class definition

class DenseSubMatrix : public DenseMatrixBase {

private:

  /// The parent matrix that contains this submatrix.
  DenseMatrix& _parent_matrix;
  /// The row offset into the parent matrix.
  int _i_off;
  /// The column offset into the parent matrix.
  int _j_off;

public:

  /**
   * Constructor.  Creates a dense submatrix of the matrix
   * \p parent.  The submatrix has dimensions \f$(m \times n)\f$,
   * and the \f$(0,0) entry of the submatrix is located
   * at the \f$(ioff,joff)\f$ location in the parent matrix.
   */
  DenseSubMatrix(DenseMatrix& parent,
                  const int ioff=0,
                  const int joff=0,
                  const int m=0,
                  const int n=0) :
    DenseMatrixBase(m,n), _parent_matrix(parent) {
    this->reposition (ioff, joff, m, n);
  };

  /// Copy Constructor.
  DenseSubMatrix(const DenseSubMatrix& other_matrix)
    : DenseMatrixBase(other_matrix._m, other_matrix._n),
      _parent_matrix(other_matrix._parent_matrix) {
    _i_off = other_matrix._i_off;
    _j_off = other_matrix._j_off;
  }
  /// Destructor
  virtual ~DenseSubMatrix() {}

  /// Set every element in the submatrix to 0.
  virtual void zero() {
    for (int i=0; i<this->m(); i++)
      for (int j=0; j<this->n(); j++)
        _parent_matrix(i + this->i_off(),j + this->j_off()) = 0.;
  }
  // =====================================
  // Return functions
  // =====================================
  /// @returns a reference to the parent matrix.
  DenseMatrix& parent () {
    return _parent_matrix;
  }

  ///  @returns the \p (i,j) element of the submatrix
  double operator() (const int i,const int j) const;
  ///   @returns the \p (i,j) element of the submatrix as a writeable reference.
  double & operator() (const int i,const int j) {
    assert (i < this->m());
    assert (j < this->n());
    assert (i + this->i_off() < _parent_matrix.m());
    assert (j + this->j_off() < _parent_matrix.n());
    return _parent_matrix (i + this->i_off(),j + this->j_off());
  }


  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  virtual double el(const int i,const int j) const {
    return (*this)(i,j);
  }
  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  virtual double & el(const int i,const int j)     {
    return (*this)(i,j);
  }

  /// @returns the row offset into the parent matrix.
  int i_off() const {
    return _i_off;
  }
  /// @returns the column offset into the parent matrix.
  int j_off() const {
    return _j_off;
  }

  // ==================================================
  // Operations
  // ==================================================
  /// Performs the operation: (*this) <- M2 * (*this)
  virtual void left_multiply (const DenseMatrixBase& M2);

  /// Performs the operation: (*this) <- (*this) * M3
  virtual void right_multiply (const DenseMatrixBase& M3);

  /// Changes the location of the submatrix in the parent matrix.
  void reposition(const int ioff,  const int joff,
                  const int m,  const int n);

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const int i,const int j,
                const double val,DenseSubVector& rhs) {
    this->parent().condense(this->i_off()+i,
                            this->j_off()+j,
                            val, rhs.parent());
  }

};


// --------------------------------------------------
// Constructor








inline
double DenseSubMatrix::operator () (const int i,
                                     const int j) const {
  assert (i < this->m());
  assert (j < this->n());
  assert (i + this->i_off() < _parent_matrix.m());
  assert (j + this->j_off() < _parent_matrix.n());

  return _parent_matrix (i + this->i_off(),
                         j + this->j_off());
}

inline
void DenseSubMatrix::reposition(const int ioff,
                                 const int joff,
                                 const int m,
                                 const int n) {
  _i_off = ioff;
  _j_off = joff;
  this->_m = m;
  this->_n = n;

  // Make sure we still fit in the parent matrix.
  assert ((this->i_off() + this->m()) <= _parent_matrix.m());
  assert ((this->j_off() + this->n()) <= _parent_matrix.n());
}





} //end namespace femus



#endif // #ifndef __dense_matrix_h__

