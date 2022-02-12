/*=========================================================================

 Program: FEMUS
 Module: DenseMatrixBase
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_DenseMatrixBase_hpp__
#define __femus_algebra_DenseMatrixBase_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <cassert>
#include <iostream>


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class DenseVectorBase;

/**
 * Defines an abstract dense matrix base class for use in Finite Element-type
 * computations.  Specialized dense matrices, for example DenseSubMatrices,
 * can be derived from this class.
 */
class DenseMatrixBase {

protected:

  // ===============================
  // Data
  // =============================
  int _m; ///< The row dimension.
  int _n;///< The column dimension.

  // ==================================
  // Construction/Destr
  // ===============================
  /// Constructor.  Creates a dense matrix of dimension \p m by \p n.
  /// Protected so that there is no way the user can create one.
  DenseMatrixBase(const  int m=0,const  int n=0) : _m(m), _n(n) {};

public:
  ///Destructor. Empty.
  virtual ~DenseMatrixBase() {};
/// Set every element in the matrix to 0.
  virtual void zero() = 0;

  // ===============================
  // Return functions
  // ===============================
  /// @returns the \p (i,j) element of the matrix.
  virtual double el(const  int i,const  int j) const = 0;
  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  virtual double & el(const  int i, const  int j) = 0;

  /// @returns the row-dimension of the matrix.
  int m() const {
    return _m;
  }
  /// @returns the column-dimension of the matrix.
  int n() const {
    return _n;
  }

  // ===============================
  // Algebra
  // ===============================
  /// Performs the operation: (*this) <- M2 * (*this)
  virtual void left_multiply (const DenseMatrixBase& M2) = 0;
  /// Performs the operation: (*this) <- (*this) * M3
  virtual void right_multiply (const DenseMatrixBase& M3) = 0;
  /// Adds \p a to every element; matrix += a mat
  void   add (const double a, const DenseMatrixBase& mat);

  // ===============================
  // print
  // ===============================
  /// Pretty-print the matrix to \p stdout.
  void print(std::ostream& os) const;
  /// Formatted print as above but allows you to do DenseMatrix K; std::cout << K << std::endl;
  friend std::ostream& operator << (std::ostream& os, const DenseMatrixBase& m) {
    m.print(os);
    return os;
  }
  /// Prints the matrix entries with more decimal places in scientific notation.
  void print_scientific(std::ostream& os) const;



protected:
  ///  Performs the computation M1 = M2 * M3 where: M1 = (m x n) M2 = (m x p) M3 = (p x n)
  void multiply (DenseMatrixBase& M1, const DenseMatrixBase& M2,const DenseMatrixBase& M3);

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const  int i,const  int j,
                const double val,DenseVectorBase& rhs);

};


// ===========================================================
//  INLINE FUNCTIONS
// ===========================================================
inline
void DenseMatrixBase::add (const double factor,const DenseMatrixBase& mat) {
  assert (this->m() == mat.m());
  assert (this->n() == mat.n());
  for ( int j=0; j<this->n(); j++)
    for ( int i=0; i<this->m(); i++)
      this->el(i,j) += factor*mat.el(i,j);
}



} //end namespace femus



#endif // #ifndef __dense_matrix_base_h__

