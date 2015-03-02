/*=========================================================================

 Program: FEMUS
 Module: DenseSubVector
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_DenseSubvector_hpp__
#define __femus_algebra_DenseSubvector_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "DenseVector.hpp"


namespace femus {



// ===========================================
// Defines a dense subvector for use in Finite Element-type computations.
// Useful for storing element load vectors  before summation
// into a global vector, particularly when you have systems of equations.
//=====================================

// ------------------------------------------------------------
// DenseSubVector class definition
class DenseSubVector : public DenseVectorBase {
  // =====================================
  // DATA
  // =====================================
private:

  /// The parent vector that contains this subvector.
  DenseVector& _parent_vector;
  /// The length of this subvector.
  unsigned int _n;
  /// The offset into the parent vector.
  unsigned int _i_off;

public:
  // =====================================
  // CONSTRUCTOR /DESTRUCTOR
  // =====================================
  /// Constructor.  Creates a dense subvector of the vector
  DenseSubVector(DenseVector& parent,const unsigned int ioff=0, const unsigned int n=0);

  /// Destructor.  Does nothing.
  virtual ~DenseSubVector() {}

  /// @returns a reference to the parent vector.
  DenseVector& parent () {
    return _parent_vector;
  }

  /// Set every element in the subvector to 0.
  virtual void zero();

  // ====================================
  // RETURN FUNCTIONS
  // ====================================

  /// @returns the \p (i,j) element of the subvector.
  double operator() (const unsigned int i) const;

  /// @returns the \p (i,j) element of the subvector as a writeable reference.
  double & operator() (const unsigned int i);

  /// @returns the \p (i) element of the vector.
  virtual double el(const unsigned int i) const {
    return (*this)(i);
  }

  /// @returns the \p (i) element of the vector as a writeable reference.
  virtual double & el(const unsigned int i)     {
    return (*this)(i);
  }

  ///  @returns the size of the subvector.
  virtual unsigned int size() const {
    return _n;
  }

  ///  @returns the row offset into the parent vector.
  unsigned int i_off() const {
    return _i_off;
  }

  /// Changes the location of the subvector in the parent vector.
  void reposition(const unsigned int ioff,const unsigned int n);


};



// ------------------------------------------------------------
// Dense Vector member functions
inline
DenseSubVector::DenseSubVector(DenseVector& parent,
                                 const unsigned int ioff,
                                 const unsigned int n) :
  _parent_vector(parent) {
  reposition (ioff, n);
}



inline void DenseSubVector::reposition(const unsigned int ioff,
                                        const unsigned int n) {
  _i_off = ioff;
  _n = n;
  // Make sure we still fit in the parent vector.
  assert ((this->i_off() + this->size()) <= _parent_vector.size());
}

inline void DenseSubVector::zero() {
  for (unsigned int i=0; i<this->size(); i++)    _parent_vector (i + this->i_off()) = 0.;
}

inline double DenseSubVector::operator () (const unsigned int i) const {
  assert (i < this->size());
  assert (i + this->i_off() < _parent_vector.size());
  return _parent_vector (i + this->i_off());
}

// =========================================
inline double & DenseSubVector::operator () (const unsigned int i) {
  assert (i < this->size());
  assert (i + this->i_off() < _parent_vector.size());

  return _parent_vector (i + this->i_off());
}




} //end namespace femus



#endif // #ifndef __dense_vector_h__

