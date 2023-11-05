/*=========================================================================

 Program: FEMUS
 Module: DenseVector
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_DenseVector_hpp__
#define __femus_algebra_DenseVector_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "DenseVectorBase.hpp"

// C++ includes
#include <vector>
#include <cassert>
#include <cmath>


namespace femus {



/**
 * Defines a dense vector for use in Finite Element-type computations.
 * This class is to basically compliment the \p DenseMatix class.  It
 * has additional capabilities over the \p std::vector that make it
 * useful for finite elements, particulary for sys of equations.
 */

// ------------------------------------------------------------
// DenseVector class definition
class DenseVector : public DenseVectorBase {
private:

  /// The actual data values, stored as a 1D array.
  std::vector<double> _val;

public:

  /// Constructor.  Creates a dense vector of dimension \p n.
  explicit
  DenseVector(const unsigned int n=0);
  /// Copy-constructor.
  DenseVector (const DenseVector& other_vector);

  /// Copy-constructor, from a \p std::vector.
  DenseVector (const std::vector<double>& other_vector);

  /// Destructor.  Does nothing.
  ~DenseVector() {}


  ///  * Resize the vector. Sets all elements to 0.
  void resize (const unsigned int n);
  /// Set every element in the vector to 0.
  virtual void zero();

  // ============================
  // RETURN FUNCTIONS
  // ===========================

  /// Puts the principal subvector of size \p sub_n
  void get_principal_subvector (unsigned int sub_n, DenseVector& dest) const;
  ///Access to the values array
  std::vector<double>& get_values() {
    return _val;
  }
  /// Access to the values array.
  const std::vector<double>& get_values() const {
    return _val;
  }

  ///  * @returns the size of the vector.
  virtual unsigned int size() const {
    return _val.size();
  }

  /// @returns the \p (i) element of the vector.
  double operator() (const unsigned int i) const;
  /// * @returns the \p (i,j) element of the vector as a writeable reference.
  double & operator() (const unsigned int i);

  ///  * @returns the \p (i) element of the vector.
  virtual double el(const unsigned int i) const {
    return (*this)(i);
  }
  ///  * @returns the \p (i) element of the vector as a writeable reference.
  virtual double & el(const unsigned int i)     {
    return (*this)(i);
  }

  // ============================
  // ALGEBRA
  // ===========================

  ///* Assignment operator.
  DenseVector& operator = (const DenseVector& other_vector);
  ///  * STL-like swap method
  void swap(DenseVector& other_vector);
  /// Multiplies every element in the vector by \p factor
  void scale (const double factor);

  /// * Multiplies every element in the vector by \p factor.
  DenseVector& operator*= (const double factor);

  ///  * Adds \p factor times \p vec to this vector. T += a * vec
  void  add (const double a,  const DenseVector& vec);

  ///* Evaluate dot product with \p vect
  double dot (const DenseVector &vec) const;

  /// Tests if \p vec is exactly equal to this vector.
  bool operator== (const DenseVector &vec) const;
  /// Tests if \p vec is not exactly equal to this vector.
  bool operator!= (const DenseVector &vec) const;

  /// Adds \p vec to this vector.
  DenseVector& operator+= (const DenseVector &vec);

  /// Subtracts \p vec from this vector.
  DenseVector& operator-= (const DenseVector &vec);

  /// * @returns the minimum element in the vector.
  double min () const;

  /// @returns the maximum element in the vector.
  double max () const;

  /// @returns the \f$l_1\f$-norm of the vector,
  double l1_norm () const;
  /// @returns the \f$l_2\f$-norm of the vector
  double l2_norm () const;
  /// @returns the maximum absolute value of the
  double linfty_norm () const;
};



// ------------------------------------------------------------
// DenseVector member functions

inline DenseVector::DenseVector(const unsigned int n) : _val (n, 0.) {}


// ======================================================
inline DenseVector::DenseVector (const DenseVector& other_vector) :
  DenseVectorBase() {
  const std::vector<double> &other_vals = other_vector.get_values();
  _val.clear();
  _val.reserve(other_vals.size());

  for (unsigned int i=0; i<other_vals.size(); i++)   _val.push_back(other_vals[i]);
}

// =================================
inline DenseVector::DenseVector (const std::vector<double>& other_vector) :
  _val(other_vector) {  }


// =============================================
inline DenseVector& DenseVector::operator = (const DenseVector& other_vector) {
  //  _val = other_vector._val;
  const std::vector<double> &other_vals = other_vector.get_values();
  _val.clear();
  _val.reserve(other_vals.size());
  for (unsigned int i=0; i<other_vals.size(); i++)    _val.push_back(other_vals[i]);
  return *this;
}

// ======================================
inline void DenseVector::swap(DenseVector& other_vector) {
  _val.swap(other_vector._val);
}
// =================================================
inline void DenseVector::resize(const unsigned int n) {
  _val.resize(n);
  zero();
}
// ==================
inline void DenseVector::zero() {
  std::fill (_val.begin(),_val.end(),0.);
}
// =================================
inline double DenseVector::operator () (const unsigned int i) const {
  assert (i < _val.size());
  return _val[i];
}
// ====================================
inline double & DenseVector::operator () (const unsigned int i) {
  assert (i < _val.size());
  return _val[i];
}
// ====================================
inline void DenseVector::scale (const double factor) {
  for (unsigned int i=0; i<_val.size(); i++)    _val[i] *= factor;
}
// =======================================
inline DenseVector& DenseVector::operator*= (const double factor) {
  this->scale(factor);
  return *this;
}
// ======================================
inline void DenseVector::add (const double factor, const DenseVector& vec) {
  assert(this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    (*this)(i) += factor*vec(i);
}
// =========================================
inline double DenseVector::dot (const DenseVector& vec) const {
  assert(this->size() == vec.size());
  double val = 0.;
  for (unsigned int i=0; i<this->size(); i++)    val += (*this)(i)*vec(i);
  return val;
}
// ===========================================
inline bool DenseVector::operator== (const DenseVector& vec) const {
  assert (this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    if ((*this)(i) != vec(i))
      return false;

  return true;
}

// ===============================================
inline bool DenseVector::operator!= (const DenseVector& vec) const {
  assert (this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    if ((*this)(i) != vec(i))
      return true;
  return false;
}

// ==================================
inline DenseVector& DenseVector::operator+= (const DenseVector& vec) {
  assert (this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    (*this)(i) += vec(i);
  return *this;
}

// ==============================
inline DenseVector& DenseVector::operator-= (const DenseVector& vec) {
  assert (this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    (*this)(i) -= vec(i);
  return *this;
}

// ===================================
inline double DenseVector::min () const {
  assert (this->size());
  double my_min = (*this)(0);
  for (unsigned int i=1; i!=this->size(); i++)    {
    double current = (*this)(i);
    my_min = (my_min < current? my_min : current);
  }
  return my_min;
}

// ==================================
inline double DenseVector::max () const {
  assert (this->size());
  double my_max = (*this)(0);
  for (unsigned int i=1; i!=this->size(); i++) {
    double current = (*this)(i);
    my_max = (my_max > current? my_max : current);
  }
  return my_max;
}

// =======================
inline double DenseVector::l1_norm () const {
  double my_norm = 0.;
  for (unsigned int i=0; i!=this->size(); i++)  my_norm += fabs((*this)(i));
  return my_norm;
}

// ==================================================
inline double DenseVector::l2_norm () const {
  double my_norm = 0.;
  for (unsigned int i=0; i!=this->size(); i++)   my_norm +=((*this)(i))*((*this)(i));
  return sqrt(my_norm);
}

// ======================================================
inline double DenseVector::linfty_norm () const {
  if (!this->size())    return 0.;
  double my_norm = fabs((*this)(0));

  for (unsigned int i=1; i!=this->size(); i++)    {
    double current = fabs((*this)(i));
    my_norm = (my_norm > current? my_norm : current);
  }
  return sqrt(my_norm);
}

// ==================================================================
inline void DenseVector::get_principal_subvector (unsigned int sub_n,
    DenseVector& dest) const {
  assert( sub_n <= this->size() );
  dest.resize(sub_n);
  for(unsigned int i=0; i<sub_n; i++)  dest(i) = (*this)(i);
}


} //end namespace femus



#endif // #ifndef __dense_vector_h__

