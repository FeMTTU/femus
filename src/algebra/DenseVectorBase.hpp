/*=========================================================================

 Program: FEMUS
 Module: DenseVectorBase
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_DenseVectorBase_hpp__
#define __femus_algebra_DenseVectorBase_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <iostream>


namespace femus {


// Local Includes

// =====================================================
// Defines an abstract dense vector base class for use in
// Finite Element-type computations. Specialized dense vectors,
// for example DenseSubVectors, can be derived from this class.
// ================================================

// ======================================
// DenseVectorBase class definition

class DenseVectorBase {
public:
  /// Constructor.  Empty.
  DenseVectorBase() {}

  /// Destructor.  Does nothing.
  virtual ~DenseVectorBase() {}

  /// Set every element in the vector to 0.  Needs to
  virtual void zero() = 0;

  // ============================
  // RETURN FUNCTIONS
  // ===========================
  /// * @returns the \p (i) element of the vector.
  virtual double el(const unsigned int i) const = 0;
  /// * @returns the \p (i) element of the vector as a writeable reference.
  virtual double & el(const unsigned int i) = 0;

  /// @returns the size of the vector.
  virtual unsigned int size() const = 0;

  // ============================
  // PRINT FUNCTIONS
  // ===========================
  ///  Pretty-print the vector to \p stdout.
  void print(std::ostream& os=std::cout) const;
  /// Same as above, but allows you to print using the
  friend std::ostream& operator << (std::ostream& os, const DenseVectorBase& v) {
    v.print(os);
    return os;
  }
  /// Prints the entries of the vector with additional
  void print_scientific(std::ostream& os) const;
};





} //end namespace femus



#endif // #ifndef __dense_vector_base_h__

