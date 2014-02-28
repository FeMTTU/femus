#ifndef __dense_vector_baseM_h__
#define __dense_vector_baseM_h__

#include <iostream>
// Local Includes

#include "Typedefs_conf.hpp"
// =====================================================
// Defines an abstract dense vector base class for use in
// Finite Element-type computations. Specialized dense vectors,
// for example DenseSubVectors, can be derived from this class.
// ================================================
 

// ======================================
// DenseVectorBase class definition

class DenseVectorBaseM
{
public:
  /// Constructor.  Empty.
  DenseVectorBaseM() {}
  
  /// Destructor.  Does nothing.    
  virtual ~DenseVectorBaseM() {}

  /// Set every element in the vector to 0.  Needs to
  virtual void zero() = 0;
  
  // ============================
  // RETURN FUNCTIONS
  // ===========================
  /// * @returns the \p (i) element of the vector.
  virtual Real el(const unsigned int i) const = 0;
  /// * @returns the \p (i) element of the vector as a writeable reference.
  virtual Real & el(const unsigned int i) = 0;

  /// @returns the size of the vector.
  virtual unsigned int size() const = 0; 
  
  // ============================
  // PRINT FUNCTIONS
  // ===========================
  ///  Pretty-print the vector to \p stdout.
  void print(std::ostream& os=std::cout) const;
  /// Same as above, but allows you to print using the
  friend std::ostream& operator << (std::ostream& os, const DenseVectorBaseM& v){
    v.print(os);    return os;
  }  
  /// Prints the entries of the vector with additional
  void print_scientific(std::ostream& os) const;
};




#endif // #ifndef __dense_vector_base_h__

