#ifndef __dense_matrix_baseM_h__
#define __dense_matrix_baseM_h__

 
#include "Typedefs_conf.hpp" 
 
// C++ includes
#include <cassert>
#include <iostream>

// Forward Delcarations
class DenseVectorBaseM;
// template<typename T> class DenseVectorBase<T>;



/**
 * Defines an abstract dense matrix base class for use in Finite Element-type
 * computations.  Specialized dense matrices, for example DenseSubMatrices,
 * can be derived from this class.
 */ 
class DenseMatrixBaseM
{
  
protected:
  
  // ===============================
  // Data
  // =============================
  unsigned int _m; ///< The row dimension.
  unsigned int _n;///< The column dimension.
  
  // ==================================
   // Construction/Destr
  // ===============================
  /// Constructor.  Creates a dense matrix of dimension \p m by \p n.
  /// Protected so that there is no way the user can create one.
  DenseMatrixBaseM(const unsigned int m=0,const unsigned int n=0) : _m(m), _n(n) {};
  
public: 
  ///Destructor. Empty.     
  virtual ~DenseMatrixBaseM() {};
 /// Set every element in the matrix to 0. 
  virtual void zero() = 0;
  
  // ===============================
  // Return functions
  // ===============================
  /// @returns the \p (i,j) element of the matrix.
  virtual Real el(const unsigned int i,const unsigned int j) const = 0;
  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  virtual Real & el(const unsigned int i, const unsigned int j) = 0;
  
   /// @returns the row-dimension of the matrix.
  unsigned int m() const { return _m; }
  /// @returns the column-dimension of the matrix.
  unsigned int n() const { return _n; }

  // ===============================
  // Algebra
  // ===============================	 
  /// Performs the operation: (*this) <- M2 * (*this) 
  virtual void left_multiply (const DenseMatrixBaseM& M2) = 0;
  /// Performs the operation: (*this) <- (*this) * M3
  virtual void right_multiply (const DenseMatrixBaseM& M3) = 0;
  /// Adds \p a to every element; matrix += a mat 
  void   add (const Real a, const DenseMatrixBaseM& mat);
  
  // ===============================
  // print
  // =============================== 
  /// Pretty-print the matrix to \p stdout.
  void print(std::ostream& os) const;
  /// Formatted print as above but allows you to do DenseMatrix K; std::cout << K << std::endl;
  friend std::ostream& operator << (std::ostream& os, const DenseMatrixBaseM& m){m.print(os);return os;}
  /// Prints the matrix entries with more decimal places in scientific notation.
  void print_scientific(std::ostream& os) const;
  
 
   
protected:
  ///  Performs the computation M1 = M2 * M3 where: M1 = (m x n) M2 = (m x p) M3 = (p x n)
  void multiply (DenseMatrixBaseM& M1, const DenseMatrixBaseM& M2,const DenseMatrixBaseM& M3);

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
   void condense(const unsigned int i,const unsigned int j,
		const Real val,DenseVectorBaseM& rhs);

};


// ===========================================================
//  INLINE FUNCTIONS
// ===========================================================
inline
void DenseMatrixBaseM::add (const Real factor,const DenseMatrixBaseM& mat)
{
  assert (this->m() == mat.m());assert (this->n() == mat.n());
  for (unsigned int j=0; j<this->n(); j++)
    for (unsigned int i=0; i<this->m(); i++)
      this->el(i,j) += factor*mat.el(i,j);
}


#endif // #ifndef __dense_matrix_base_h__

