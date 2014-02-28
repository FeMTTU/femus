#ifndef __dense_submatrixM_h__
#define __dense_submatrixM_h__

// C++ includes
#include <cassert>

#include "Typedefs_conf.hpp"

// Local Includes
#include "dense_matrix_baseM.hpp"
#include "dense_matrixM.hpp"
#include "dense_subvectorM.hpp"



/**
 * Defines a dense submatrix for use in Finite Element-type computations.
 * Useful for storing element stiffness matrices before summation
 * into a global matrix, particularly when you have systems of equations.
 *
 * @author Benjamin S. Kirk, 2003
 */ 

// ------------------------------------------------------------
// DenseSubMatrixM class definition

class DenseSubMatrixM : public DenseMatrixBaseM
{
  
private:
  
  /// The parent matrix that contains this submatrix.
  DenseMatrixM& _parent_matrix;
  /// The row offset into the parent matrix.
  unsigned int _i_off;
  /// The column offset into the parent matrix.
  unsigned int _j_off;
  
public:

  /**
   * Constructor.  Creates a dense submatrix of the matrix
   * \p parent.  The submatrix has dimensions \f$(m \times n)\f$,
   * and the \f$(0,0) entry of the submatrix is located
   * at the \f$(ioff,joff)\f$ location in the parent matrix.
   */
  DenseSubMatrixM(DenseMatrixM& parent,
		 const unsigned int ioff=0,
		 const unsigned int joff=0,
		 const unsigned int m=0,
		 const unsigned int n=0) : 
     DenseMatrixBaseM(m,n), _parent_matrix(parent){
  this->reposition (ioff, joff, m, n);
};

  /// Copy Constructor.   
DenseSubMatrixM(const DenseSubMatrixM& other_matrix)
  : DenseMatrixBaseM(other_matrix._m, other_matrix._n),
    _parent_matrix(other_matrix._parent_matrix){
  _i_off = other_matrix._i_off;  _j_off = other_matrix._j_off; 
}
  /// Destructor     
  virtual ~DenseSubMatrixM() {};

    /// Set every element in the submatrix to 0.
  virtual void zero(){
  for (unsigned int i=0; i<this->m(); i++)
    for (unsigned int j=0; j<this->n(); j++)
      _parent_matrix(i + this->i_off(),j + this->j_off()) = 0.;
}
  // =====================================
  // Return functions 
  // =====================================
  /// @returns a reference to the parent matrix.
  DenseMatrixM& parent () { return _parent_matrix; }
  
  ///  @returns the \p (i,j) element of the submatrix
  Real operator() (const unsigned int i,const unsigned int j) const;
  ///   @returns the \p (i,j) element of the submatrix as a writeable reference.
  Real & operator() (const unsigned int i,const unsigned int j){
  assert (i < this->m());  assert (j < this->n());
  assert (i + this->i_off() < _parent_matrix.m());  assert (j + this->j_off() < _parent_matrix.n());
  return _parent_matrix (i + this->i_off(),j + this->j_off());
}

  
  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  virtual Real el(const unsigned int i,const unsigned int j) const { return (*this)(i,j); }
  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  virtual Real & el(const unsigned int i,const unsigned int j)     { return (*this)(i,j); } 
  
   /// @returns the row offset into the parent matrix.
  unsigned int i_off() const { return _i_off; }
  /// @returns the column offset into the parent matrix.
  unsigned int j_off() const { return _j_off; }

  // ==================================================
  // Operations
  // ==================================================
  /// Performs the operation: (*this) <- M2 * (*this) 
  virtual void left_multiply (const DenseMatrixBaseM& M2);

  /// Performs the operation: (*this) <- (*this) * M3
  virtual void right_multiply (const DenseMatrixBaseM& M3);
  
  /// Changes the location of the submatrix in the parent matrix. 
  void reposition(const unsigned int ioff,  const unsigned int joff,
		  const unsigned int m,  const unsigned int n);

  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   */
  void condense(const unsigned int i,const unsigned int j,
		const Real val,DenseSubVectorM& rhs){
    this->parent().condense(this->i_off()+i,
			    this->j_off()+j,
			    val, rhs.parent());
  }

};


// -------------------------------------------------- 
// Constructor



  



   
inline
Real DenseSubMatrixM::operator () (const unsigned int i,
				  const unsigned int j) const
{
  assert (i < this->m());
  assert (j < this->n());
  assert (i + this->i_off() < _parent_matrix.m());
  assert (j + this->j_off() < _parent_matrix.n());
  
  return _parent_matrix (i + this->i_off(),
			 j + this->j_off());
}

inline
   void DenseSubMatrixM::reposition(const unsigned int ioff,
				   const unsigned int joff,
				   const unsigned int m,
				   const unsigned int n)
{				   
  _i_off = ioff;
  _j_off = joff;
  this->_m = m;
  this->_n = n;

  // Make sure we still fit in the parent matrix.
  assert ((this->i_off() + this->m()) <= _parent_matrix.m());
  assert ((this->j_off() + this->n()) <= _parent_matrix.n());
}




#endif // #ifndef __dense_matrix_h__

