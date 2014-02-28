#ifndef __dense_subvectorM_h__
#define __dense_subvectorM_h__


#include "Typedefs_conf.hpp"
#include "dense_vectorM.hpp"

// ===========================================
// Defines a dense subvector for use in Finite Element-type computations.
// Useful for storing element load vectors  before summation
// into a global vector, particularly when you have systems of equations.
 //=====================================

// ------------------------------------------------------------
// DenseSubVector class definition
class DenseSubVectorM : public DenseVectorBase
{
  // =====================================
  // DATA
  // =====================================
  private:
    
  /// The parent vector that contains this subvector.
  DenseVectorM& _parent_vector;
  /// The length of this subvector.
  unsigned int _n;
  /// The offset into the parent vector.
  unsigned int _i_off;
  
public:
  // =====================================
  // CONSTRUCTOR /DESTRUCTOR
  // =====================================
  /// Constructor.  Creates a dense subvector of the vector
  DenseSubVectorM(DenseVectorM& parent,const unsigned int ioff=0, const unsigned int n=0);

  /// Destructor.  Does nothing.     
  virtual ~DenseSubVectorM() {}
  
  /// @returns a reference to the parent vector.
  DenseVectorM& parent () { return _parent_vector; }
  
  /// Set every element in the subvector to 0.
  virtual void zero();

  // ====================================
  // RETURN FUNCTIONS
  // ====================================
  
  /// @returns the \p (i,j) element of the subvector.
  Real operator() (const unsigned int i) const;

  /// @returns the \p (i,j) element of the subvector as a writeable reference.
  Real & operator() (const unsigned int i);

  /// @returns the \p (i) element of the vector.
  virtual Real el(const unsigned int i) const { return (*this)(i); }

  /// @returns the \p (i) element of the vector as a writeable reference.
  virtual Real & el(const unsigned int i)     { return (*this)(i); }
  
  ///  @returns the size of the subvector.
  virtual unsigned int size() const { return _n; }

  ///  @returns the row offset into the parent vector.
  unsigned int i_off() const { return _i_off; }

  /// Changes the location of the subvector in the parent vector. 
  void reposition(const unsigned int ioff,const unsigned int n);
  

};



// ------------------------------------------------------------
// Dense Vector member functions
inline
DenseSubVectorM::DenseSubVectorM(DenseVectorM& parent,
				  const unsigned int ioff,
				  const unsigned int n) :
  _parent_vector(parent)
{
  reposition (ioff, n);
}



inline void DenseSubVectorM::reposition(const unsigned int ioff,
				   const unsigned int n)
{				   
  _i_off = ioff;
  _n = n;
  // Make sure we still fit in the parent vector.
  assert ((this->i_off() + this->size()) <= _parent_vector.size());
}

inline void DenseSubVectorM::zero(){
  for (unsigned int i=0; i<this->size(); i++)    _parent_vector (i + this->i_off()) = 0.;
}

inline Real DenseSubVectorM::operator () (const unsigned int i) const
{
  assert (i < this->size());
  assert (i + this->i_off() < _parent_vector.size());
  return _parent_vector (i + this->i_off());
}

// =========================================
inline Real & DenseSubVectorM::operator () (const unsigned int i)
{
  assert (i < this->size());
  assert (i + this->i_off() < _parent_vector.size());
  
  return _parent_vector (i + this->i_off());
}



#endif // #ifndef __dense_vector_h__

