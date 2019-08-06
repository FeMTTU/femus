/*=========================================================================

 Program: FEMuS
 Module: NumericVector
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMuS
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_NumericVector_hpp__
#define __femus_algebra_NumericVector_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "SolverPackageEnum.hpp"
#include "ParalleltypeEnum.hpp"
#include "FemusConfig.hpp"

// C++ includes
#include <vector>
#include <set>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <iostream>


namespace femus {


//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class DenseVector;
class DenseSubVector;
class SparseMatrix;

/**
 * Numeric vector. Provides a uniform interface
 * to vector storage schemes for different linear
 * algebra libraries.
 */

class NumericVector {


public:
  
  /** Dummy-Constructor. Dimension=0 */
  explicit
  NumericVector (const ParallelType = AUTOMATIC);

  /** Constructor. Set dimension to \p n and initialize all elements with zero. */
  explicit
  NumericVector (const  int n,
                  const ParallelType = AUTOMATIC);

  /** Constructor. Set local dimension to \p n_local, the global dimension to \p n */
  NumericVector (const  int n,
                  const  int n_local,
                  const ParallelType = AUTOMATIC);

  /** Constructor. Set local dimension to \p n_local, the global dimension to \p n */
  NumericVector (const  int N,
                  const  int n_local,
                  const std::vector< int>& ghost,
                  const ParallelType = AUTOMATIC);

  /** Builds a \p NumericVector using the linear solver package */
  /** specified by \p solver_package */
  static std::unique_ptr<NumericVector>
  build(const SolverPackage solver_package = LSOLVER);
  
  /** Creates a copy of this vector and returns it in an \p AutoPtr. */
  virtual std::unique_ptr<NumericVector > clone () const = 0;

  /** Destructor, deallocates memory. */
  virtual ~NumericVector () {
    clear ();
  }
  
  /** @returns the \p NumericalVectorM to a pristine state. */
  virtual void clear () {
    _is_closed= false;
    _is_initialized = false;
  }

  /** Call the assemble functions */
  virtual void close () = 0;
  
  /**
   * Change the dimension of the vector to \p N. The reserved memory for
   * this vector remains unchanged if possible, to make things faster, but
   * this may waste some memory, so take this in the back of your head.
   * However, if \p N==0 all memory is freed, i.e. if you want to resize
   * the vector and release the memory not needed, you have to first call
   * \p init(0) and then \p init(N). 
   */
  virtual void init (const  int,
                     const  int,
                     const bool = false,
                     const ParallelType = AUTOMATIC) = 0;

  /** call init with n_local = N, */
  virtual void init (const  int,
                     const bool = false,
                     const ParallelType = AUTOMATIC) = 0;

  /** Create a vector that holds tha local indices plus those specified
   * in the \p ghost argument. 
   */
  virtual void init (const  int /*N*/,
                     const  int /*n_local*/,
                     const std::vector< int>& /*ghost*/,
                     const bool /*fast*/ = false,
                     const ParallelType = AUTOMATIC) = 0;

  /** Creates a vector that has the same dimension and storage type as
   *\p other, including ghost dofs.
   */
  virtual void init (const NumericVector& other,
                     const bool fast = false) = 0;


///  Creates the subvector "subvector" from the indices in the
/// "rows" array.  Similar to the create_submatrix routine for
/// the SparseMatrix class, it is currently only implemented for
/// PetscVectors.
  virtual void create_subvector(NumericVector&,const std::vector< int>&) const {
    std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
    exit(0);
  }

  // =====================================
  // SETTING FUNCTIONS
  // =====================================
  /// v(i) = value
  virtual void set (const  int i, const double value) = 0;
  /// v(i) += value
  virtual void add (const  int i, const double value) = 0;

/// Set all entries to zero. Equivalent to \p v = 0
  virtual void zero () = 0;
  /// \f$U(0-N) = s\f$: fill all components.
  virtual NumericVector & operator= (const double s) = 0;
  ///  \f$U = V\f$: copy all components.
  virtual NumericVector & operator= (const NumericVector &V) = 0;
  ///  \f$U = V\f$: copy all components.
  virtual NumericVector & operator= (const std::vector<double> &v) = 0;

  /// \f$ U=v \f$ where v is a std::vector<double> andspecify WHERE to insert it
  virtual void insert (const std::vector<double>& v,
                       const std::vector< int>& dof_indices) = 0;
  /// \f$U=V\f$, and specify WHERE to insert
  virtual void insert (const NumericVector& V,
                       const std::vector< int>& dof_indices) = 0;
  /// \f$ U=V \f$ insert
  virtual void insert (const DenseVector& V,
                       const std::vector< int>& dof_indices) = 0;
  /// \f$ U=V \f$ and specify WHERE to insert it
  virtual void insert (const DenseSubVector& V,
                       const std::vector< int>& dof_indices) = 0;
  // =====================================
  // RETURN FUNCTIONS
  // =====================================
  /// @returns true if the vector has been initialized,
  virtual bool initialized() const {
    return _is_initialized;
  }
  /// @returns true if the vector is closed and ready
  virtual bool closed() const {
    return _is_closed;
  }

  /// @returns the type (SERIAL, PARALLEL, GHOSTED) of the vector.
  ParallelType type() const {
    return _type;
  }
  /// @returns the type (SERIAL, PARALLEL, GHOSTED) of the vector.
  ParallelType & type() {
    return _type;
  }

  /// @returns the minimum element in the vector.
  virtual double min () const = 0;
  /// @returns the maximum element in the vector.
  virtual double max () const = 0;
  /// returns the sum of the elements in a vector
  virtual double sum() const = 0;

  /// @returns the \f$l_1\f$-norm of the vector, i.e.
  virtual double l1_norm () const = 0;
  /// @returns the \f$l_2\f$-norm of the vector, i.e.
  virtual double l2_norm () const = 0;
  /// @returns the maximum absolute value of the
  virtual double linfty_norm () const = 0;

  /// @returns the \f$l_1\f$-norm of the vector, i.e.
  virtual double subset_l1_norm (const std::set< int> & indices);
  /// @returns the \f$l_2\f$-norm of the vector, i.e.
  virtual double subset_l2_norm (const std::set< int> & indices);
  /// @returns the maximum absolute value of the
  virtual double subset_linfty_norm (const std::set< int> & indices);


  /// @returns dimension of the vector.
  virtual  int size () const = 0;
  /// @returns the local size of the vector (index_stop-index_start).
  virtual  int local_size() const = 0;
  /// @returns the index of the first vector element
  virtual  int first_local_index() const = 0;
  /// @returns the index+1 of the last vector element
  virtual  int last_local_index() const = 0;

  ///Access components, returns \p U(i).
  virtual double operator() (const  int i) const = 0;
  /// @returns the element \p U(i)
  virtual double el(const  int i) const {
    return (*this)(i);
  }

  /**
   * Access multiple components at once.  \p values will be resized,
   * if necessary, and filled.  The default implementation calls \p
   * operator() for each index, but some implementations may supply
   * faster methods here.
   */
  virtual void get(const std::vector< int>& index, std::vector<double>& values) const;

  // =====================================
  // algebra FUNCTIONS
  // =====================================

  /// Addition operator. Fast equivalent to \p U.add(1, V).
  virtual NumericVector & operator += (const NumericVector &V) = 0;
  /// Subtraction operator. Fast equivalent to \p U.add(-1, V).
  virtual NumericVector & operator -= (const NumericVector &V) = 0;
  /// Multiplication operator. Equivalent to \p U.scale(a)
  NumericVector & operator *= (const double a) {
    this->scale(a);
    return *this;
  }
  /// Division operator. Equivalent to \p U.scale(1./a)
  NumericVector & operator /= (const double a) {
    this->scale(1./a);
    return *this;
  }

  /// \f$U(0-LIBMESH_DIM)+=s\f$. Addition of \p s to all components. Note
  virtual void add (const double s) = 0;
  /// \f$U+=V\f$: Simple vector addition, equal to the
  virtual void add (const NumericVector& V) = 0;
  /// \f$U+=a*V\f$. Simple vector addition, equal to the
  virtual void add (const double a, const NumericVector& v) = 0;

  /// \f$ U+=v \f$ where \p v is a std::vector !!!fast
  virtual void add_vector_blocked(const std::vector<double>& v,
			const std::vector< int>& dof_indices) =0;
      
  virtual void add_vector_blocked(const std::vector<double>& v,
                                  const std::vector< unsigned>& dof_indices) = 0;
            
  virtual void insert_vector_blocked(const std::vector<double>& v,
			const std::vector< int>& dof_indices) =0;
  
  /// \f$ U+=v \f$ where v is a DenseVector
  virtual void add_vector (const std::vector<double>& v,
                           const std::vector< int>& dof_indices) = 0;
  /// \f$U+=V\f$, where U and V are type
  virtual void add_vector (const NumericVector& V,
                           const std::vector< int>& dof_indices) = 0;
  /// \f$U+=A*V\f$, add the product A*v
  virtual void add_vector (const NumericVector& /*A*/,const SparseMatrix& /*V*/) = 0;
  virtual void resid (const NumericVector &/*rhs_in*/,const NumericVector& /*A*/,const SparseMatrix& /*V*/) = 0;
  virtual void matrix_mult (const NumericVector &vec_in,const SparseMatrix &mat_in) = 0;
  virtual void matrix_mult_transpose(const NumericVector &vec_in,const SparseMatrix &mat_in) = 0; 
  /// \f$U+=A*V\f$, add the product of a \p ShellMatrix \p
//   void add_vector (const NumericVector& v,
// 		   const ShellMatrix<double>& a);
//
  /// \f$ U+=V \f$ where U and V are type
  virtual void add_vector (const DenseVector& V,
                           const std::vector<unsigned int>& dof_indices) = 0;

  /// Scale each element
  virtual void scale (const double factor) = 0;
  /// v = abs(v)... that is, each entry in v is replaced
  virtual void abs() = 0;
  /// Computes the dot product, p = U.V
  virtual double dot(const NumericVector&) const = 0;

  /// Exchanges the values/sizes of two vectors.
  virtual void swap (NumericVector &v);

  // =====================================
  // PARALLEL OPERATIONS
  // =====================================

  /// Creates a copy of the global vector in the local vector \p v_local.
  virtual void localize (std::vector<double>& v_local) const = 0;
  /// Same, but fills a \p NumericVector instead
  virtual void localize (NumericVector& v_local) const = 0;

  /// Creates a local vector \p v_local
  virtual void localize (NumericVector& v_local,
                         const std::vector< int>& send_list) const = 0;
  /// Updates a local vector with selected values from neighboring
  virtual void localize (const  int first_local_idx,
                         const  int last_local_idx,
                         const std::vector< int>& send_list) = 0;
  /// Creates a local copy of the global vector
  virtual void localize_to_one (std::vector<double>& v_local,
                                const  int proc_id=0) const = 0;
  /// Creates a local copy of the global vector
  virtual void localize_to_all (std::vector<double>& v_local) const = 0;
  /// @returns \p -1 when \p this is equivalent to \p other_vector,
  virtual int compare (const NumericVector &other_vector,
                       const double threshold = 1.e-20) const;
  /// Computes the pointwise (i.e. component-wise) product of \p vec1
  virtual void pointwise_mult (const NumericVector& vec1,
                               const NumericVector& vec2) = 0;


  // PRINTING FUNCTIONS
  /** Prints the local contents of the vector to the screen. */
  virtual void print(std::ostream& os=std::cout) const;
  
  /** Prints the global contents of the vector to the screen. */
  virtual void print_global(std::ostream& os=std::cout) const;
  
  /** Same as above but allows you to use stream syntax. */
  friend std::ostream& operator << (std::ostream& os, const NumericVector& v) {
    v.print_global(os);
    return os;
  }
  
  virtual void BinaryPrint(const char* fileName){
    std::cout << "BinaryPrint is not available for this vector type\n";
    abort();
  };
  
  virtual void BinaryLoad(const char* fileName){
    std::cout << "BinaryLoad is not available for this vector type\n";
    abort();
  };
  
protected:

  // member data
  /** Flag to see if the Numeric assemble routines have been called yet */
  bool _is_closed;
  
  /** Flag to tell if init has been called yet */
  bool _is_initialized;
  
  /** Type of vector */
  ParallelType _type;

};

/**
 *----------------------- Inline functions ----------------------------------
 */


inline NumericVector::NumericVector (const ParallelType type) :
  _is_closed(false),  _is_initialized(false),  _type(type) {}
  

inline NumericVector::NumericVector (const  int /*n*/,
                                       const ParallelType type) :
  _is_closed(false),_is_initialized(false), _type(type) {
  std::cout<< "Abstract base class! ";
  exit(0); // Abstract base class!
}


inline NumericVector::NumericVector (const int /*n*/,const int /*n_local*/,
                                       const ParallelType type) :
  _is_closed(false),  _is_initialized(false),  _type(type) {
  std::cout<< "Abstract base class! ";
  exit(0); // Abstract base class!
}


inline NumericVector::NumericVector (const int /*n*/,const int /*n_local*/,
                                       const std::vector<int>& /*ghost*/,
                                       const ParallelType type) :
  _is_closed(false),  _is_initialized(false),  _type(type) {
  std::cout<< "Abstract base class! ";
  exit(0); // Abstract base class!
}


inline void NumericVector::get(const std::vector<int>& index,
                                std::vector<double>& values) const {
  const  int num = index.size();
  values.resize(num);
  for( int i=0; i<num; i++) values[i] = (*this)(index[i]);
}


inline void  NumericVector::swap (NumericVector &v) {
  std::swap(_is_closed, v._is_closed);
  std::swap(_is_initialized, v._is_initialized);
  std::swap(_type, v._type);
}


inline void NumericVector::print(std::ostream& os) const {
  assert (this->initialized());
  os << "Size\tglobal =  " << this->size()
     << "\t\tlocal =  " << this->local_size() << std::endl;
  os << "#\tValue" << std::endl;
  for ( int i=this->first_local_index(); i<this->last_local_index(); i++)
    os << i << "\t" << (*this)(i) << std::endl;
}


inline void NumericVector::print_global(std::ostream& os) const {
  assert (this->initialized());
  std::vector<double> v(this->size());
  this->localize(v);
  // Right now we only want one copy of the output
#ifdef HAVE_MPIM   
  if(static_cast<int>(libMeshPrivateData::_processor_id)) return ;
#else
  if(0)return;
#endif
  os << "Size\tglobal =  " << this->size() << std::endl;
  os << "#\tValue" << std::endl;
  for (unsigned int i=0; i!=v.size(); i++) os << i << "\t" << v[i] << std::endl;
  
}


} //end namespace femus



#endif  // #ifdef __numeric_vector_h__
