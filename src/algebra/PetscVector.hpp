/*=========================================================================

 Program: FEMUS
 Module: PetscVector
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_PetscVector_hpp__
#define __femus_algebra_PetscVector_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

// C++ includes
#include <map>
#include <vector>
#include <cstdio>
// Local includes
#include "NumericVector.hpp"
#include "PetscMacro.hpp"
#include "ParalleltypeEnum.hpp"
#include "Casts.hpp"

/// Petsc include files.
EXTERN_C_FOR_PETSC_BEGIN
# include <petscvec.h>
EXTERN_C_FOR_PETSC_END

#include <mpi.h>


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class SparseMatrix;

/**
 * Petsc vector. Provides a nice interface to the
 * Petsc C-based data structures for parallel vectors.
 */

class PetscVector : public NumericVector {

public:
  // =====================================
  // Constructor /Destructor
  // =====================================
  ///  Dummy-Constructor. Dimension=0
  explicit
  PetscVector (const ParallelType type = AUTOMATIC);
  /// Constructor. Set dimension to \p n and initialize all elements with zero.
  explicit
  PetscVector (const int n,
                const ParallelType type = AUTOMATIC);
  /// Constructor. Set local dimension to \p n_local, the global dimension
  /// to \p n, and initialize all elements with zero.
  PetscVector (const int n,
                const int n_local,
                const ParallelType type = AUTOMATIC);
  /// Constructor. Set local dimension to \p n_local, the global dimension to \p n
  PetscVector (const int N,
                const int n_local,
                const std::vector<int>& ghost,
                const ParallelType type = AUTOMATIC);
  /// Constructor.  Creates a PetscVector assuming a valid PETSc Vec object
  PetscVector(Vec v);

  /// Destructor, deallocates memory
  ~PetscVector ();


  /// Creates a copy of this vector and returns it in an \p AutoPtr.
  std::unique_ptr<NumericVector > clone () const;


  /// Call the assemble functions
  void close ();
  /// This function returns the \p PetscVector to a pristine state.
  void clear ();


  /// Change the dimension of the vector to \p N. doublehe reserved memory for
  /// this vector remains unchanged if possible, to make things faster, but
  /// this may waste some memory, so take this in the back of your head.
  /// However, if \p N==0 all memory is freed, i.e. if you want to resize
  /// the vector and release the memory not needed, you have to first call
  /// \p init(0) and then \p init(N).
  void init (const int N,
             const int n_local,
             const bool         fast=false,
             const ParallelType type=AUTOMATIC);

  /// call init with n_local = N,
  void init (const int N,
             const bool         fast=false,
             const ParallelType type=AUTOMATIC);

  /// Create a vector that holds tha local indices plus those specified in the \p ghost argument.
  void init (const int /*N*/,
             const int /*n_local*/,
             const std::vector<int>& /*ghost*/,
             const bool /*fast*/ = false,
             const ParallelType = AUTOMATIC);
  /// Creates a vector that has the same dimension
  void init (const NumericVector& other,
             const bool fast = false);

  // ===========================
  // SETTING FUNCTIONS
  // ===========================

  /// v(i) = value
  void set (const int i, const double value);
  /// v(i) += value
  void add (const int i, const double value);

  /// Set all entries to zero. Equivalent to \p v = 0
  void zero ();

  /// \f$U(0-N) = s\f$: fill all components.
  NumericVector & operator= (const double s);
  ///  \f$U = V\f$: copy all components.
  NumericVector & operator= (const NumericVector &V);
  /// \f$U = V\f$: copy all components.
  PetscVector & operator= (const PetscVector &V);
  ///  \f$U = V\f$: copy all components.
  NumericVector & operator= (const std::vector<double> &v);

  /// \f$ U=v \f$ where v is a std::vector<double>
  void insert (const std::vector<double>& v,const std::vector<int>& dof_indices);
  /// \f$U=V\f$, where U and V are type NumericVector and you
  void insert(const NumericVector& V,const std::vector<int>& dof_indices);
  /// \f$ U=V \f$ where V is type DenseVector
  void insert (const DenseVector& V, const std::vector<int>& dof_indices);
  /// \f$ U=V \f$ where V is type DenseSubVector
  void insert (const DenseSubVector& V, const std::vector<int>& dof_indices);

  // ===========================
  // RETURN FUNCTIONS
  // ===========================
  /// Returns the raw PETSc vector context pointer.  Note this is generally
  /// not required in user-level code. Just don't do anything crazy like calling VecDestroy()!
  Vec vec () {
    assert (_vec != NULL);
    return _vec;
  }

  double min () const; ///< This function returns the minimum element in the vector.
  double max () const; ///< This function returns the maximum element in the vector.
  double sum () const; ///< This function returns the sum of values in a vector

  double l1_norm () const;     ///< This function returns the \f$l_1\f$-norm of the vector
  double l2_norm () const;     ///< This function returns the \f$l_2\f$-norm of the vector
  double linfty_norm () const; ///< This function returns the maximum absolute value of the elements of this vector


  int size () const;         ///< This function returns dimension of the vector
  int local_size() const;    ///< This function returns the local size of the vector (index_stop-index_start)
  int first_local_index() const;  ///< This function returns the index of the first vector element
  int last_local_index() const;   ///< This function returns the index of the last vector element

  /// Maps the global index \p i to the corresponding global index. If
  /// the index is not a ghost cell, this is done by subtraction the
  /// number of the first local index.  If it is a ghost cell, it has
  /// to be looked up in the map.
  int map_global_to_local_index(const int i) const;

  /// Access components, returns \p U(i).
  double operator() (const int i) const;


  /// Access multiple components at once. Overloaded method that
  /// should be faster (probably much faster) than calling \p
  /// operator() individually for each index.
  void get(const std::vector<int>& index, std::vector<double>& values) const;

  // ===========================
  // ALGEBRA FUNCTIONS
  // ===========================

  /// Addition operator. Fast equivalent to \p U.add(1, V).
  NumericVector & operator += (const NumericVector &V);
  /// Subtraction operator. Fast equivalent to \p U.add(-1, V).
  NumericVector & operator -= (const NumericVector &V);

  /// \f$U(0-LIBMESH_DIM)+=s\f$. Addition of \p s to all components.
  void add (const double s);
  /// \f$ U+=V \f$. Simple vector addition, equal to the
  void add (const NumericVector& V);
  /// \f$ U+=a*V \f$. Simple vector addition, equal to the
  void add (const double a, const NumericVector&  v);

  /// \f$ U+=v \f$ where \p v is a std::vector !!!fast
  void add_vector_blocked(const std::vector<double>& v,
			  const std::vector< int>& dof_indices);
  
  void add_vector_blocked(const std::vector<double>& v,
                          const std::vector< unsigned>& dof_indices);

  void insert_vector_blocked(const std::vector<double>& values,
                             const std::vector< int>& dof_indices);
  
  /// \f$ U+=v \f$ where \p v is a std::vector
  void add_vector (const std::vector<double>& v,
                   const std::vector<int>& dof_indices);
  
  /// \f$ U+=V \f$ where U and V are type  \p NumericVector and you
  void add_vector (const NumericVector& V,
                   const std::vector<int>& dof_indices);
  ///\f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
  void add_vector (const NumericVector &V,const SparseMatrix &A);
  ///\f$U+=V \f$ where U and V are type DenseVector
  void add_vector (const DenseVector& V,const std::vector<unsigned int>& dof_indices);
  void resid (const NumericVector & rhs_in,const NumericVector& x,const SparseMatrix& A);
  /// \f$U+=A*V\f$, add the product A*v
  void matrix_mult (const NumericVector &vec_in,const SparseMatrix &mat_in);

  void matrix_mult_transpose (const NumericVector &vec_in,const SparseMatrix &mat_in);
  /// Scale each element of the vector by the given factor.
  void scale (const double factor);

  /// v = abs(v)... that is, each entry in v is replaced
  void abs();

  /// Computes the dot product, p = U.V
  double dot(const NumericVector& V) const;
  /// Computes the pointwise (i.e. component-wise) product of \p vec1
  /// and \p vec2 and stores the result in \p *this.
  void pointwise_mult (const NumericVector& vec1,
                       const NumericVector& vec2);


  // ===========================
  // PARALLEL OPERATIONS
  // ===========================

  /// Creates a copy of the global vector in the local vector \p v_local.
  void localize (std::vector<double>& v_local) const;
  /// Same, but fills a \p NumericVector instead of a \p std::vector.
  void localize (NumericVector& v_local) const;

  /// Creates a local vector \p v_local containing
  void localize (NumericVector& v_local,const std::vector<int>& send_list) const;
  /// Updates a local vector with selected values from neighboring processors
  void localize (const int first_local_idx,
                 const int last_local_idx,
                 const std::vector<int>& send_list);
  /// Creates a local copy of the global vector in
  /// \p v_local only on processor \p proc_id.  By
  /// default the data is sent to processor 0.
  void localize_to_one (std::vector<double>& v_local,
                        const int proc_id=0) const;

  void localize_to_all (std::vector<double>& v_local) const;

  /// Creates a "subvector" from this vector using the rows indices of the "rows" array.
  void create_subvector(NumericVector& subvector,
                        const std::vector<int>& rows) const;

  /// Swaps the raw PETSc vector context pointers.
  void swap (NumericVector &v);

  void BinaryPrint(const char* fileName);
  void BinaryLoad(const char* fileName);

protected:

  /// Queries the array (and the local form if the vector is ghosted) from Petsc.
  void _get_array(void) const;
  ///  Restores the array (and the local form if the vector is ghosted) to Petsc.
  void _restore_array(void) const;

private:

  /// Actual Petsc vector datatype
  Vec _vec;

  /// If \p true, the actual Petsc array of the values of the vector is currently accessible.
  /// doublehat means that the members \p _local_form and \p _values are valid.
  mutable bool _array_is_present;

  /// Petsc vector datatype to hold the local form of a ghosted vector.
  /// doublehe contents of this field are only valid if the vector is
  /// ghosted and \p _array_is_present is \p true.
  mutable Vec _local_form;

  /// Pointer to the actual Petsc array of the values of the vector.
  /// doublehis pointer is only valid if \p _array_is_present is \p true.
  mutable PetscScalar* _values;

  /// doubleype for map that maps global to local ghost cells.
  typedef std::map<int,int> GlobalToLocalMap;

  /// Map that maps global to local ghost cells (will be empty if not in ghost cell mode
  GlobalToLocalMap _global_to_local_map;

  /// doublehis boolean value should only be set to false
  /// for the constructor which takes a PETSc Vec object.
  bool _destroy_vec_on_exit;

#ifndef NDEBUG
  ///Size of the local form, for being used in assertations.  doublehe
  /// contents of this field are only valid if the vector is ghosted
  /// and \p _array_is_present is \p true.
  mutable int _local_size;
#endif

};


/**
 * ----------------------- Inline functions ----------------------------------
 */

inline PetscVector::PetscVector (const ParallelType type)
  : _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
    _global_to_local_map(),
    _destroy_vec_on_exit(true) {
  this->_type = type;
}

inline PetscVector::PetscVector (const int n,const ParallelType type)
  : _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
    _global_to_local_map(),
    _destroy_vec_on_exit(true) {
  this->init(n, n, false, type);
}

inline PetscVector::PetscVector (const int n,
                                 const int n_local,
                                 const ParallelType type)
  : _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
    _global_to_local_map(),
    _destroy_vec_on_exit(true) {
  this->init(n, n_local, false, type);
}

inline PetscVector::PetscVector (const int n,
                                 const int n_local,
                                 const std::vector<int>& ghost,
                                 const ParallelType type)
  : _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
    _global_to_local_map(),
    _destroy_vec_on_exit(true) {
  this->init(n, n_local, ghost, false, type);
}

inline PetscVector::PetscVector (Vec v)
  : _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
    _global_to_local_map(),
    _destroy_vec_on_exit(false) {
  this->_vec = v;
  this->_is_closed = true;
  this->_is_initialized = true;

  /* We need to ask PETSc about the (local to global) ghost value
     mapping and create the inverse mapping out of it.  */
  int ierr=0;
  int petsc_local_size=0;
  ierr = VecGetLocalSize(_vec, &petsc_local_size);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

//   // Get the vector type from PETSc.
//   const VecType type;
//   ierr = VecGetType(_vec, &type);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);

    // Get the vector type from PETSc.
  // As of Petsc 3.0.0, the VecType #define lost its const-ness, so we
  // need to have it in the code
#if PETSC_VERSION_LESS_THAN(3,0,0) || !PETSC_VERSION_LESS_THAN(3,4,0)
  // Pre-3.0 and petsc-dev (as of October 2012) use non-const versions
  VecType ptype;
#else
  const VecType ptype;
#endif
  ierr = VecGetType(_vec, &ptype);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  if((strcmp(ptype,VECSHARED) == 0) || (strcmp(ptype,VECMPI) == 0)) {
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,1,1)
    ISLocalToGlobalMapping mapping = _vec->mapping;
#else
    ISLocalToGlobalMapping mapping;
    ierr = VecGetLocalToGlobalMapping(_vec, &mapping);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
//     ISLocalToGlobalMapping mapping;
//     ierr = VecGetLocalToGlobalMapping(_vec, &mapping);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);

    // If is a sparsely stored vector, set up our new mapping
    if (mapping) {
      const unsigned int local_size = static_cast<unsigned int>(petsc_local_size);
      const unsigned int ghost_begin = static_cast<unsigned int>(petsc_local_size);
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,4,0)
      const numeric_index_type ghost_end = static_cast<numeric_index_type>(mapping->n);
#else
      PetscInt n;
      ierr = ISLocalToGlobalMappingGetSize(mapping, &n);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      const unsigned int ghost_end = static_cast<unsigned int>(n);
#endif
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,1,1)
      const PetscInt *indices = mapping->indices;
#else
      const PetscInt *indices;
      ierr = ISLocalToGlobalMappingGetIndices(mapping,&indices);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
      for(unsigned int i=ghost_begin; i<ghost_end; i++)
        _global_to_local_map[indices[i]] = i-local_size;
      this->_type = GHOSTED;
#if !PETSC_VERSION_RELEASE || !PETSC_VERSION_LESS_THAN(3,1,1)
      ierr = ISLocalToGlobalMappingRestoreIndices(mapping, &indices);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
//       const unsigned int ghost_end = static_cast<unsigned int>(mapping->n);
//       const int *indices = mapping->indices;
//       ierr = ISLocalToGlobalMappingGetIndices(mapping,&indices);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
//       for(unsigned int i=ghost_begin; i<ghost_end; i++)
//         _global_to_local_map[indices[i]] = i-local_size;
//       this->_type = GHOSTED;
    }
    else
      this->_type = PARALLEL;
  }
  else
    this->_type = SERIAL;

  this->close();
}


inline PetscVector::~PetscVector () {
  this->clear ();
}

inline void PetscVector::init (const int n,
                                const int n_local,
                                const bool fast,
                                const ParallelType type) {
  int ierr=0;
  int petsc_n=static_cast<int>(n);
  int petsc_n_local=static_cast<int>(n_local);
  // Clear initialized vectors
  if (this->initialized())  this->clear();
  if (type == AUTOMATIC)    {
    if (n == n_local)   this->_type = SERIAL;
    else    this->_type = PARALLEL;
  } else    this->_type = type;

  assert ((this->_type==SERIAL && n==n_local) || this->_type==PARALLEL);

  // create a sequential vector if on only 1 processor
  if (this->_type == SERIAL) {
    ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n, &_vec);
    CHKERRABORT(PETSC_COMM_SELF,ierr);
    ierr = VecSetFromOptions (_vec);
    CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  // otherwise create an MPI-enabled vector
  else if (this->_type == PARALLEL) {
    assert (n_local <= n);
    ierr = VecCreateMPI (MPI_COMM_WORLD, petsc_n_local, petsc_n, &_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecSetFromOptions (_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else   {
    std::cout << "Not good" <<std::endl;
    abort();
  }
  this->_is_initialized = true;
  this->_is_closed = true;
  if (fast == false)  this->zero ();
}


inline void PetscVector::init (const int n,
                                const bool fast,
                                const ParallelType type) {
  this->init(n,n,fast,type);
}


inline void PetscVector::init (const int n,
                               const int n_local,
                               const std::vector<int>& ghost,
                               const bool fast,
                               const ParallelType type) {
  int ierr=0;
  PetscInt petsc_n=static_cast<int>(n);
  PetscInt petsc_n_local=static_cast<int>(n_local);
  PetscInt petsc_n_ghost=static_cast<int>(ghost.size());

  // If the mesh is not disjoint, every processor will either have
  // all the dofs, none of the dofs, or some non-zero dofs at the
  // boundary between processors.
  //
  // However we can't assert this, because someone might want to
  // construct a GHOSTED vector which doesn't include neighbor element
  // dofs.  Boyce tried to do so in user code, and we're going to want
  // to do so in System::project_vector().
  //
  // libmesh_assert(n_local == 0 || n_local == n || !ghost.empty());

  assert(sizeof(PetscInt) == sizeof(int));
  // If the mesh is disjoint, the following assertion will fail.
  // If the mesh is not disjoint, every processor will either have
  // all the dofs, none of the dofs, or some non-zero dofs at the
  // boundary between processors.
  //assert(n_local == 0 || n_local == n || !ghost.empty());

  PetscInt* petsc_ghost = ghost.empty() ? PETSC_NULL :
                     const_cast<int*>(reinterpret_cast<const PetscInt*>(&ghost[0]));

  // Clear initialized vectors
  if (this->initialized())   this->clear();

  assert(type == AUTOMATIC || type == GHOSTED);
  this->_type = GHOSTED;

  /* Make the global-to-local ghost cell map.  */
  for (int i=0; i<(int)ghost.size(); i++){
    _global_to_local_map[ghost[i]] = i;
  }

  /* Create vector.  */
  ierr = VecCreateGhost (MPI_COMM_WORLD, petsc_n_local, petsc_n,
                         petsc_n_ghost, petsc_ghost, &_vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = VecSetFromOptions (_vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  this->_is_initialized = true;
  this->_is_closed = true;
  if (fast == false)
    this->zero ();
}


inline void PetscVector::init(const NumericVector& other, const bool fast) {
  // Clear initialized vectors
  if (this->initialized())   this->clear();
  const PetscVector& v = libmeshM_cast_ref<const PetscVector&>(other);
  // Other vector should restore array.
  if (v.initialized())    {
    v._restore_array();
  }

  this->_global_to_local_map = v._global_to_local_map;
  this->_is_closed      = v._is_closed;
  this->_is_initialized = v._is_initialized;
  this->_type = v._type;

  if (v.size() != 0)   {
    int ierr = 0;
    ierr = VecDuplicate (v._vec, &this->_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  if (fast == false)   this->zero ();
}


inline void PetscVector::close () {
  this->_restore_array();
  int ierr=0;

  ierr = VecAssemblyBegin(_vec);  					CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecAssemblyEnd(_vec);  					CHKERRABORT(MPI_COMM_WORLD,ierr);

  if (this->type() == GHOSTED) {
    ierr = VecGhostUpdateBegin(_vec,INSERT_VALUES,SCATTER_FORWARD);  	CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostUpdateEnd(_vec,INSERT_VALUES,SCATTER_FORWARD);  	CHKERRABORT(MPI_COMM_WORLD,ierr);

  }
  this->_is_closed = true;
}


inline void PetscVector::clear () {
  if (this->initialized())    this->_restore_array();

  if ((this->initialized()) && (this->_destroy_vec_on_exit))    {
    int ierr=0;
    ierr = VecDestroy(&_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  this->_is_closed = this->_is_initialized = false;
  _global_to_local_map.clear();
}


inline void PetscVector::zero (){
    this->_restore_array();
    int ierr=0;
    PetscScalar z=0.;

    if (this->type() != GHOSTED) {
        ierr = VecSet (_vec, z);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    else {
        /* Vectors that include ghost values require a special
        handling.  */
        Vec loc_vec; ierr = VecGhostGetLocalForm (_vec,&loc_vec);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
        ierr = VecSet (loc_vec, z);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
        ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
}


inline std::unique_ptr<NumericVector > PetscVector::clone () const {
  std::unique_ptr<NumericVector> cloned_vector (new PetscVector);
  cloned_vector->init(*this, true);
  *cloned_vector = *this;
  return cloned_vector;
}


inline int PetscVector::size () const {
  assert (this->initialized());
  int ierr=0, petsc_size=0;
  if (!this->initialized())   return 0;
  ierr = VecGetSize(_vec, &petsc_size);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<int>(petsc_size);
}


inline
int PetscVector::local_size () const {
  assert (this->initialized());
  int ierr=0, petsc_size=0;
  ierr = VecGetLocalSize(_vec, &petsc_size);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<int>(petsc_size);
}


inline int PetscVector::first_local_index () const {
  assert (this->initialized());
  int ierr=0, petsc_first=0, petsc_last=0;
  ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<int>(petsc_first);
}


inline int PetscVector::last_local_index () const {
  assert (this->initialized());
  int ierr=0, petsc_first=0, petsc_last=0;
  ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<int>(petsc_last);
}


inline int PetscVector::map_global_to_local_index (const int i) const {
  assert (this->initialized());

  int ierr=0, petsc_first=0, petsc_last=0;
  ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  const int first = static_cast<int>(petsc_first);
  const int last = static_cast<int>(petsc_last);

  if ((i>=first) && (i<last))    {
    return i-first;
  }

  GlobalToLocalMap::const_iterator it = _global_to_local_map.find(i);
  assert (it!=_global_to_local_map.end());
  return it->second+last-first;
}


inline double PetscVector::operator() (const int i) const {
  this->_get_array();
  const int local_index = this->map_global_to_local_index(i);
#ifndef NDEBUG
    if (this->type() == GHOSTED) assert(local_index<_local_size);
#endif
  return static_cast<double>(_values[local_index]);

//   double value;
//   VecGetValues(_vec,1,&i,&value);
//
//   return value;

}


inline void PetscVector::get(const std::vector<int>& index, std::vector<double>& values) const {
  this->_get_array();

  const int num = index.size();
  values.resize(num);

  for (int i=0; i<num; i++) {
    const int local_index = this->map_global_to_local_index(index[i]);
#ifndef NDEBUG
    if (this->type() == GHOSTED) assert(local_index<_local_size);
#endif
    values[i] = static_cast<double>(_values[local_index]);
  }
}

inline double PetscVector::min () const {
  this->_restore_array();
  int index=0, ierr=0;
  PetscReal min=0.;
  ierr = VecMin (_vec, &index, &min);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // this return value is correct: VecMin returns a PetscReal
  return static_cast<double>(min);
}


inline double PetscVector::max() const {
  this->_restore_array();
  int index=0, ierr=0;
  PetscReal max=0.;
  ierr = VecMax (_vec, &index, &max);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // this return value is correct: VecMax returns a PetscReal
  return static_cast<double>(max);
}


inline void PetscVector::swap (NumericVector &other) {
  NumericVector::swap(other);
  PetscVector& v = libmeshM_cast_ref<PetscVector&>(other);
  std::swap(_vec, v._vec);
  std::swap(_destroy_vec_on_exit, v._destroy_vec_on_exit);
  std::swap(_global_to_local_map, v._global_to_local_map);
  std::swap(_array_is_present, v._array_is_present);
  std::swap(_local_form, v._local_form);
  std::swap(_values, v._values);
}


inline void PetscVector::_get_array(void) const {
  assert (this->initialized());
  if (!_array_is_present) {
    int ierr=0;
    if (this->type() != GHOSTED) {
      ierr = VecGetArray(_vec, &_values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    } else {
      ierr = VecGhostGetLocalForm (_vec,&_local_form);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGetArray(_local_form, &_values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
#ifndef NDEBUG
      int local_size = 0;
      ierr = VecGetLocalSize(_local_form, &local_size);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      _local_size = static_cast<int>(local_size);
#endif
    }
    _array_is_present = true;
  }
}


inline void PetscVector::_restore_array(void) const {
  assert (this->initialized());
  if (_array_is_present) {
    int ierr=0;
    if (this->type() != GHOSTED) {
      ierr = VecRestoreArray (_vec, &_values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      _values = NULL;
    } else	{
      ierr = VecRestoreArray (_local_form, &_values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      _values = NULL;
      ierr = VecGhostRestoreLocalForm (_vec,&_local_form);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      _local_form = NULL;
#ifndef NDEBUG
      _local_size = 0;
#endif
    }
    _array_is_present = false;
  }
}

} //end namespace femus


#endif
#endif

