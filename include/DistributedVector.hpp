#ifndef __distributed_vector_h__
#define __distributed_vector_h__

#include "ParalleltypeEnum.hpp"

// C++ includes
#include <vector>
#include <algorithm>
// #include <complex>
#include <limits>

// Local includes
#include "NumericVector.hpp"
#include "Parallel.hpp"

// ======================================
// Distributed vector. Provides an interface for simple
// parallel, distributed vectors. Offers some collective
// communication capabilities.  Note that the class will
// sill function without MPI, but only on one processor.
// This lets us keep the parallel details behind the scenes.
// =====================================================
class DistributedVector : public NumericVector {
  // ==============================================
private:

  /// Actual vector datatype to hold vector entries
  std::vector<double> _values;
  /// The global vector size
  int _global_size;
  /// The local vector size
  int _local_size;

  /// The first component stored locally
  int _first_local_index;
  /// The last component (+1) stored locally
  int _last_local_index;


public:

  ///  Dummy-Constructor. Dimension=0
  explicit
  DistributedVector (const ParallelType = AUTOMATIC);

  /// Constructor. Set dimension to \p n and initialize all elements with zero.
  explicit
  DistributedVector (const  int n,const ParallelType type = AUTOMATIC);

  /// Constructor. Set local dimension to \p n_local,
  DistributedVector (const  int n,
                      const  int n_local,
                      const ParallelType type = AUTOMATIC);

  /// Constructor. Set local dimension to \p n_local
  DistributedVector (const  int N,
                      const  int n_local,
                      const std::vector< int>& ghost,
                      const ParallelType type = AUTOMATIC);

  /// Destructor, deallocates memory. Made virtual to allow
  ~DistributedVector ();

  /// Call the assemble functions
  void close ();

  /// @returns the \p DistributedVector to a pristine state.
  void clear ();

  /// Set all entries to zero. Equivalent to \p v = 0, but more obvious and faster.
  void zero ();

  /// Creates a copy of this vector and returns it in an \p AutoPtr.
  std::auto_ptr<NumericVector > clone () const;

  ///  Change the dimension of the vector to \p N. The reserved memory for
  // this vector remains unchanged if possible, to make things faster, but
  // this may waste some memory, so take this in the back of your head.
  // However, if \p N==0 all memory is freed, i.e. if you want to resize
  // the vector and release the memory not needed, you have to first call
  // \p init(0) and then \p init(N).  On \p fast==false, the vector is filled by zeros.
  void init (const  int N,
             const  int n_local,
             const bool         fast=false,
             const ParallelType type=AUTOMATIC);

  /// call init with n_local = N,
  void init (const  int N,
             const bool         fast=false,
             const ParallelType type=AUTOMATIC);

  /// Create a vector that holds tha local indices plus those specified in the \p ghost argument.
  virtual void init (const  int /*N*/,
                     const  int /*n_local*/,
                     const std::vector< int>& /*ghost*/,
                     const bool /*fast*/ = false,
                     const ParallelType = AUTOMATIC);

  /// Creates a vector that has the same dimension and storage type as \p other, including ghost dofs.
  virtual void init (const NumericVector& other,
                     const bool fast = false) {
    this->init(other.size(),other.local_size(),fast,other.type());
  }
  // ===================================
  //  SETTINGS
  // ===================================

  /// v(i) = value
  void set (const  int i, const double value);
  /// v(i) += value
  void add (const  int i, const double value);

  /// \f$U(0-N) = s\f$: fill all components.
  NumericVector & operator= (const double s);
  ///  \f$U = V\f$: copy all components.
  NumericVector & operator= (const NumericVector &V);

  ///  \f$U = V\f$: copy all components.
  DistributedVector & operator= (const DistributedVector &V);
  ///  \f$U = V\f$: copy all components.
  NumericVector & operator= (const std::vector<double> &v);

  /// \f$ U=v \f$ where v is a std::vect
  virtual void insert (const std::vector<double>& v,
                       const std::vector< int>& dof_indices);
  /// \f$ U=v \f$ where v is a NumericVect
  virtual void insert (const NumericVector& V,
                       const std::vector< int>& dof_indices);
  /// \f$ U=V \f$ where V is type
  // DenseVector and you want to specify WHERE to insert it
  virtual void insert (const DenseVector& V,
                       const std::vector< int>& dof_indices);
  /// \f$ U=V \f$ where V is type  DenseSubVector
  virtual void insert (const DenseSubVector& V,
                       const std::vector< int>& dof_indices);

// ===================================
//  RETURN FUNCTIONS
// ===================================

  /// @returns the minimum element in the vector.
  double min () const;
  /// @returns the maximum element in the vector.
  double max () const;
  /// @returns the sum of all values in the vector
  double sum() const;
  /// @returns the \f$l_1\f$-norm of the vector
  double l1_norm () const;
  /// @returns the \f$l_2\f$-norm of the vector
  double l2_norm () const;
  /// @returns the maximum absolute value of the elements
  double linfty_norm () const;

  /// @returns dimension of the vector.
  int size () const;
  /// @returns the local size of the vector (index_stop-index_start)
  int local_size() const;
  /// @returns the index of the first vector element actually stored on this processor
  int first_local_index() const;
  /// @returns the index of the last vector element
  int last_local_index() const;

  /// Access components, returns \p U(i).
  double operator() (const  int i) const;

  // ===================================
  //  ALGEBRA
  // ===================================

  /// Addition operator. Fast equivalent to \p U.add(1, V).
  NumericVector & operator += (const NumericVector &V);
  /// Subtraction operator. Fast equivalent to \p U.add(-1, V).
  NumericVector & operator -= (const NumericVector &V);

  /// Addition of \p s (scalar) to all components. Note
  void add (const double s);
  /// \f$U+=V\f$ Simple vector addition
  void add (const NumericVector& V);
  /// \f$U+=a*V\f$. Simple vector addition
  void add (const double a, const NumericVector& v);

  /// \f$U+=v\f$ where v is a \p std::vector<double>
  void add_vector (const std::vector<double>& v,
                   const std::vector< int>& dof_indices);
  /// \f$U+=V\f$ where U and V are type  NumericVector and you
  // want to specify WHERE to add the \p NumericVector V
  void add_vector (const NumericVector& V,
                   const std::vector< int>& dof_indices);
  /// \f$U+=A*V\f$. Add the product A*V
  void add_vector (const NumericVector&/* V*/, const SparseMatrix&/* A*/) {
    std::cout <<"error";
    exit(0);
  }
  void add_vector (const NumericVector& /*A*/,const SparseRectangularMatrix& /*V*/) {
    std::cout <<"error";
    exit(0);
  }

  void resid (const NumericVector &/*rhs_in*/,const NumericVector& /*x_in*/,
              const SparseMatrix& /*A_in*/) {
    std::cout <<"error";
    exit(0);
  }
  void matrix_mult (const NumericVector &/*vec_in*/,const SparseRectangularMatrix &/*mat_in*/) {
    std::cout <<"error";
    exit(0);
  }

  void matrix_mult_transpose(const NumericVector &vec_in,const SparseRectangularMatrix &mat_in){
    std::cout <<"error";
    exit(0);
  }
  
  
  /// \f$U+=V\f$ with DenseVector V
  void add_vector (const DenseVector& V,
                   const std::vector< int>& dof_indices);



  /// Scale each element of the vector by the given factor.
  void scale (const double factor);
  /// v = abs(v)... that is, each entry in v is replaced by its absolute value.
  virtual void abs();

  /// Computes the dot product, p = U.V
  virtual double dot(const NumericVector& V) const;

  /// Computes the pointwise (i.e. component-wise) product of \p vec1
  /// and \p vec2 and stores the result in \p *this.
  virtual void pointwise_mult (const NumericVector& vec1,
                               const NumericVector& vec2);
  /// Swaps the vector data and metadata
  virtual void swap (NumericVector &v);

  // ===================================
  //  DISTRIBUTED FUNCTIONS
  // ===================================
  /// Creates a copy of the global vector in the local vector \p v_local.
  void localize (std::vector<double>& v_local) const;
  /// Same, but fills a \p NumericVector instead of a \p std::vector.
  void localize (NumericVector& v_local) const;

  /// Creates a local vector \p v_local containing only
  /// information relevant to this processor, as defined by the \p send_list.
  void localize (NumericVector& v_local,
                 const std::vector< int>& send_list) const;
  ///  Updates a local vector with selected values from neighboring
  /// processors, as defined by \p send_list.
  void localize (const  int first_local_idx,
                 const  int last_local_idx,
                 const std::vector< int>& send_list);

  /// Creates a local copy of the global vector in \p v_local only
  /// on processor \p proc_id.  By default the data is sent to processor 0.
  /// This method is useful for outputting data from one processor.
  void localize_to_one (std::vector<double>& v_local,
                        const  int proc_id=0) const;

  // PRINT
/// Print the contents of the matrix in hdf5 sparse matrix format.
  void print_hdf5(const std::string /*name="NULL"*/) const {
    std::cout<<"Not implemented \n";
    abort();
  }
  void print_personal(std::ostream& os=std::cout) const {
    os<<"Not implemented \n";
    abort();
  }
};


//--------------------------------------------------------------------------
// DistributedVector inline methods
// ============================================
inline DistributedVector::DistributedVector (const ParallelType type) :
  _global_size      (0),  _local_size       (0),
  _first_local_index(0),  _last_local_index (0) {
  this->_type = type;
}
// ============================================
inline DistributedVector::DistributedVector (const  int n,
    const ParallelType type) {
  this->init(n, n, false, type);
}
// ============================================
inline DistributedVector::DistributedVector (const  int n,
    const  int n_local,
    const ParallelType type) {
  this->init(n, n_local, false, type);
}
// ============================================
inline DistributedVector::DistributedVector (const  int n,
    const  int n_local,
    const std::vector< int>& ghost,
    const ParallelType type) {
  this->init(n, n_local, ghost, false, type);
}
// ============================================
inline DistributedVector::~DistributedVector () {
  this->clear ();
}

// ============================================
inline void DistributedVector::init (const  int n,
                                      const  int n_local,
                                      const bool fast,
                                      const ParallelType type) {

  // This function must be run on all processors at once
  parallel_onlyM();
  assert (n_local <= n);

  if (type == AUTOMATIC)    {
    if (n == n_local)    this->_type = SERIAL;
    else  this->_type = PARALLEL;
  } else    this->_type = type;
  assert ((this->_type==SERIAL && n==n_local) || this->_type==PARALLEL);

  // Clear the data structures if already initialized
  if (this->initialized())    this->clear();
  // Initialize data structures
  _values.resize(n_local);
  _local_size  = n_local;
  _global_size = n;
  _first_local_index = 0;

#ifdef HAVE_MPI

  int n_proc=0, proc_id=0;
  MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size (MPI_COMM_WORLD, &n_proc);

  std::vector<int> local_sizes     (n_proc, 0);
  local_sizes[proc_id] = n_local;
  Parallel::sum(local_sizes);

  // _first_local_index is the sum of _local_size
  // for all processor ids less than ours
  for (int p=0; p<proc_id; p++)    _first_local_index += local_sizes[p];

#  ifdef DEBUG
  // Make sure all the local sizes sum up to the global
  // size, otherwise there is big trouble!
  int sum=0;
  for (int p=0; p<n_proc; p++)    sum += local_sizes[p];
  assert (sum == static_cast<int>(n));

#  endif

#else
  // No other options without MPI!
  if (n != n_local)    {
    std::cerr << "ERROR:  MPI is required for n != n_local!"
              << std::endl;
    abort();
  }

#endif

  _last_local_index = _first_local_index + n_local;
  // Set the initialized flag
  this->_is_initialized = true;
  // Zero the components unless directed otherwise
  if (!fast)    this->zero();
}
// ============================================
inline
void DistributedVector::init (const  int n,
                               const  int n_local,
                               const std::vector< int>& /*ghost*/,
                               const bool fast,
                               const ParallelType type) {
  // TODO: we shouldn't ignore the ghost sparsity pattern
  this->init(n, n_local, fast, type);
}
// ============================================
/* Default implementation for solver packages for which ghosted
   vectors are not yet implemented.  */
// void DistributedVector::init (const NumericVector& other,
//                                const bool fast)
// {
//     this->init(other.size(),other.local_size(),fast,other.type());
// }
// ============================================
inline void DistributedVector::init (const  int n,
                                      const bool fast,
                                      const ParallelType type) {
  this->init(n,n,fast,type);
}
// ============================================
inline void DistributedVector::close () {
  assert (this->initialized());
  this->_is_closed = true;
}
// ============================================
inline void DistributedVector::clear () {
  _values.clear();
  _global_size =_local_size =_first_local_index=_last_local_index = 0;
  this->_is_closed = this->_is_initialized = false;
}
// ============================================
inline void DistributedVector::zero () {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  std::fill (_values.begin(),_values.end(),0.);
}
// ============================================
inline std::auto_ptr<NumericVector > DistributedVector::clone () const {
  std::auto_ptr<NumericVector > cloned_vector (new DistributedVector);
  cloned_vector->init(*this, true);
  *cloned_vector = *this;
  return cloned_vector;
}
// ============================================
inline  int DistributedVector::size () const {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  return _global_size;
}
// ============================================
inline  int DistributedVector::local_size () const {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  return _local_size;
}
// ============================================
inline  int DistributedVector::first_local_index () const {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  return _first_local_index;
}

// ============================================
inline  int DistributedVector::last_local_index () const {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  return _last_local_index;
}

// ============================================
inline double DistributedVector::operator() (const  int i) const {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  assert ( ((i >= first_local_index()) &&
            (i <  last_local_index())) );

  return _values[i - _first_local_index];
}

// ============================================
inline void DistributedVector::set (const  int i, const double value) {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  assert (i<size());
  assert (i-first_local_index() < local_size());

  _values[i - _first_local_index] = value;
}

// ============================================
inline void DistributedVector::add (const  int i, const double value) {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  assert (i<size());
  assert (i-first_local_index() < local_size());

  _values[i - _first_local_index] += value;
}

// ============================================
inline double DistributedVector::min () const {
  // This function must be run on all processors at once
  parallel_onlyM();

  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_min = _values.size() ?
                     (double)(_values[0]) : std::numeric_limits<double>::max();
  for ( int i = 1; i < (int)_values.size(); ++i)
    local_min = std::min((double)(_values[i]), local_min);

  Parallel::min(local_min);

  return local_min;
}

// ============================================
inline double DistributedVector::max() const {
  // This function must be run on all processors at once
  parallel_onlyM();

  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_max = _values.size() ?
                     (double)(_values[0]) : -std::numeric_limits<double>::max();
  for ( int i = 1; i < (int)_values.size(); ++i)
    local_max = std::max((double)(_values[i]), local_max);

  Parallel::max(local_max);

  return local_max;
}

// ============================================
inline void DistributedVector::swap (NumericVector &other) {
  DistributedVector& v = dynamic_cast<DistributedVector&>(other);
  std::swap(_global_size, v._global_size);
  std::swap(_local_size, v._local_size);
  std::swap(_first_local_index, v._first_local_index);
  std::swap(_last_local_index, v._last_local_index);
  // This should be O(1) with any reasonable STL implementation
  std::swap(_values, v._values);
}


#endif  // #ifdef __distributed_vector_h__
