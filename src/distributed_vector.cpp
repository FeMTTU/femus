// C++ includes
#include <cmath> // for std::abs
#include <cassert>

// Local Includes
#include "DistributedVector.hpp"
#include "DenseVector.hpp"
#include "DenseSubvector.hpp"
#include "Parallel.hpp"



//--------------------------------------------------------------------------
// DistributedVector methods

// =====================================
double DistributedVector::sum () const {
  // This function must be run on all processors at once
  parallel_onlyM();

  assert (this->initialized());
  assert ((int)(int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  double local_sum = 0.;
  for ( int i=0; i<local_size(); i++)  local_sum += _values[i];
  Parallel::sum(local_sum);
  return local_sum;
}

// ========================================
double DistributedVector::l1_norm () const {
  // This function must be run on all processors at once
  parallel_onlyM();
  assert (this->initialized());
  assert ((int)(int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  double local_l1 = 0.;
  for ( int i=0; i<local_size(); i++)  local_l1 += std::abs(_values[i]);
  Parallel::sum(local_l1);
  return local_l1;
}

// ==========================================
double DistributedVector::l2_norm () const {
  // This function must be run on all processors at once
  parallel_onlyM();
  assert (this->initialized());
  assert ((int)(int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  double local_l2 = 0.;
  for ( int i=0; i<local_size(); i++)  local_l2 += (_values[i]*_values[i]);
  Parallel::sum(local_l2);
  return std::sqrt(local_l2);
}

// =============================================
double DistributedVector::linfty_norm () const {
  // This function must be run on all processors at once
  parallel_onlyM();
  assert (this->initialized());
  assert ((int)(int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  double local_linfty = 0.;
  for ( int i=0; i<local_size(); i++)
    local_linfty  = std::max(local_linfty,
                             static_cast<double>(std::abs(_values[i]))
                            ); // Note we static_cast so that both
  // types are the same, as required
  // by std::max
  Parallel::max(local_linfty);
  return local_linfty;
}

// ===============================================
NumericVector& DistributedVector::operator += (const NumericVector& v) {
  assert (this->closed());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  add(1., v);
  return *this;
}

// ============================================
NumericVector& DistributedVector::operator -= (const NumericVector& v) {
  assert (this->closed());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  add(-1., v);
  return *this;
}

// ============================================
void DistributedVector::add_vector (const std::vector<double>& v,
                                     const std::vector< int>& dof_indices) {
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<(int)(int)v.size(); i++)    add (dof_indices[i], v[i]);
}

// ============================================
void DistributedVector::add_vector (const NumericVector& V,
                                     const std::vector< int>& dof_indices) {
  assert ((int)V.size() == (int)dof_indices.size());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  for ( int i=0; i<(int)V.size(); i++)  add (dof_indices[i], V(i));
}

// ============================================
void DistributedVector::add_vector (const DenseVector& V,
                                     const std::vector< int>& dof_indices) {
  assert ((int)V.size() == dof_indices.size());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<(int)V.size(); i++)    add (dof_indices[i], V(i));
}

// ============================================
void DistributedVector::add (const double v) {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<local_size(); i++)    _values[i] += v;
}

// ============================================
void DistributedVector::add (const NumericVector& v) {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  add (1., v);
}

// ============================================
void DistributedVector::add (const double a, const NumericVector& v) {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  add(a, v);
}

// ============================================
void DistributedVector::insert (const std::vector<double>& v,
                                 const std::vector< int>& dof_indices) {
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<(int)v.size(); i++)    this->set (dof_indices[i], v[i]);
}

// ============================================
void DistributedVector::insert (const NumericVector& V,
                                 const std::vector< int>& dof_indices) {
  assert ((int)V.size() == (int)dof_indices.size());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<(int)V.size(); i++)    this->set (dof_indices[i], V(i));
}

// ============================================
void DistributedVector::insert (const DenseVector& V,
                                 const std::vector< int>& dof_indices) {
  assert ((int)V.size() == dof_indices.size());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<(int)V.size(); i++) this->set (dof_indices[i], V(i));
}

// ============================================
void DistributedVector::insert (const DenseSubVector& V,
                                 const std::vector< int>& dof_indices) {
  assert ((int)V.size() == dof_indices.size());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<(int)V.size(); i++)  this->set (dof_indices[i], V(i));
}

// ============================================
void DistributedVector::scale (const double factor) {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<local_size(); i++)   _values[i] *= factor;
}

// ============================================
void DistributedVector::abs() {
  assert (this->initialized());
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<local_size(); i++) this->set(i,std::abs(_values[i]));
}

// ============================================
double DistributedVector::dot (const NumericVector& V) const {
  // This function must be run on all processors at once
  parallel_onlyM();

  // Make sure the NumericVector passed in is really a DistributedVector
  const DistributedVector* v = dynamic_cast<const DistributedVector*>(&V);

  // Make sure that the two vectors are distributed in the same way.
  assert ( this->first_local_index() == v->first_local_index() );
  assert ( this->last_local_index()  == v->last_local_index()  );

  // The result of dotting together the local parts of the vector.
  double local_dot = 0;

  for ( int i=0; i<this->local_size(); i++)
    local_dot += this->_values[i] * v->_values[i];

  // The local dot products are now summed via MPI
  Parallel::sum(local_dot);
  return local_dot;
}

// ============================================
NumericVector&
DistributedVector::operator = (const double s) {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<local_size(); i++) _values[i] = s;
  return *this;
}

// ============================================
NumericVector&
DistributedVector::operator = (const NumericVector& v_in) {
  // Make sure the NumericVector passed in is really a DistributedVector
  const DistributedVector* v = static_cast<const DistributedVector*>(&v_in);
  *this = *v;
  return *this;
}

// ============================================
DistributedVector&
DistributedVector::operator = (const DistributedVector& v) {
  this->_is_initialized    = v._is_initialized;
  this->_is_closed         = v._is_closed;

  _global_size       = v._global_size;
  _local_size        = v._local_size;
  _first_local_index = v._first_local_index;
  _last_local_index  = v._last_local_index;

  if (v.local_size() == this->local_size()) {
    _values = v._values;
  } else  {
    abort();
  }

  return *this;
}

// ============================================
NumericVector&
DistributedVector::operator = (const std::vector<double>& v) {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  if ((int)v.size() == local_size())    _values = v;

  else if ((int)v.size() == size())
    for ( int i=first_local_index(); i<last_local_index(); i++)
      _values[i-first_local_index()] = v[i];
  else    {
    abort();
  }

  return *this;
}

// ============================================
void DistributedVector::localize (NumericVector& v_local_in) const {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  DistributedVector* v_local = dynamic_cast<DistributedVector*>(&v_local_in);
  v_local->_first_local_index = 0;
  v_local->_global_size = v_local->_local_size = v_local->_last_local_index = size();
  v_local->_is_initialized =    v_local->_is_closed = true;

  // Call localize on the vector's values.  This will help
  // prevent code duplication
  localize (v_local->_values);

#ifndef LIBMESH_HAVE_MPI
  assert (local_size() == size());
#endif
}

// ============================================
void DistributedVector::localize (NumericVector& v_local_in,
                                   const std::vector< int>&) const {
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  // TODO: We don't yet support the send list; this is inefficient:
  localize (v_local_in);
}

// ============================================
void DistributedVector::localize (const  int first_local_idx,
                                   const  int last_local_idx,
                                   const std::vector< int>& send_list) {
  // Only good for serial vectors
  assert (this->size() == this->local_size());
  assert (last_local_idx > first_local_idx);
  assert ((int)send_list.size() <= this->size());
  assert (last_local_idx < this->size());

  const  int size       = this->size();
  const  int local_size = (last_local_idx - first_local_idx + 1);

  // Don't bother for serial cases
  if ((first_local_idx == 0) &&      (local_size == size))    return;

  // Build a parallel vector, initialize it with the local
  // parts of (*this)
  DistributedVector parallel_vec;

  parallel_vec.init (size, local_size, true, PARALLEL);

  // Copy part of *this into the parallel_vec
  for ( int i=first_local_idx; i<=last_local_idx; i++)
    parallel_vec._values[i-first_local_idx] = _values[i];

  // localize like normal
  parallel_vec.localize (*this, send_list);
}

// ============================================
void DistributedVector::localize (std::vector<double>& v_local) const {
  // This function must be run on all processors at once
  parallel_onlyM();

  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local = this->_values;
  Parallel::allgather (v_local);

#ifndef LIBMESH_HAVE_MPI
  assert (local_size() == size());
#endif
}

// ============================================
void DistributedVector::localize_to_one (std::vector<double>& v_local,
    const  int pid) const {
  // This function must be run on all processors at once
  parallel_onlyM();

  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local = this->_values;
//   Parallel::gather (pid, v_local);

#ifndef LIBMESH_HAVE_MPI
  assert (local_size() == size());
#endif
}

// ============================================
void DistributedVector::pointwise_mult (const NumericVector&,
    const NumericVector&)
//void DistributedVector::pointwise_mult (const NumericVector& vec1,
//					   const NumericVector& vec2)
{
  std::cout << " Not implemented";
  abort();
}

