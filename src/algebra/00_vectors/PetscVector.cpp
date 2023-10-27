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

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

#include <sstream>
#include "hdf5.h"

// Local Includes
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include "DenseVector.hpp"
#include "DenseSubvector.hpp"
#include "Parallel.hpp"
#include "PetscMacro.hpp"
#include "Casts.hpp"
#include <numeric>


namespace femus {



// ============================================
///< This function returns the sum of values in a vector
double PetscVector::sum() const {
  this->_restore_array();
  assert(this->closed());
  int ierr=0;
  PetscScalar value=0.;
  ierr = VecSum(_vec, &value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// ====================================================
/// This function returns the \f$l_1\f$-norm of the vector
double PetscVector::l1_norm() const {
  this->_restore_array();
  assert(this->closed());
  int ierr=0;
  PetscReal value=0.;
  ierr = VecNorm(_vec, NORM_1, &value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// =============================================
/// This function returns the \f$l_2\f$-norm of the vector
double PetscVector::l2_norm() const {
  this->_restore_array();
  assert(this->closed());
  int ierr=0;
  PetscReal value=0.;
  ierr = VecNorm(_vec, NORM_2, &value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// ============================================
/// This function returns the maximum absolute value of the elements of this vector
double PetscVector::linfty_norm() const {
  this->_restore_array();
  assert(this->closed());
  int ierr=0;
  PetscReal value=0.;
  ierr = VecNorm(_vec, NORM_INFINITY, &value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// =======================================
NumericVector& PetscVector::operator += (const NumericVector& v) {
  this->_restore_array();
  assert(this->closed());
  this->add(1., v);
  return *this;
}

// ============================================================
NumericVector& PetscVector::operator -= (const NumericVector& v) {
  this->_restore_array();
  assert(this->closed());
  this->add(-1., v);
  return *this;
}

// =============================================================
void PetscVector::set(const  int i, const double value) {
  this->_restore_array();
  assert(i<size());
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = VecSetValues(_vec, 1, &i_val, &petsc_value, INSERT_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->_is_closed = false;
}

// ============================================================
void PetscVector::add(const  int i, const double value) {
  this->_restore_array();
  assert(i<size());
  int ierr=0;
  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues(_vec, 1, &i_val, &petsc_value, ADD_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->_is_closed = false;
}


// ===============================================================
void PetscVector::add_vector_blocked(const std::vector<double>& values,
				     const std::vector< int>& dof_indices) {
  //this->_restore_array();
  int dof_size = dof_indices.size();
  assert(values.size() == dof_size);

  int ierr = VecSetValues(_vec,dof_size,&dof_indices[0],&values[0],ADD_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

}

// ===============================================================
void PetscVector::add_vector_blocked(const std::vector<double>& values,
                                     const std::vector< unsigned >& dof_indices) {
  //this->_restore_array();
  int dof_size = dof_indices.size();
  assert(values.size() == dof_size);
  
  int ierr = VecSetValues(_vec,dof_size,(int*)&dof_indices[0],&values[0],ADD_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
                                     }

void PetscVector::insert_vector_blocked(const std::vector<double>& values,
				     const std::vector< int>& dof_indices) {
  //this->_restore_array();
  int dof_size = dof_indices.size();
  assert(values.size() == dof_size);

  int ierr = VecSetValues(_vec,dof_size,&dof_indices[0],&values[0],INSERT_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

}

// ===============================================================
void PetscVector::add_vector(const std::vector<double>& v,
                              const std::vector< int>& dof_indices) {
  this->_restore_array();
  assert(v.size() == dof_indices.size());
  for (int i=0; i<(int)v.size(); i++)    this->add(dof_indices[i], v[i]);
}

// ===========================================================
void PetscVector::add_vector(const NumericVector& V,
                              const std::vector< int>& dof_indices) {
  assert((int)V.size() == (int)dof_indices.size());
  for (int i=0; i<(int)V.size(); i++) this->add(dof_indices[i], V(i));
}

// =========================================================
void PetscVector::add_vector(const NumericVector& V_in,
                              const SparseMatrix& A_in) {
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector* V = static_cast<const PetscVector*>(&V_in);
  const PetscMatrix* A = static_cast<const PetscMatrix*>(&A_in);
  int ierr=0;
  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.
  ierr = MatMultAdd(const_cast<PetscMatrix*>(A)->mat(), V->_vec, _vec, _vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}
// ====================================================
void PetscVector::add_vector(const DenseVector& V,
                              const std::vector<unsigned int>& dof_indices) {
  assert((int)V.size() == dof_indices.size());
  for (int i=0; i<(int)V.size(); i++)  this->add(dof_indices[i], V(i));
}
// ====================================================
void PetscVector::matrix_mult(const NumericVector &vec_in,const SparseMatrix &mat_in) {
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector* v = static_cast<const PetscVector*>(&vec_in);
  const PetscMatrix* A = static_cast<const PetscMatrix*>(&mat_in);
  int ierr=0;
  A->close();
  ierr = MatMult(const_cast<PetscMatrix*>(A)->mat(),v->_vec,_vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->close();
  return;
}

/// This function computes the residual r=P^t x
// =========================================================

void PetscVector::matrix_mult_transpose(const NumericVector &vec_in,const SparseMatrix &mat_in) {
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector* v = static_cast<const PetscVector*>(&vec_in);
  const PetscMatrix* A = static_cast<const PetscMatrix*>(&mat_in);
  int ierr=0;
  A->close();
  ierr = MatMultTranspose(const_cast<PetscMatrix*>(A)->mat(),v->_vec,_vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->close();
  return;
}

/// This function computes the residual r=b-Ax
// =========================================================
void PetscVector::resid(const NumericVector& b_in,const NumericVector& x_in,
                         const SparseMatrix& A_in) {
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVector* b = static_cast<const PetscVector*>(&b_in);
  const PetscVector* x = static_cast<const PetscVector*>(&x_in);
  const PetscMatrix* A = static_cast<const PetscMatrix*>(&A_in);
  int ierr=0; /* A->close();*/
  // residual computation r=b-Ax
  ierr = MatMult(const_cast<PetscMatrix*>(A)->mat(), x->_vec, _vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecAYPX(_vec,-1,b->_vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// ====================================================
void PetscVector::add(const double v_in) {
  this->_restore_array();
  int ierr=0;
  PetscScalar* values;
  const PetscScalar v = static_cast<PetscScalar>(v_in);

  if (this->type() != GHOSTED)    {
    const int n   = static_cast<int>(this->local_size());
    const int fli = static_cast<int>(this->first_local_index());

    for (int i=0; i<n; i++) {
      ierr = VecGetArray(_vec, &values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      int ig = fli + i;
      PetscScalar value = (values[i] + v);
      ierr = VecRestoreArray(_vec, &values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetValues(_vec, 1, &ig, &value, INSERT_VALUES);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  } else {
    /* Vectors that include ghost values require a special
    handling.  */
    Vec loc_vec;
    ierr = VecGhostGetLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    int n=0;
    ierr = VecGetSize(loc_vec, &n);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    for (int i=0; i<n; i++) {
      ierr = VecGetArray(loc_vec, &values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      PetscScalar value = (values[i] + v);
      ierr = VecRestoreArray(loc_vec, &values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetValues(loc_vec, 1, &i, &value, INSERT_VALUES);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }

    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  this->_is_closed = false;
}

// ==============================================
void PetscVector::add(const NumericVector& v) {
  this->add(1., v);
}

// ====================================================
void PetscVector::add(const double a_in, const NumericVector& v_in) {
  this->_restore_array();
  int ierr = 0;
  PetscScalar a = static_cast<PetscScalar>(a_in);

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector* v = static_cast<const PetscVector*>(&v_in);
  v->_restore_array();
  assert(this->size() == v->size());
  if (this->type() != GHOSTED) {
    ierr = VecAXPY(_vec, a, v->_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else {
    Vec loc_vec;
    Vec v_loc_vec;
    ierr = VecGhostGetLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostGetLocalForm(v->_vec,&v_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecAXPY(loc_vec, a, v_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(v->_vec,&v_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
}


// ========================================================
void PetscVector::insert(const std::vector<double>& v,
                          const std::vector< int>& dof_indices) {
  assert(v.size() == dof_indices.size());
  for (int i=0; i<(int)v.size(); i++)    this->set(dof_indices[i], v[i]);
}

// =======================================================
void PetscVector::insert(const NumericVector& V,
                          const std::vector< int>& dof_indices) {
  assert((int)V.size() == (int)dof_indices.size());
  for (int i=0; i<(int)V.size(); i++)   this->set(dof_indices[i], V(i));
}

// ========================================================
void PetscVector::insert(const DenseVector& V,
                          const std::vector< int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for (int i=0; i<(int)V.size(); i++)  this->set(dof_indices[i], V(i));
}

// =========================================================
void PetscVector::insert(const DenseSubVector& V,
                          const std::vector< int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for (int i=0; i<(int)V.size(); i++)  this->set(dof_indices[i], V(i));
}

// ================================================
void PetscVector::scale(const double factor_in) {
  this->_restore_array();
  int ierr = 0;
  PetscScalar factor = static_cast<PetscScalar>(factor_in);

  if (this->type() != GHOSTED)    {
    ierr = VecScale(_vec, factor);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else   {
    Vec loc_vec;
    ierr = VecGhostGetLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScale(loc_vec, factor);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
}

// ================================================
void PetscVector::abs() {
  this->_restore_array();
  int ierr = 0;
  if (this->type() != GHOSTED)   {
    ierr = VecAbs(_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else    {
    Vec loc_vec;
    ierr = VecGhostGetLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecAbs(loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
}

// ===========================================
double PetscVector::dot(const NumericVector& V) const {
  this->_restore_array();
  // Error flag and Return value
  int ierr = 0;
  PetscScalar value=0.;
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector* v = static_cast<const PetscVector*>(&V);
  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecDot(this->_vec, v->_vec, &value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// ===============================================
NumericVector& PetscVector::operator = (const double s_in) {
  this->_restore_array();
  assert(this->closed());

  int ierr = 0;
  PetscScalar s = static_cast<PetscScalar>(s_in);

  if (this->size() != 0)    {
    if (this->type() != GHOSTED)  {
      ierr = VecSet(_vec, s);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    } else {
      Vec loc_vec;
      ierr = VecGhostGetLocalForm(_vec,&loc_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSet(loc_vec, s);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  }

  return *this;
}



// =========================================================
NumericVector& PetscVector::operator = (const NumericVector& v_in) {
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVector* v = static_cast<const PetscVector*>(&v_in);
  *this = *v;
  return *this;
}

// =====================================================
PetscVector& PetscVector::operator = (const PetscVector& v) {
  this->_restore_array();
  v._restore_array();
  assert(this->_type == v._type);
  assert(this->size() == (int)v.size());
  assert(this->local_size() == v.local_size());
  assert(this->_global_to_local_map == v._global_to_local_map);
  if ((int)v.size() != 0)    {
    int ierr = 0;
    if (this->type() != GHOSTED)  {
      ierr = VecCopy(v._vec, this->_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    } else {
      Vec loc_vec;
      Vec v_loc_vec;
      ierr = VecGhostGetLocalForm(_vec,&loc_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostGetLocalForm(v._vec,&v_loc_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecCopy(v_loc_vec, loc_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostRestoreLocalForm(v._vec,&v_loc_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  }
  return *this;
}

// =====================================================
NumericVector& PetscVector::operator = (const std::vector<double>& v) {
  this->_restore_array();
  const  int nl   = this->local_size();
  const  int ioff = this->first_local_index();
  int ierr=0;
  PetscScalar* values;
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == (int)v.size())    {
    ierr = VecGetArray(_vec, &values);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    for (int i=0; i<nl; i++) {
      values[i] =  static_cast<PetscScalar>(v[i+ioff]);
    }
    ierr = VecRestoreArray(_vec, &values);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else    {
    assert(this->local_size() == (int)v.size());
    ierr = VecGetArray(_vec, &values);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    for (int i=0; i<nl; i++)  values[i] = static_cast<PetscScalar>(v[i]);
    ierr = VecRestoreArray(_vec, &values);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  // Make sure ghost dofs are up to date
  if (this->type() == GHOSTED)   this->close();

  return *this;
}

// ============================================================
/// This function localizes data from any processor to a local vector
void PetscVector::localize(
  NumericVector& v_local_in  //  local vector
) const { // ===================================================

  this->_restore_array();
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector* v_local = static_cast<PetscVector*>(&v_local_in);
  assert(v_local != NULL);
  assert(v_local->size() == this->size());
  int ierr = 0;
  const int n = this->size();
  IS is;
  VecScatter scatter;

  // Create idx, idx[i] = i;
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
//   Utility::iota (idx.begin(), idx.end(), 0);

  // Create the index set & scatter object
  ierr = ISCreateGeneral(MPI_COMM_WORLD, n, &idx[0], PETSC_USE_POINTER, &is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = VecScatterCreate(_vec,   is,v_local->_vec, is,&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Perform the scatter
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterEnd(scatter, _vec, v_local->_vec,
                       INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Clean up
  ierr = ISDestroy(&is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTED)   v_local->close();
  return;

}

// ==========================================================
/// This function localizes data from any processor to a local vector
void PetscVector::localize(
  NumericVector& v_local_in,         //  local vector
  const std::vector< int>& send_list  //  trasfer data vactor
) const { // ===========================================
  // Workaround for a strange bug at large-scale.
  // If we have ghosting, PETSc lets us just copy the solution, and
  // doing so avoids a segfault?
  if (v_local_in.type() == GHOSTED &&  this->type() == PARALLEL) {
    v_local_in = *this;
    return;
  }
  this->_restore_array();
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector* v_local = static_cast<PetscVector*>(&v_local_in);

  assert(v_local != NULL);
  assert((int)v_local->size() == (int)this->size());
  assert((int)send_list.size() <= (int)v_local->size());

  int ierr=0;
  const  int n_sl = send_list.size();
  IS is;
  VecScatter scatter;

  std::vector<int> idx(n_sl + this->local_size());
  for (int i=0; i<n_sl; i++)   idx[i] = static_cast<int>(send_list[i]);
  for (int i = 0; i != this->local_size(); ++i)   idx[n_sl+i] = i + this->first_local_index();

  // Create the index set & scatter object
  if (idx.empty())  ierr = ISCreateGeneral(MPI_COMM_WORLD,n_sl+this->local_size(),
                             PETSC_NULL, PETSC_USE_POINTER, &is);
  else  ierr = ISCreateGeneral(MPI_COMM_WORLD,n_sl+this->local_size(),
                                 &idx[0],  PETSC_USE_POINTER,&is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterCreate(_vec,is,v_local->_vec, is,&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Perform the scatter
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterEnd(scatter, _vec, v_local->_vec,
                       INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Clean up
  ierr = ISDestroy(&is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTED)  v_local->close();
  return;

}

// ================================================================
void PetscVector::localize(
  const  int first_local_idx,        // first index
  const  int last_local_idx,         // last index
  const std::vector< int>& send_list // index list
) { // ==================================

  assert(send_list.size() <= this->size());
  assert(last_local_idx+1 <= this->size());
  // -------------------------------------------
  const unsigned int size       = this->size();
  const unsigned int local_size = (last_local_idx - first_local_idx + 1);
  int ierr=0;
  int npi=1;

  // But we do need to stay in sync for degenerate cases
  MPI_Comm_size(MPI_COMM_WORLD, &npi);
  if (npi == 1)    return;
  // Build a parallel vector, initialize it with the local parts of (*this)
  PetscVector parallel_vec;
  parallel_vec.init(size, local_size, true, PARALLEL);

  {
    // Copy part of *this into the parallel_vec -------------
    IS is;
    VecScatter scatter;
    std::vector<int> idx(local_size);
    std::iota(idx.begin(), idx.end(), first_local_idx);

    // Create the index set & scatter object
    ierr = ISCreateGeneral(MPI_COMM_WORLD, local_size,
                           local_size ? &idx[0] : NULL, PETSC_USE_POINTER, &is);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterCreate(_vec,is, parallel_vec._vec, is, &scatter);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Perform the scatter
    ierr = VecScatterBegin(scatter, _vec, parallel_vec._vec,
                           INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterEnd(scatter, _vec, parallel_vec._vec,
                         INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Clean up
    ierr = ISDestroy(&is);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterDestroy(&scatter);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } // -------------------------------------------------------

  // localize like normal
  parallel_vec.close();
  parallel_vec.localize(*this, send_list);
  this->close();

  return;
}



// =============================================================
void PetscVector::localize(std::vector<double>& v_local) const {
  this->_restore_array();

  // This function must be run on all processors at once
  parallel_only();

  int ierr=0;
  PetscScalar *values;
  const int n = this->size();
  const int nl = this->local_size();
  v_local.clear();
  v_local.resize(n, 0.);

  ierr = VecGetArray(_vec, &values);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  int ioff = first_local_index();
  for (int i=0; i<nl; i++)   v_local[i+ioff] = static_cast<double>(values[i]);
  ierr = VecRestoreArray(_vec, &values);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  Parallel::sum(v_local);   //TODO must become Parallel
}

//  ===========================================================
void PetscVector::localize_to_one(
  std::vector<double>& v_local,
  const  int pid
) const {
  this->_restore_array();
  int ierr=0;
  PetscScalar *values;
  const int n  = size();
  const int nl = local_size();
  v_local.resize(n);
  // only one processor
  if (n == nl) {
    ierr = VecGetArray(_vec, &values);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    for (int i=0; i<n; i++) v_local[i] = static_cast<double>(values[i]);
    ierr = VecRestoreArray(_vec, &values);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  // otherwise multiple processors
  else {
    int ioff = this->first_local_index();
    std::vector<double> local_values(n, 0.);
    {
      ierr = VecGetArray(_vec, &values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      for (int i=0; i<nl; i++) {
	local_values[i+ioff] = static_cast<double>(values[i]);
      }
      ierr = VecRestoreArray(_vec, &values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    ierr = MPI_Reduce(&local_values[0], &v_local[0],n,MPI_DOUBLE,MPI_SUM,pid,MPI_COMM_WORLD);
  }
}

//  ===========================================================
void PetscVector::localize_to_all(
  std::vector<double>& v_local) const {
  this->_restore_array();
  int ierr=0;
  PetscScalar *values;
  const int n  = size();
  const int nl = local_size();
  v_local.resize(n);
  // only one processor
  if (n == nl) {
    ierr = VecGetArray(_vec, &values);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    for (int i=0; i<n; i++) v_local[i] = static_cast<double>(values[i]);
    ierr = VecRestoreArray(_vec, &values);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  // otherwise multiple processors
  else {
    int ioff = this->first_local_index();
    std::vector<double> local_values(n, 0.);
    {
      ierr = VecGetArray(_vec, &values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      for (int i=0; i<nl; i++) {
	local_values[i+ioff] = static_cast<double>(values[i]);
      }
      ierr = VecRestoreArray(_vec, &values);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    for(int iproc=0;iproc<nprocs;iproc++){
      ierr = MPI_Reduce(&local_values[0], &v_local[0],n,MPI_DOUBLE,MPI_SUM,iproc,MPI_COMM_WORLD);
    }
  }
}

// =======================================================
void PetscVector::pointwise_mult(const NumericVector& vec1,
                                  const NumericVector& vec2) {
  this->_restore_array();

  int ierr = 0;

  // Convert arguments to PetscVector*.
  const PetscVector* vec1_petsc = static_cast<const PetscVector*>(&vec1);
  const PetscVector* vec2_petsc = static_cast<const PetscVector*>(&vec2);

  // Call PETSc function.
  if (this->type() != GHOSTED) {
    ierr = VecPointwiseMult(this->vec(),
                            const_cast<PetscVector*>(vec1_petsc)->vec(),
                            const_cast<PetscVector*>(vec2_petsc)->vec());
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else {
    Vec loc_vec;
    Vec v1_loc_vec;
    Vec v2_loc_vec;
    ierr = VecGhostGetLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostGetLocalForm(const_cast<PetscVector*>(vec1_petsc)->vec(),&v1_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostGetLocalForm(const_cast<PetscVector*>(vec2_petsc)->vec(),&v2_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecPointwiseMult(loc_vec,v1_loc_vec,v2_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecGhostRestoreLocalForm(const_cast<PetscVector*>(vec1_petsc)->vec(),&v1_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(const_cast<PetscVector*>(vec2_petsc)->vec(),&v2_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }

}

// ==============================================================
/// This function creates a subvector
void PetscVector::create_subvector(
  NumericVector& subvector,     //  subvector
  const std::vector< int>& rows  //  row index vector
) const {// ======================
  this->_restore_array();

  // PETSc data structures
  IS parent_is, subvector_is;
  VecScatter scatter;
  int ierr = 0;

  // Make sure the passed in subvector is really a PetscVector
  PetscVector* petsc_subvector = static_cast<PetscVector*>(&subvector);

  // If the petsc_subvector is already initialized, we assume that the
  // user has already allocated the *correct* amount of space for it.
  // If not, we use the appropriate PETSc routines to initialize it.
  if (!petsc_subvector->initialized())  {
    // Initialize the petsc_subvector to have enough space to hold
    // the entries which will be scattered into it.  Note: such an
    // init() function (where we let PETSc decide the number of local
    // entries) is not currently offered by the PetscVector
    // class.  Should we differentiate here between sequential and
    // parallel vector creation based on libMesh::n_processors() ?
    ierr = VecCreateMPI(MPI_COMM_WORLD,
                        PETSC_DECIDE,          // n_local
                        rows.size(),           // n_global
                        &(petsc_subvector->_vec));
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecSetFromOptions(petsc_subvector->_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Mark the subvector as initialized
    petsc_subvector->_is_initialized = true;
  } else {
    petsc_subvector->_restore_array();
  }

  // Use iota to fill an array with entries [0,1,2,3,4,...rows.size()]
  std::vector<int> idx(rows.size());
  std::iota(idx.begin(), idx.end(), 0);
//   Utility::iota (idx.begin(), idx.end(), 0);

  // Construct index sets
  ierr = ISCreateGeneral(MPI_COMM_WORLD,rows.size(),(int*) &rows[0],
                         PETSC_USE_POINTER,&parent_is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = ISCreateGeneral(MPI_COMM_WORLD,rows.size(),(int*) &idx[0],
                         PETSC_USE_POINTER,&subvector_is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Construct the scatter object
  ierr = VecScatterCreate(this->_vec,parent_is,petsc_subvector->_vec,
                          subvector_is,&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Actually perform the scatter
  ierr = VecScatterBegin(scatter, this->_vec, petsc_subvector->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterEnd(scatter, this->_vec, petsc_subvector->_vec,
                       INSERT_VALUES,SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Clean up
  ierr = ISDestroy(&parent_is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = ISDestroy(&subvector_is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return;
}
  
  void PetscVector::BinaryPrint(const char* fileName){
    
    PetscViewer binv;
    PetscViewerBinaryOpen(MPI_COMM_WORLD,fileName, FILE_MODE_WRITE, &binv);
    VecView(_vec, binv);
    PetscViewerDestroy(&binv);
   
    
  }
  
  void PetscVector::BinaryLoad(const char* fileName){
    
    PetscViewer binv;
    PetscViewerBinaryOpen(MPI_COMM_WORLD,fileName, FILE_MODE_READ, &binv);
    VecLoad(_vec, binv);
    PetscViewerDestroy(&binv);
    this->close();
  }

} //end namespace femus



#endif
