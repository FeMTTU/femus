#include "FemusExtLib_conf.hpp"

#ifdef FEMUS_HAVE_PETSC

#include <sstream>

#include "hdf5.h"

#include "Typedefs_conf.hpp"
// C++ includes
// Local Includes
#include "petsc_vectorM.hpp"
#include "petsc_matrixM.hpp"
#include "DenseSubvector.hpp"
#include "DenseVector.hpp"
#include "parallelM.hpp"
#include "petsc_macroM.hpp"

#include "Casts.hpp"  // TODO #include "utility.h"  

//-----------------------------------------------------------------------
// PetscVector members

// void PetscVectorM::init (const NumericVectorM& v, const bool fast)
// {
//   libmesh_error();
//   init (v.local_size(), v.size(), fast);
//   vec = libmesh_cast_ref<const PetscVectorM&>(v).vec;
// }

// ============================================
Real PetscVectorM::sum () const{
  this->_restore_array(); assert(this->closed());
  int ierr=0;  PetscScalar value=0.;
  ierr = VecSum (_vec, &value);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<Real>(value);
}

// ====================================================
Real PetscVectorM::l1_norm () const{
  this->_restore_array();  assert(this->closed());
  int ierr=0;  PetscReal value=0.;
  ierr = VecNorm (_vec, NORM_1, &value); CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<Real>(value);
}

// =============================================
Real PetscVectorM::l2_norm () const{
  this->_restore_array();  assert(this->closed());
  int ierr=0;  PetscReal value=0.;
  ierr = VecNorm (_vec, NORM_2, &value);CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<Real>(value);
}

// ============================================
Real PetscVectorM::linfty_norm () const{
  this->_restore_array(); assert(this->closed());
  int ierr=0;  PetscReal value=0.;
  ierr = VecNorm (_vec, NORM_INFINITY, &value); CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<Real>(value);
}

// =======================================
NumericVectorM& PetscVectorM::operator += (const NumericVectorM& v){
  this->_restore_array(); assert(this->closed());
  this->add(1., v);
  return *this;
}

// ============================================================
NumericVectorM& PetscVectorM::operator -= (const NumericVectorM& v){
  this->_restore_array();  assert(this->closed());
  this->add(-1., v);
  return *this;
}

// =============================================================
void PetscVectorM::set (const unsigned int i, const Real value){
  this->_restore_array();  assert(i<size());
  int ierr=0;  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, INSERT_VALUES);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->_is_closed = false;
}

// ============================================================
void PetscVectorM::add (const unsigned int i, const Real value){
  this->_restore_array();  assert(i<size());
  int ierr=0;  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues (_vec, 1, &i_val, &petsc_value, ADD_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->_is_closed = false;
}

// ===============================================================
void PetscVectorM::add_vector (const std::vector<Real>& v,
				 const std::vector<unsigned int>& dof_indices){
  this->_restore_array(); assert (v.size() == dof_indices.size());
  for (unsigned int i=0; i<v.size(); i++)    this->add (dof_indices[i], v[i]);
}

// ===========================================================
void PetscVectorM::add_vector (const NumericVectorM& V,
				 const std::vector<unsigned int>& dof_indices){
  assert (V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++) this->add (dof_indices[i], V(i));
}

// =========================================================
void PetscVectorM::add_vector (const NumericVectorM& V_in,
				 const SparseMatrixM& A_in){
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVectorM* V = static_cast<const PetscVectorM*>(&V_in);
  const PetscMatrixM* A = static_cast<const PetscMatrixM*>(&A_in);
  int ierr=0;
  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.  
  ierr = MatMultAdd(const_cast<PetscMatrixM*>(A)->mat(), V->_vec, _vec, _vec);
         CHKERRABORT(MPI_COMM_WORLD,ierr); 
}

// ====================================================
void PetscVectorM::add_vector (const DenseVector& V,
				 const std::vector<unsigned int>& dof_indices){
  assert (V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)  this->add (dof_indices[i], V(i));
}
// ====================================================
void PetscVectorM::matrix_mult (const NumericVectorM &vec_in,const SparseMatrixM &mat_in){
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVectorM* v = static_cast<const PetscVectorM*>(&vec_in);
  const PetscMatrixM* A = static_cast<const PetscMatrixM*>(&mat_in);
  int ierr=0;
  A->close();

  ierr = MatMult(const_cast<PetscMatrixM*>(A)->mat(),v->_vec,_vec);
   CHKERRABORT(MPI_COMM_WORLD,ierr);
   return;
}

/// This function computes the residual r=b-Ax
// =========================================================
void PetscVectorM::resid (const NumericVectorM& b_in,const NumericVectorM& x_in,
				 const SparseMatrixM& A_in){
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVectorM* b = static_cast<const PetscVectorM*>(&b_in);
  const PetscVectorM* x = static_cast<const PetscVectorM*>(&x_in);
  const PetscMatrixM* A = static_cast<const PetscMatrixM*>(&A_in);
  int ierr=0; /* A->close();*/
  // residual computation r=b-Ax
   ierr = MatMult(const_cast<PetscMatrixM*>(A)->mat(), x->_vec, _vec);CHKERRABORT(MPI_COMM_WORLD,ierr);
   ierr = VecAYPX( _vec,-1,b->_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// ====================================================
void PetscVectorM::add (const Real v_in){
  this->_restore_array();  int ierr=0;
  PetscScalar* values;  const PetscScalar v = static_cast<PetscScalar>(v_in);  

  if(this->type() != GHOSTEDM)    {
      const int n   = static_cast<int>(this->local_size());
      const int fli = static_cast<int>(this->first_local_index());
      
      for (int i=0; i<n; i++)	{
	  ierr = VecGetArray (_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
	  int ig = fli + i;  PetscScalar value = (values[i] + v);
	  ierr = VecRestoreArray (_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecSetValues (_vec, 1, &ig, &value, INSERT_VALUES); CHKERRABORT(MPI_COMM_WORLD,ierr); 
	}
    }
  else    {
      /* Vectors that include ghost values require a special
	 handling.  */
      Vec loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);

      int n=0;
      ierr = VecGetSize(loc_vec, &n);  CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      for (int i=0; i<n; i++){
	  ierr = VecGetArray (loc_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
	  PetscScalar value = (values[i] + v);
	  ierr = VecRestoreArray (loc_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecSetValues (loc_vec, 1, &i, &value, INSERT_VALUES);  CHKERRABORT(MPI_COMM_WORLD,ierr); 
	}

      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  this->_is_closed = false;
}

// ==============================================
void PetscVectorM::add (const NumericVectorM& v){
  this->add (1., v);
}

// ====================================================
void PetscVectorM::add (const Real a_in, const NumericVectorM& v_in){
  this->_restore_array();
  int ierr = 0;  PetscScalar a = static_cast<PetscScalar>(a_in);

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVectorM* v = static_cast<const PetscVectorM*>(&v_in);
  v->_restore_array();  assert(this->size() == v->size());
  if(this->type() != GHOSTEDM)  {
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//       // 2.2.x & earlier style
//       ierr = VecAXPY(&a, v->_vec, _vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else  // 2.3.x & later style
      ierr = VecAXPY(_vec, a, v->_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
    }
  else    {
      Vec loc_vec;     Vec v_loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostGetLocalForm (v->_vec,&v_loc_vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//       // 2.2.x & earlier style
//       ierr = VecAXPY(&a, v_loc_vec, loc_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else    // 2.3.x & later style
      ierr = VecAXPY(loc_vec, a, v_loc_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
      ierr = VecGhostRestoreLocalForm (v->_vec,&v_loc_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
}


// ========================================================
void PetscVectorM::insert (const std::vector<Real>& v,
			     const std::vector<unsigned int>& dof_indices){
  assert (v.size() == dof_indices.size());
  for (unsigned int i=0; i<v.size(); i++)    this->set (dof_indices[i], v[i]);
}

// =======================================================
void PetscVectorM::insert (const NumericVectorM& V,
			     const std::vector<unsigned int>& dof_indices){
  assert (V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)   this->set (dof_indices[i], V(i));
}

// ========================================================
void PetscVectorM::insert (const DenseVector& V,
			     const std::vector<unsigned int>& dof_indices){
  assert (V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)  this->set (dof_indices[i], V(i));
}

// =========================================================
void PetscVectorM::insert (const DenseSubVector& V,
			     const std::vector<unsigned int>& dof_indices){
  assert (V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)  this->set (dof_indices[i], V(i));
}

// ================================================
void PetscVectorM::scale (const Real factor_in){
  this->_restore_array();  int ierr = 0;
  PetscScalar factor = static_cast<PetscScalar>(factor_in);
  
  if(this->type() != GHOSTEDM)    {
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//       // 2.2.x & earlier style
//       ierr = VecScale(&factor, _vec);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
      // 2.3.x & later style	 
      ierr = VecScale(_vec, factor);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
    }
  else   {
      Vec loc_vec; ierr = VecGhostGetLocalForm (_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);

// #if PETSC_VERSION_LESS_THAN(2,3,0)
//       // 2.2.x & earlier style
//       ierr = VecScale(&factor, loc_vec);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
      // 2.3.x & later style	 
      ierr = VecScale(loc_vec, factor); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
}

// ================================================
void PetscVectorM::abs(){
  this->_restore_array();
  int ierr = 0;
  if(this->type() != GHOSTEDM)   {
      ierr = VecAbs(_vec);
      CHKERRABORT(MPI_COMM_WORLD,ierr);  
    }
  else    {
      Vec loc_vec; ierr = VecGhostGetLocalForm (_vec,&loc_vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecAbs(loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);  
      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
}

// ===========================================
Real PetscVectorM::dot (const NumericVectorM& V) const{
  this->_restore_array();
  // Error flag and Return value
  int ierr = 0; PetscScalar value=0.;
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVectorM* v = static_cast<const PetscVectorM*>(&V);
  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecDot(this->_vec, v->_vec, &value);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<Real>(value);
}

// ===============================================
NumericVectorM& PetscVectorM::operator = (const Real s_in){
  this->_restore_array();  assert(this->closed());

  int ierr = 0;  PetscScalar s = static_cast<PetscScalar>(s_in);

  if (this->size() != 0)    {
      if(this->type() != GHOSTEDM)	{
// #if PETSC_VERSION_LESS_THAN(2,3,0)
// 	  // 2.2.x & earlier style
// 	  ierr = VecSet(&s, _vec);
// 	  CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
	  // 2.3.x & later style	 
	  ierr = VecSet(_vec, s); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
	}
      else
	{
	  Vec loc_vec; ierr = VecGhostGetLocalForm (_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);

// #if PETSC_VERSION_LESS_THAN(2,3,0)
// 	  // 2.2.x & earlier style
// 	  ierr = VecSet(&s, loc_vec);
// 	  CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
	  // 2.3.x & later style	 
	  ierr = VecSet(loc_vec, s); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif

	  ierr = VecGhostRestoreLocalForm (_vec,&loc_vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
	}
    }
  
  return *this;
}



// =========================================================
NumericVectorM& PetscVectorM::operator = (const NumericVectorM& v_in){
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVectorM* v = static_cast<const PetscVectorM*>(&v_in);
  *this = *v;
  return *this;
}

// =====================================================
PetscVectorM& PetscVectorM::operator = (const PetscVectorM& v){
  this->_restore_array();  v._restore_array();
  assert (this->_type == v._type);  assert (this->size() == v.size());
  assert (this->local_size() == v.local_size());  assert (this->_global_to_local_map == v._global_to_local_map);
  if (v.size() != 0)    {
      int ierr = 0;
      if(this->type() != GHOSTEDM)	{
	  ierr = VecCopy (v._vec, this->_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	}
      else {
	  Vec loc_vec;	  Vec v_loc_vec;
	  ierr = VecGhostGetLocalForm (_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecGhostGetLocalForm (v._vec,&v_loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecCopy (v_loc_vec, loc_vec);	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecGhostRestoreLocalForm (v._vec,&v_loc_vec);	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	}
    }
  return *this;
}

// =====================================================
NumericVectorM& PetscVectorM::operator = (const std::vector<Real>& v){
  this->_restore_array();
  const unsigned int nl   = this->local_size();
  const unsigned int ioff = this->first_local_index();
  int ierr=0;
  PetscScalar* values;
      
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == v.size())    {
      ierr = VecGetArray (_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
      for (unsigned int i=0; i<nl; i++)	values[i] =  static_cast<PetscScalar>(v[i+ioff]);
      ierr = VecRestoreArray (_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else    {
      assert (this->local_size() == v.size());
      ierr = VecGetArray (_vec, &values);   CHKERRABORT(MPI_COMM_WORLD,ierr);
      for (unsigned int i=0; i<nl; i++)	values[i] = static_cast<PetscScalar>(v[i]);
      ierr = VecRestoreArray (_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  // Make sure ghost dofs are up to date
  if (this->type() == GHOSTEDM)   this->close();

  return *this;
}

// ============================================================
void PetscVectorM::localize (NumericVectorM& v_local_in) const  {
  
  this->_restore_array();

  // Make sure the NumericVector passed in is really a PetscVector
  PetscVectorM* v_local = static_cast<PetscVectorM*>(&v_local_in);
  
  assert (v_local != NULL);  assert (v_local->size() == this->size());
  int ierr = 0;  const int n = this->size();
  IS is;  VecScatter scatter;
  // Create idx, idx[i] = i;
  std::vector<int> idx(n); 
  iota (idx.begin(), idx.end(), 0);//   Utility::iota (idx.begin(), idx.end(), 0);
  // Create the index set & scatter object
  ierr = ISCreateGeneral(MPI_COMM_WORLD, n, &idx[0], PETSC_USE_POINTER, &is); CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterCreate(_vec, is, v_local->_vec, is, &scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Perform the scatter
  // API argument order change in PETSc 2.3.3
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,INSERT_VALUES, SCATTER_FORWARD);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterEnd  (scatter, _vec, v_local->_vec,INSERT_VALUES, SCATTER_FORWARD);
         CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Clean up
  ierr = ISDestroy (&is);CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTEDM)  v_local->close();
}

// ==========================================================
void PetscVectorM::localize (NumericVectorM& v_local_in,
			       const std::vector<unsigned int>& send_list) const{
  this->_restore_array();
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVectorM* v_local = static_cast<PetscVectorM*>(&v_local_in);

  assert (v_local != NULL);  assert (v_local->size() == this->size());
  assert (send_list.size() <= v_local->size());
  
  int ierr=0;  const unsigned int n_sl = send_list.size();
  IS is;  VecScatter scatter;

  std::vector<int> idx(n_sl + this->local_size());
  for (unsigned int i=0; i<n_sl; i++)   idx[i] = static_cast<int>(send_list[i]);
  for (unsigned int i = 0; i != this->local_size(); ++i)   idx[n_sl+i] = i + this->first_local_index();
  
  // Create the index set & scatter object
  if (idx.empty())  ierr = ISCreateGeneral(MPI_COMM_WORLD,n_sl+this->local_size(), PETSC_NULL, PETSC_USE_POINTER, &is);
  else    ierr = ISCreateGeneral(MPI_COMM_WORLD,n_sl+this->local_size(), &idx[0], PETSC_USE_POINTER, &is);
           CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterCreate(_vec,is,v_local->_vec, is,&scatter);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Perform the scatter
// #if PETSC_VERSION_LESS_THAN(2,3,3)
// 	 
//   ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
// 			 SCATTER_FORWARD, scatter);
//          CHKERRABORT(MPI_COMM_WORLD,ierr);
//   
//   ierr = VecScatterEnd  (_vec, v_local->_vec, INSERT_VALUES,
// 			 SCATTER_FORWARD, scatter);
//          CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
// #else
  // API argument order change in PETSc 2.3.3
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,INSERT_VALUES, SCATTER_FORWARD);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  ierr = VecScatterEnd  (scatter, _vec, v_local->_vec,INSERT_VALUES, SCATTER_FORWARD);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
	
  // Clean up
  ierr = ISDestroy (&is);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Make sure ghost dofs are up to date
  if (v_local->type() == GHOSTEDM)  v_local->close();
}

// ================================================================
void PetscVectorM::localize (const unsigned int first_local_idx,
			       const unsigned int last_local_idx,
			       const std::vector<unsigned int>& send_list){
  this->_restore_array();
  // Only good for serial vectors.
  // assert (this->size() == this->local_size());
  assert (last_local_idx > first_local_idx); assert (send_list.size() <= this->size());
  assert (last_local_idx < this->size());
  
  const unsigned int size       = this->size();
  const unsigned int local_size = (last_local_idx - first_local_idx + 1);
  int ierr=0;  
  
  // Don't bother for serial cases
  if ((first_local_idx == 0) &&      (local_size == size))   return;
  
  // Build a parallel vector, initialize it with the local parts of (*this)
  PetscVectorM parallel_vec; parallel_vec.init (size, local_size, true, PARALLELM);

  // Copy part of *this into the parallel_vec
  {
    IS is;   VecScatter scatter;

    // Create idx, idx[i] = i+first_local_idx;
    std::vector<int> idx(local_size);
    iota (idx.begin(), idx.end(), first_local_idx);  // Utility::iota (idx.begin(), idx.end(), first_local_idx);

    // Create the index set & scatter object
    ierr = ISCreateGeneral(MPI_COMM_WORLD, local_size, &idx[0], PETSC_USE_POINTER, &is);CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterCreate(_vec,is,parallel_vec._vec, is, &scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);

    // Perform the scatter
// #if PETSC_VERSION_LESS_THAN(2,3,3)
// 
//     ierr = VecScatterBegin(_vec, parallel_vec._vec, INSERT_VALUES,
// 			   SCATTER_FORWARD, scatter);
//            CHKERRABORT(MPI_COMM_WORLD,ierr);
//   
//     ierr = VecScatterEnd  (_vec, parallel_vec._vec, INSERT_VALUES,
// 			   SCATTER_FORWARD, scatter);
//            CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
// #else
	   
      // API argument order change in PETSc 2.3.3
    ierr = VecScatterBegin(scatter, _vec, parallel_vec._vec,INSERT_VALUES, SCATTER_FORWARD);
           CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterEnd  (scatter, _vec, parallel_vec._vec,INSERT_VALUES, SCATTER_FORWARD);
           CHKERRABORT(MPI_COMM_WORLD,ierr);
	   
// #endif
    // Clean up
    ierr = ISDestroy (&is); CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);
  }

  // localize like normal
  parallel_vec.close();  parallel_vec.localize (*this, send_list);
  this->close();
}



// =============================================================
void PetscVectorM::localize (std::vector<Real>& v_local) const{
  this->_restore_array();

  // This function must be run on all processors at once
  parallel_onlyM();

  int ierr=0; PetscScalar *values;
  const int n = this->size();  const int nl = this->local_size();
  v_local.clear();  v_local.resize(n, 0.);

  ierr = VecGetArray (_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);

  unsigned int ioff = first_local_index();
  for (int i=0; i<nl; i++)   v_local[i+ioff] = static_cast<Real>(values[i]);
  ierr = VecRestoreArray (_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
  ParallelM::sum(v_local);   //TODO must become ParallelM
}

//  ===========================================================
void PetscVectorM::localize_to_one (std::vector<Real>& v_local,
					 const unsigned int pid) const{
  this->_restore_array();
  int ierr=0; PetscScalar *values;
  const int n  = size(); const int nl = local_size();
 
  v_local.resize(n);
  // only one processor
  if (n == nl)    {      
      ierr = VecGetArray (_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
      for (int i=0; i<n; i++)v_local[i] = static_cast<Real>(values[i]);
      ierr = VecRestoreArray (_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  // otherwise multiple processors
  else    {
      unsigned int ioff = this->first_local_index();
      std::vector<Real> local_values (n, 0.);
      {
	ierr = VecGetArray (_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	for (int i=0; i<nl; i++) local_values[i+ioff] = static_cast<Real>(values[i]);
	ierr = VecRestoreArray (_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
      }
      MPI_Reduce (&local_values[0], &v_local[0],n,MPI_REAL,MPI_SUM,pid,MPI_COMM_WORLD);
    }
}

// =======================================================
void PetscVectorM::pointwise_mult (const NumericVectorM& vec1,
				     const NumericVectorM& vec2){
  this->_restore_array();

  int ierr = 0;

  // Convert arguments to PetscVector*.
  const PetscVectorM* vec1_petsc = static_cast<const PetscVectorM*>(&vec1);
  const PetscVectorM* vec2_petsc = static_cast<const PetscVectorM*>(&vec2);

  // Call PETSc function.

#if PETSC_VERSION_LESS_THAN(2,3,1)

  std::cout << "This method has been developed with PETSc 2.3.1.  "
	    << "No one has made it backwards compatible with older "
	    << "versions of PETSc so far; however, it might work "
	    << "without any change with some older version." << std::endl;
  libmesh_error();

#else
  
  if(this->type() != GHOSTEDM)
    {
      ierr = VecPointwiseMult(this->vec(),
			      const_cast<PetscVectorM*>(vec1_petsc)->vec(),
			      const_cast<PetscVectorM*>(vec2_petsc)->vec());
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  else
    {
	  Vec loc_vec;
	  Vec v1_loc_vec;
	  Vec v2_loc_vec;
	  ierr = VecGhostGetLocalForm (_vec,&loc_vec);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecGhostGetLocalForm (const_cast<PetscVectorM*>(vec1_petsc)->vec(),&v1_loc_vec);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecGhostGetLocalForm (const_cast<PetscVectorM*>(vec2_petsc)->vec(),&v2_loc_vec);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);

	  ierr = VecPointwiseMult(loc_vec,v1_loc_vec,v2_loc_vec);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);

	  ierr = VecGhostRestoreLocalForm (const_cast<PetscVectorM*>(vec1_petsc)->vec(),&v1_loc_vec);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecGhostRestoreLocalForm (const_cast<PetscVectorM*>(vec2_petsc)->vec(),&v2_loc_vec);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
      
#endif

}


// ======================================================

/// Print the contents of the vector to hdf5 file,
void PetscVectorM::print_hdf5(std::string file) const {
   this->_restore_array(); assert (this->closed());
  uint ndim=this->size();
  std::ostringstream name;
  name.str("");    name << file.c_str() << ".h5";
  hid_t fileP = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  // Value storage
//   double *val=new double[ndim];  for (uint i=0;i<ndim;i++)  val[i]= V__GetCmp(&_vec,i+1);

// print matrix values
  hsize_t dimsf[2]; dimsf[0]=1;    dimsf[1] = 1;
  int *nd=new int[1];nd[0]=ndim;
  name.str(""); name << "DIM" ;
  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,nd);
  H5Sclose(dataspace);  H5Dclose(dataset);

  dimsf[0]=ndim;
  name.str(""); name << "VAL" ;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_DOUBLE,
                      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, &(_vec));
  H5Sclose(dataspace);  H5Dclose(dataset);

 
  H5Fclose(fileP);
 /* delete []val; */delete []nd;
  return;
  
// void PetscVectorM::print_hdf5(const std::string name) const{
//   this->_restore_array(); assert (this->closed());
//   
//   int ierr=0; 
//   PetscViewer petsc_viewer;
//   ierr = PetscViewerCreate (MPI_COMM_WORLD, &petsc_viewer);
//          CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//   // Create an ASCII file containing the matrix if a filename was provided.  
//   if (name != "NULL"){
//       ierr = PetscViewerASCIIOpen( MPI_COMM_WORLD,
// 				   name.c_str(),&petsc_viewer);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       ierr = PetscViewerSetFormat (petsc_viewer,
// 				   PETSC_VIEWER_ASCII_MATLAB);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       ierr = VecView (_vec, petsc_viewer);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//     }
// 
//   // Otherwise the matrix will be dumped to the screen.
//   else {
//       ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
// 				   PETSC_VIEWER_ASCII_MATLAB);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       ierr = VecView (_vec, PETSC_VIEWER_STDOUT_WORLD);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//     }
//   // Destroy the viewer.
//   ierr = PetscViewerDestroy (petsc_viewer);
//          CHKERRABORT(MPI_COMM_WORLD,ierr);
 }

// ========================================================
void PetscVectorM::print_personal(std::ostream& /*name*/) const{
  this->_restore_array();
  assert (this->closed());
  
  int ierr=0; 
  PetscViewer petsc_viewer;
  ierr = PetscViewerCreate(MPI_COMM_WORLD,&petsc_viewer);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
				   PETSC_VIEWER_ASCII_MATLAB);
             CHKERRABORT(MPI_COMM_WORLD,ierr);
  
      ierr = VecView (_vec, PETSC_VIEWER_STDOUT_WORLD);
             CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Destroy the viewer.
  ierr = PetscViewerDestroy (&petsc_viewer);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
}



// ==============================================================
void PetscVectorM::create_subvector(NumericVectorM& subvector,
				      const std::vector<unsigned int>& rows) const{
  this->_restore_array();

  // PETSc data structures
  IS parent_is, subvector_is;
  VecScatter scatter;
  int ierr = 0;
  
  // Make sure the passed in subvector is really a PetscVector
  PetscVectorM* petsc_subvector = static_cast<PetscVectorM*>(&subvector);
  
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
			  &(petsc_subvector->_vec)); CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecSetFromOptions (petsc_subvector->_vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
      // Mark the subvector as initialized
      petsc_subvector->_is_initialized = true;
    }
  else {
      petsc_subvector->_restore_array();
    }
  
  // Use iota to fill an array with entries [0,1,2,3,4,...rows.size()]
  std::vector<int> idx(rows.size());
  iota (idx.begin(), idx.end(), 0);
//   Utility::iota (idx.begin(), idx.end(), 0);

  // Construct index sets
  ierr = ISCreateGeneral(MPI_COMM_WORLD,
			 rows.size(),
			 (int*) &rows[0],
			  PETSC_USE_POINTER,
			 &parent_is); CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = ISCreateGeneral(MPI_COMM_WORLD,
			 rows.size(),
			 (int*) &idx[0],
			  PETSC_USE_POINTER,
			 &subvector_is); CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Construct the scatter object
  ierr = VecScatterCreate(this->_vec,
			  parent_is,
			  petsc_subvector->_vec,
			  subvector_is,
			  &scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Actually perform the scatter
#if PETSC_VERSION_LESS_THAN(2,3,3)
  ierr = VecScatterBegin(this->_vec,
			 petsc_subvector->_vec,
			 INSERT_VALUES,
			 SCATTER_FORWARD,
			 scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = VecScatterEnd(this->_vec,
		       petsc_subvector->_vec,
		       INSERT_VALUES,
		       SCATTER_FORWARD,
		       scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);
#else
  // API argument order change in PETSc 2.3.3
  ierr = VecScatterBegin(scatter,
			 this->_vec,
			 petsc_subvector->_vec,
			 INSERT_VALUES,
			 SCATTER_FORWARD); CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = VecScatterEnd(scatter,
		       this->_vec,
		       petsc_subvector->_vec,
		       INSERT_VALUES,
		       SCATTER_FORWARD); CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
  
  // Clean up 
  ierr = ISDestroy(&parent_is);       CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = ISDestroy(&subvector_is);    CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter); CHKERRABORT(MPI_COMM_WORLD,ierr); 

}




//------------------------------------------------------------------
// Explicit instantiations
// template class PetscVector<Number>;



#endif // #ifdef LIBMESH_HAVE_PETSC
