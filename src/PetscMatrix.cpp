#include "FEMTTUConfig.h"

#if HAVE_PETSC == 1

// Local includes
#include "Parallel.hpp"
#include "PetscMatrix.hpp"
#include "PetscVector.hpp"
#include "DenseMatrix.hpp"
#include <mpi.h>
#include <hdf5.h>
#include <sstream>
using namespace std;

//-----------------------------------------------------------------------
// PetscMatrix members
void PetscMatrix::init(const  int m,
                        const  int n,
                        const  int m_l,
                        const  int n_l,
                        const  int nnz,
                        const  int noz) {
  // Set matrix dimension
  _m=m;
  _n=n;
  _m_l=m_l;
  _n_l=n_l;

// _len_row.resize(_m_l);
// _len_diag_row.resize(_m_l);
  // We allow 0x0 matrices now
  //if ((m==0) || (n==0))
  //  return;

  // Clear initialized matrices
  if (this->initialized())  this->clear();
  this->_is_initialized = true;

  // processor info
  int proc_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  int ierr     = 0;
  int m_global = static_cast<int>(m);
  int n_global = static_cast<int>(n);
  int m_local  = static_cast<int>(m_l);
  int n_local  = static_cast<int>(n_l);
  int n_nz     = static_cast<int>(nnz);
  int n_oz     = static_cast<int>(noz);

// create a sequential matrix on one processor
  if (numprocs == 1) {
    assert ((m_l == m) && (n_l == n));

    // Create matrix.  Revisit later to do preallocation and make more efficient
    ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global,
                            n_nz, PETSC_NULL, &_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = MatSetFromOptions (_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }

  else {
    parallel_onlyM();
 
    ierr = MatCreate(MPI_COMM_WORLD, &_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = MatSetSizes(_mat, m_l, n_l, m, n);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = MatSetType(_mat, MATMPIAIJ); // Automatically chooses seqaij or mpiaij
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = MatMPIAIJSetPreallocation(_mat, nnz, PETSC_NULL, noz, PETSC_NULL);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    
  }
  
  this->zero ();
}

// =====================================0
void PetscMatrix::update_sparsity_pattern(
  int m_global,                          // # global rows
  int n_global,                          // # global columns
  int m_local,                           // # local rows (local proc)
  int n_local,                           // # local columns (local proc)
  const std::vector< int>  n_nz, // # diagoanl entries
  const std::vector< int>  n_oz  // # offset entries
) {

  // Clear initialized matrices
  if (this->initialized())    this->clear();
  this->_is_initialized = true;

  // processor info
  int proc_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  int ierr     = 0;

  // create a sequential matrix on one processor -------------------
  if (numprocs == 1)    {
    assert((m_local == m_global) && (n_local == n_global));
    if (n_nz.empty())
      ierr = MatCreateSeqAIJ(MPI_COMM_WORLD, m_global, n_global,
                             PETSC_DEFAULT, (int*) PETSC_NULL, &_mat);
    else
      ierr = MatCreateSeqAIJ(MPI_COMM_WORLD, m_global, n_global,
                             PETSC_DEFAULT, (int*) &n_nz[0], &_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = MatSetFromOptions(_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else    { // multi processors ----------------------------
    parallel_onlyM();             //TODO
    if (n_nz.empty()) {

//old piece of code with Petsc lib version 3.1
//       ierr = MatCreateMPIAIJ(MPI_COMM_WORLD,
//                              m_local, n_local,
//                              m_global, n_global,
//                              PETSC_NULL, (int*) PETSC_NULL,
//                              PETSC_NULL, (int*) PETSC_NULL,
//                              &_mat
//                             );

      ierr = MatCreate(MPI_COMM_WORLD, &_mat);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = MatSetType(_mat, MATMPIAIJ);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = MatMPIAIJSetPreallocation(_mat, 0, 0, 0, 0);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

    }

    else {

//old piece of code with Petsc lib version 3.1
//       ierr = MatCreateMPIAIJ(MPI_COMM_WORLD,
//                              m_local, n_local,
//                              m_global, n_global,
//                              PETSC_NULL, (int*) &n_nz[0],
//                              PETSC_NULL, (int*) &n_oz[0],
//                              &_mat
//                             );

      ierr = MatCreate(MPI_COMM_WORLD, &_mat);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = MatSetType(_mat, MATMPIAIJ); // Automatically chooses seqaij or mpiaij
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = MatMPIAIJSetPreallocation(_mat, 0, (int*)&n_nz[0], 0, (int*)&n_oz[0]);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

    }

    ierr = MatSetFromOptions(_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } // ----------------------------------------------------

  this->zero();
  return;
}

// =====================================0
// void PetscMatrix::update_sparsity_pattern(const Graph &sparsity_pattern) {
//
//   // dimension
//   const  int m   =  sparsity_pattern._m; // global rows;
//   const  int n   =  sparsity_pattern._n; // global columns
//   const  int n_l =sparsity_pattern._nl;  // local rows
//   const  int m_l =sparsity_pattern._ml;  // local columns
//   const  int ml_start =sparsity_pattern._ml_start;
//
//   // vectors n_nz (diagonal) and n_oz (offset) -------------------------------
//   std::vector< int> n_nz(m_l);   std::vector< int> n_oz(m_l);
//   for (int i=0;i<m_l;i++) {
//     int len=sparsity_pattern[ml_start+i].size()-1;
//     n_oz[i]=sparsity_pattern[ml_start+i][len];
//     n_nz[i]=len-n_oz[i];
//   }
//   update_sparsity_pattern(m,n,m_l,n_l,n_oz,n_nz);
//   return;
// }

// ==========================================
/// This function sets all the entries to 0
void PetscMatrix::zero() {
  assert(this->initialized());
  semiparallel_onlyM();
  int ierr=0;
  ierr = MatZeroEntries(_mat);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// =========================================================================
/// This function sets all row entries to 0 then puts diag_value in the diagonal entry
void PetscMatrix::zero_rows(std::vector<int> & rows, double diag_value) {
  assert(this->initialized());
  semiparallel_onlyM();
  int ierr=0;
  if (!rows.empty()) ierr = MatZeroRows(_mat, rows.size(), &rows[0], diag_value,0,0);   // add,0,0,)    !!!!
  else   ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value,0,0);                       // add,0,0,)    !!!!
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// =================================================
void PetscMatrix::clear() {
  int ierr=0;
  if ((this->initialized()) && (this->_destroy_mat_on_exit)) {
    semiparallel_onlyM();
    ierr = MatDestroy(&_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    this->_is_initialized = false;
  }
}

// ============================================
double PetscMatrix::l1_norm() const {
  assert(this->initialized());
  semiparallel_onlyM();
  int ierr=0;
  PetscReal petsc_value;
  double value;
  assert(this->closed());
  ierr = MatNorm(_mat, NORM_1, &petsc_value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  value = static_cast<double>(petsc_value);
  return value;
}

// ==========================================
double PetscMatrix::linfty_norm() const {
  assert(this->initialized());
  semiparallel_onlyM();
  int ierr=0;
  PetscReal petsc_value;
  double value;
  assert(this->closed());
  ierr = MatNorm(_mat, NORM_INFINITY, &petsc_value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  value = static_cast<double>(petsc_value);
  return value;
}

// ===================================================
void PetscMatrix::print_personal(
  std::ostream& os // pointer stream
) const {// =========================================
  assert(this->initialized());
  std::cout << "\n PetscMatrix::print_personal \n";
#ifndef NDEBUG
  if (os != std::cout)
    std::cerr << "Warning! PETSc can only print to std::cout!" << std::endl;
#endif
  int ierr=0;
  ierr = MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}
// =====================================================
void PetscMatrix::print_hdf5(const std::string /*name*/) const {




//   hid_t status =0;
//
// //   std::ostringstream name;
// //   name << femus_dir << "/" << input_dir << f_matrix  << Level1 << ext_h5;
//   // print quadratic  --------------------------------------------------
// //    hid_t fileM= H5Fopen(name.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
//    hid_t fileM = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
// //   hid_t fileM = H5Fopen(name.c_str(),  H5P_DEFAULT,H5P_DEFAULT);
//   // matrix dimensions
//   hsize_t dimsf[2];     dimsf[0]=2;  dimsf[1] = 1;
//    int rowcln[2];   rowcln[0]=_n_l;    rowcln[1]=_m_l;
//
//
//   hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
//   hid_t dataset = H5Dcreate(fileM,"DIM",H5T_NATIVE_INT,
//                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,rowcln);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//
//
// //    status =_mgutils.print_Ihdf5(fileM,"DIM", dimsf,rowcln);
//   // print matrix position of P
// //    std::ostringstream name1;name1 << "POS" << qq;
// //   dimsf[0]=count_q;  status =_mgutils.print_Ihdf5(fileM,name1.str().c_str(), dimsf,Mat_q);
//   // print row length of R
// //   std::ostringstream name2;name2 << "LEN";
//   dimsf[0]=_m_l;//    status =_mgutils.print_Ihdf5(fileM,name2.str().c_str(), dimsf,_len_row[0]);
//    dataspace = H5Screate_simple(2,dimsf, NULL);
//    dataset = H5Dcreate(fileM,"LEN",H5T_NATIVE_INT,
//                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&_len_row[0]);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//
//   // print row length of R
// //   std::ostringstream name3;name3 << "OFFLEN";
// //   status =_mgutils.print_Ihdf5(fileM,name3.str().c_str(), dimsf,_len_diag_row[0]);
//     dataspace = H5Screate_simple(2,dimsf, NULL);
//    dataset = H5Dcreate(fileM,"DIAGLEN",H5T_NATIVE_INT,
//                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&_len_diag_row[0]);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//   //   clean
//    H5Fclose(fileM);



// std::cout << " PetscMatrix::print_hdf5 in file   " <<  name   <<   "\n";
// #endif-
}


// =====================================================
// void PetscMatrix::read_hdf5(const std::string namefile,const int mode,const int ml_init) {
//
//
//    hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
//
//   // quadratic-quadratic  ------------------------------------------------
//   std::ostringstream name_dst;  name_dst << "DIM" <<mode;
//   int rowcln_q[2];hid_t dataset=H5Dopen(file_id,name_dst.str().c_str(), H5P_DEFAULT);
//   hid_t status=H5Dread(dataset,
//                        H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,rowcln_q);
//   int n_row_q=rowcln_q[0];
//   // matrix offset (subdomainsand levelels)
// //   int *off_nd_q= _mgmesh._off_nd;
//   // row legth
//   name_dst.str(""); name_dst << "LEN" <<mode;
//   uint *length_row=new uint[n_row_q+1];
//   dataset=H5Dopen(file_id,name_dst.str().c_str(), H5P_DEFAULT);
//   status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,length_row);
//   // matrix off diagonal
//   name_dst.str(""); name_dst <<"OFFLEN" <<mode;
//   uint *length_offrow=new uint[n_row_q+1];
//   dataset=H5Dopen(file_id,name_dst.str().c_str(), H5P_DEFAULT);
//   status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,length_offrow);
// //   // matrix off diagonal
// //   uint *pos_row=new uint[length_row[n_row_q]];dataset=H5Dopen(file_id,"POS2", H5P_DEFAULT);
// //   status=H5Dread(dataset,
// //                  H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,pos_row);
//
//     for ( int i=ml_init; i<ml_init+_m_l; i++) {
//       // row index
//       _len_row[i-ml_init]= length_row[i+1]-length_row[i];
//       int  loffrow= length_offrow[i+1]-length_offrow[i];
//       _len_diag_row[i-ml_init]=_len_row[i-ml_init]-loffrow;
//     }
//
//   delete   []length_row;
//   delete   []length_offrow;
//   H5Fclose(file_id);
//   return;
//
// }

// // ============================================================
void PetscMatrix::add_matrix(const DenseMatrix& dm,
                              const std::vector< int>& rows,
                              const std::vector< int>& cols) {
  assert(this->initialized());
  const  int m = dm.m();
  assert((int)rows.size() == m);
  const  int n = dm.n();
  assert((int)cols.size() == n);

  int ierr=0;
  // These casts are required for PETSc <= 2.1.5
  ierr = MatSetValues(_mat,
                      m, (int*) &rows[0],
                      n, (int*) &cols[0],
                      (PetscScalar*) &dm.get_values()[0],
                      ADD_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// ============================================================
/// This Petsc adds a dense matrix to a sparse matrix
void PetscMatrix::add_matrix_blocked(
  const std::vector< double > &mat_values,  // blocked matrix stored as an m*n array 
  const std::vector< int>& rows, // row vector indexes
  const std::vector< int>& cols) { // column vector indices
  // ==================================================
  assert (this->initialized());
  int ierr=0;
//   const  int m = dm.m();  //senza DenseMatrix non puoi fare questo controllo, oppure aggiungi altri argomenti
//   const  int n = dm.n();
//   assert ((int)rows.size() == m);
//   assert ((int)cols.size() == n);
  // These casts are required for PETSc <= 2.1.5
//   ierr = MatSetValuesBlocked(_mat,/*m*/rows.size(), (int*) &rows[0],
// 			          /*n*/cols.size(), (int*) &cols[0],
//                       (PetscScalar*) mat_ptr/*&dm.get_values()[0]*/,ADD_VALUES);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);  
  const  int m = (int)rows.size(); 
  const  int n = (int)cols.size();
  assert ( m*n == mat_values.size());
  
  //These casts are required for PETSc <= 2.1.5
  ierr = MatSetValuesBlocked(_mat,m, &rows[0],n, &cols[0],
                      (PetscScalar*) &mat_values[0],ADD_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  return;
}



// // ============================================================

void PetscMatrix::matrix_PtAP(const SparseMatrix &mat_P, const SparseMatrix &mat_A, const bool &mat_reuse){
  
  const PetscMatrix* A = static_cast<const PetscMatrix*>(&mat_A);
  A->close();
  
  const PetscMatrix* P = static_cast<const PetscMatrix*>(&mat_P);
  P->close();
  
  int ierr=0;
  if(mat_reuse){  
    ierr = MatPtAP(const_cast<PetscMatrix*>(A)->mat(), const_cast<PetscMatrix*>(P)->mat(), MAT_REUSE_MATRIX,1.0,&_mat);
  }
  else{
    this->clear();
    ierr = MatPtAP(const_cast<PetscMatrix*>(A)->mat(), const_cast<PetscMatrix*>(P)->mat(), MAT_INITIAL_MATRIX ,1.0,&_mat);
    this->_is_initialized = true;
  }
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// ===========================================================
/// This function either creates or re-initializes a matrix called "submatrix".
void PetscMatrix::_get_submatrix(SparseMatrix& submatrix,
                                  const std::vector< int> &rows,
                                  const std::vector< int> &cols,
                                  const bool reuse_submatrix) const {

  // Can only extract submatrices from closed matrices
  this->close();

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix* petsc_submatrix = libmeshM_cast_ptr<PetscMatrix*>(&submatrix);

  // If we're not reusing submatrix and submatrix is already initialized
  // then we need to clear it, otherwise we get a memory leak.
  if (!reuse_submatrix && submatrix.initialized())  submatrix.clear();
  // Construct row and column index sets.
  int ierr=0;
  IS isrow, iscol;

  ierr = ISCreateGeneral(MPI_COMM_WORLD,rows.size(),(int*) &rows[0],PETSC_COPY_VALUES, &isrow);    // PETSC_COPY_VALUES is my first choice; see also PETSC_OWN_POINTER, PETSC_USE_POINTER
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = ISCreateGeneral(MPI_COMM_WORLD,cols.size(),(int*) &cols[0],PETSC_COPY_VALUES, &iscol);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

//---

#if !PETSC_VERSION_LESS_THAN(3,0,1) || !PETSC_VERSION_RELEASE
  ierr = MatGetSubMatrix(_mat,
                         isrow,
                         iscol,
                         (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                         &(petsc_submatrix->_mat));
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#else
  ierr = MatGetSubMatrix(_mat,
                         isrow,
                         iscol,
                         PETSC_DECIDE,
                         (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                         &(petsc_submatrix->_mat));
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif

  // Specify that the new submatrix is initialized and close it.
  petsc_submatrix->_is_initialized = true;
  petsc_submatrix->close();

  // Clean up PETSc data structures
  ierr = ISDestroy(&isrow);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = ISDestroy(&iscol);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  return;
}

// ======================================================
void PetscMatrix::get_diagonal(NumericVector& dest) const {
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector& petsc_dest = libmeshM_cast_ref<PetscVector&>(dest);
  // Call PETSc function.
  // Needs a const_cast since PETSc does not work with const.
  int ierr =  MatGetDiagonal(const_cast<PetscMatrix*>(this)->mat(),petsc_dest.vec());
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return;
}

/// This function copies the transpose of the matrix
// ===========================================================
void PetscMatrix::get_transpose(
  SparseMatrix& dest  // output transpose matrix
) const {

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix& petsc_dest = libmeshM_cast_ref<PetscMatrix&>(dest);

  // If we aren't reusing the matrix then need to clear dest,
  // otherwise we get a memory leak
  if (&petsc_dest != this)    dest.clear();

  int ierr;
#if PETSC_VERSION_LESS_THAN(3,0,0)
  if (&petsc_dest == this)
    ierr = MatTranspose(_mat,PETSC_NULL);
  else
    ierr = MatTranspose(_mat,&petsc_dest._mat);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#else
  // FIXME - we can probably use MAT_REUSE_MATRIX in more situations
  if (&petsc_dest == this)
    ierr = MatTranspose(_mat,MAT_REUSE_MATRIX,&petsc_dest._mat);
  else
    ierr = MatTranspose(_mat,MAT_INITIAL_MATRIX,&petsc_dest._mat);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
  // Specify that the transposed matrix is initialized and close it.
  petsc_dest._is_initialized = true;
  petsc_dest.close();
}

#endif // #ifdef LIBMESH_HAVE_PETSC
