#include "FEMTTUConfig.h"

#if HAVE_PETSC == 1

// Local includes
#include "PetscMacro.hpp"
#include "PetscRectangularMatrix.hpp"
#include "DenseMatrix.hpp"
#include "PetscVector.hpp"
#include <sstream>
#include "hdf5.h"
using namespace std;


//-----------------------------------------------------------------------
// PetscMatrix members
void PetscRectangularMatrix::init (const  int m,
                          const  int n,
                          const  int m_l,
                          const  int n_l,
                          const  int nnz,
                          const  int noz) {
  _m=m;
  _n=n;
  _m_l=m_l;
  _n_l=n_l;

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

//   int m_global = static_cast<int>(m);
//   int n_global = static_cast<int>(n);
//   int m_local  = static_cast<int>(m_l);
//   int n_local  = static_cast<int>(n_l);
  //PetscInt n_nz = static_cast<PetscInt>(nnz);
  //PetscInt n_oz = static_cast<PetscInt>(noz);

  //create a sequential matrix on one processor
  if (numprocs == 1) {
    assert ((m_l == m) && (n_l == n));

    // Create matrix.  Revisit later to do preallocation and make more efficient
   ierr = MatCreateSeqAIJ (PETSC_COMM_SELF, m, n,nnz, PETSC_NULL, &_mat);
   CHKERRABORT(PETSC_COMM_SELF,ierr);
   ierr = MatSetFromOptions (_mat);
   CHKERRABORT(PETSC_COMM_SELF,ierr);
 
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
  this->zero();
}

// =====================================0
void PetscRectangularMatrix::init () {
  std::cout << "  PetscRectangularMatrix::init: Our PetscMatrix  needs dimensions";
  abort();
}
// =====================================0

// =====================================0
void PetscRectangularMatrix::update_sparsity_pattern(
  int m_global,                          // # global rows
  int n_global,                          // # global columns
  int m_local,                           // # local rows (local proc)
  int n_local,                           // # local columns (local proc)
  const std::vector< int>  n_oz, // # diagoanl entries
  const std::vector< int>  n_nz  // # offset entries
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
  if (numprocs == 1) {
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
  } else { // multi processors ----------------------------
    parallel_onlyM();     //TODO
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

    }   else {

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



// ===================================
void PetscRectangularMatrix::zero () {
  assert (this->initialized());
  semiparallel_onlyM();
  int ierr=0;
  ierr = MatZeroEntries(_mat);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// =======================================================
/// This function sets to zero the desired rows:

void PetscRectangularMatrix::zero_rows (std::vector<int> & rows, double diag_value) {
/// MatZeroRows now takes two additional optional arguments.
///  The optional arguments (x,b) can be used to specify the
/// solutions for the zeroed rows (x) and right hand side (b) to update.
/// Could be useful for setting boundary conditions
// =======================================================

  assert(this->initialized());
  semiparallel_onlyM();
  int ierr=0;

#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,1,1)
  if (!rows.empty()) ierr = MatZeroRows(_mat, rows.size(), &rows[0], diag_value);
  else  ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value);
#else

  if (!rows.empty())  ierr = MatZeroRows(_mat, rows.size(), &rows[0], diag_value, PETSC_NULL, PETSC_NULL);
  else  ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value, PETSC_NULL, PETSC_NULL);
#endif
  CHKERRABORT(MPI_COMM_WORLD,ierr);
//     assert (this->initialized());
//     semiparallel_onlyM();
//     int ierr=0;
//     if (!rows.empty())
//         ierr = MatZeroRows(_mat, rows.size(), &rows[0], diag_value);
//     else
//         ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
  return;
}



void PetscRectangularMatrix::matrix_PtAP(const SparseRectangularMatrix &mat_P, const SparseRectangularMatrix &mat_A, const bool &mat_reuse){
  
  const PetscRectangularMatrix* A = static_cast<const PetscRectangularMatrix*>(&mat_A);
  A->close();
  
  const PetscRectangularMatrix* P = static_cast<const PetscRectangularMatrix*>(&mat_P);
  P->close();
  
  int ierr=0;
  if(mat_reuse){  
    ierr = MatPtAP(const_cast<PetscRectangularMatrix*>(A)->mat(), const_cast<PetscRectangularMatrix*>(P)->mat(), MAT_REUSE_MATRIX,1.0,&_mat);
  }
  else{
    this->clear();
    ierr = MatPtAP(const_cast<PetscRectangularMatrix*>(A)->mat(), const_cast<PetscRectangularMatrix*>(P)->mat(), MAT_INITIAL_MATRIX ,1.0,&_mat);
    this->_is_initialized = true;
  }
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}




// =================================================
/// This function clears the matrix
void PetscRectangularMatrix::clear () {
  // ============================================
  int ierr=0;
  if ((this->initialized()) && (this->_destroy_mat_on_exit)) {
    //semiparallel_onlyM();
    ierr = MatDestroy (&_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    this->_is_initialized = false;
  }
}

// ============================================
double PetscRectangularMatrix::l1_norm () const {
  assert (this->initialized());

  semiparallel_onlyM();
  int ierr=0;
  PetscReal petsc_value;
  double value;
  assert (this->closed());
  ierr = MatNorm(_mat, NORM_1, &petsc_value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  value = static_cast<double>(petsc_value);
  return value;
}

// ==========================================
double PetscRectangularMatrix::linfty_norm () const {
  assert (this->initialized());
  semiparallel_onlyM();
  int ierr=0;
  PetscReal petsc_value;
  double value;
  assert (this->closed());
  ierr = MatNorm(_mat, NORM_INFINITY, &petsc_value);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  value = static_cast<double>(petsc_value);
  return value;
}



// =====================================================
void PetscRectangularMatrix::print_matlab (const std::string name) const {
  assert (this->initialized());
  semiparallel_onlyM();
  // assert (this->closed());
  this->close();
  int ierr=0;
  PetscViewer petsc_viewer;
  ierr = PetscViewerCreate (MPI_COMM_WORLD, &petsc_viewer);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.
   */
  if (name != "NULL")  {
    ierr = PetscViewerASCIIOpen( MPI_COMM_WORLD,name.c_str(),&petsc_viewer);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = PetscViewerSetFormat (petsc_viewer,PETSC_VIEWER_ASCII_MATLAB);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = MatView (_mat, petsc_viewer);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }

  // Otherwise the matrix will be dumped to the screen.
  else {
    ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  // Destroy the viewer.
  ierr = PetscViewerDestroy (&petsc_viewer);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// ============================================================
/// This Petsc adds a dense matrix to a sparse matrix
void PetscRectangularMatrix::add_matrix(
  const DenseMatrix& dm,        // dense matrix
  const std::vector< int>& rows, // row vector
  const std::vector< int>& cols) { // column vector
  // ==================================================
  assert (this->initialized());
  int ierr=0;
  const  int m = dm.m();
  const  int n = dm.n();
  assert ((int)rows.size() == m);
  assert ((int)cols.size() == n);
  // These casts are required for PETSc <= 2.1.5
  ierr = MatSetValues(_mat,m, (int*) &rows[0],n, (int*) &cols[0],
                      (PetscScalar*) &dm.get_values()[0],ADD_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return;
}

// ===========================================================
/// This Petsc function exctract submatrix from matrix
void PetscRectangularMatrix::_get_submatrix(
  SparseRectangularMatrix& submatrix,       // submatrix
  const std::vector< int> &rows,   // row vector
  const std::vector< int> &cols,   // column vector
  const bool reuse_submatrix       // flag submatrix
) const {
  // =====================================================
  this->close();  // Can only extract submatrices from closed matrices

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscRectangularMatrix* petsc_submatrix = libmeshM_cast_ptr<PetscRectangularMatrix*>(&submatrix);

  // If we're not reusing submatrix and submatrix is already initialized
  // then we need to clear it, otherwise we get a memory leak.
  if ( !reuse_submatrix && submatrix.initialized() )  submatrix.clear();
  // Construct row and column index sets.
  int ierr=0;
  IS isrow, iscol;
  ierr = ISCreateGeneral(MPI_COMM_WORLD,rows.size(),
                         (int*) &rows[0],PETSC_USE_POINTER,&isrow);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = ISCreateGeneral(MPI_COMM_WORLD,cols.size(),
                         (int*) &cols[0],PETSC_USE_POINTER,&iscol);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Extract submatrix
#if !PETSC_VERSION_LESS_THAN(3,0,1) || !PETSC_VERSION_RELEASE
  ierr = MatGetSubMatrix(_mat, isrow, iscol,
//                            PETSC_DECIDE,   //PETSC version 3.1
                         (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                         &(petsc_submatrix->_mat));
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#else
  ierr = MatGetSubMatrix(_mat,isrow, iscol,
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
void PetscRectangularMatrix::get_diagonal (NumericVector& dest) const {
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector& petsc_dest = libmeshM_cast_ref<PetscVector&>(dest);
  // Call PETSc function.
#if PETSC_VERSION_LESS_THAN(2,3,1)

  std::cout << "This method has been developed with PETSc 2.3.1.  "
            << "No one has made it backwards compatible with older "
            << "versions of PETSc so far; however, it might work "
            << "without any change with some older version." << std::endl;
  libmesh_error();
#else
  // Needs a const_cast since PETSc does not work with const.
  int ierr =  MatGetDiagonal(const_cast<PetscRectangularMatrix*>(this)->mat(),petsc_dest.vec());
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
}


// ===========================================================
void PetscRectangularMatrix::get_transpose (SparseRectangularMatrix& dest) const {
  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscRectangularMatrix& petsc_dest = libmeshM_cast_ref<PetscRectangularMatrix&>(dest);

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
// =====================================================
void PetscRectangularMatrix::print_hdf5(const std::string name) const {




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
//

//   assert(this->initialized());  semiparallel_onlyM();
//   // assert (this->closed());
//   this->close();  int ierr=0;
//   PetscViewer petsc_viewer;
//   ierr = PetscViewerCreate(MPI_COMM_WORLD, &petsc_viewer);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Otherwise the matrix will be dumped to the screen.
//   ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,
//                               PETSC_VIEWER_ASCII_MATLAB);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   ierr = MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Destroy the viewer.
//   ierr = PetscViewerDestroy(petsc_viewer);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
// #ifdef PRINT_INFO
  std::cout << " PetscMatrix::print_hdf5 in file   " <<  name   <<   "\n";
// #endif-
}


// // =====================================================
// void PetscRectangularMatrix::read_hdf5(const std::string namefile,const int mode,const int ml_init) {
//
//   hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
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
// //   H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,pos_row);
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

#endif // #ifdef LIBMESH_HAVE_PETSC
