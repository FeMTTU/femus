#include "FemusExtLib_conf.hpp"

#ifdef FEMUS_HAVE_PETSC

#include "Typedefs_conf.hpp"

// Local includes
#include "petsc_matrixM.hpp"
#include "DenseMatrix.hpp"
#include "petsc_vectorM.hpp"

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "mpi.h"
#pragma GCC diagnostic warning "-Wunused-parameter"



//-----------------------------------------------------------------------
// PetscMatrix members
void PetscMatrixM::init (const unsigned int /*m*/,
                         const unsigned int /*n*/,
                         const unsigned int /*m_l*/,
                         const unsigned int /*n_l*/,
                         const unsigned int /*nnz*/,
                         const unsigned int /*noz*/)
{
//     // We allow 0x0 matrices now
//     //if ((m==0) || (n==0))
//     //  return;
// 
//     // Clear initialized matrices
//     if (this->initialized())  this->clear();
//     this->_is_initialized = true;
// 
// 
//     int ierr     = 0;
//     int m_global = static_cast<int>(m);
//     int n_global = static_cast<int>(n);
//     int m_local  = static_cast<int>(m_l);
//     int n_local  = static_cast<int>(n_l);
//     int n_nz     = static_cast<int>(nnz);
//     int n_oz     = static_cast<int>(noz);

    // create a sequential matrix on one processor
//     if (libMesh::n_processors() == 1) {
//         libmesh_assert ((m_l == m) && (n_l == n));
// 
//         // Create matrix.  Revisit later to do preallocation and make more efficient
//         ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global,
//                                 n_nz, PETSC_NULL, &_mat);
//         CHKERRABORT(MPI_COMM_WORLD,ierr);
//         ierr = MatSetFromOptions (_mat);
//         CHKERRABORT(MPI_COMM_WORLD,ierr);
//     }
// 
//     else{
//         parallel_only();
//         ierr = MatCreateMPIAIJ (MPI_COMM_WORLD, m_local, n_local, m_global, n_global,
//                                 n_nz, PETSC_NULL, n_oz, PETSC_NULL, &_mat);
//         CHKERRABORT(MPI_COMM_WORLD,ierr);
//         ierr = MatSetFromOptions (_mat);
//         CHKERRABORT(MPI_COMM_WORLD,ierr);
//     }
//     this->zero ();
}

// =====================================0
void PetscMatrixM::init () {
    std::cout << "  PetscMatrixM::init: Our PetscMatrix  needs dimensions";
    abort();
}
// =====================================0



void PetscMatrixM::update_sparsity_pattern (const Graph & sparsity_pattern) {
    
    // Clear initialized matrices
    if (this->initialized())    this->clear();
    this->_is_initialized = true;
    int proc_id = 0;    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    
    const unsigned int m   = sparsity_pattern._m;  //this->_dof_map->n_dofs();
    const unsigned int n   = sparsity_pattern._n;
    const unsigned int n_l = sparsity_pattern._nl;
    const unsigned int m_l = sparsity_pattern._ml;
    const unsigned int ml_start = sparsity_pattern._ml_start;
    std::vector<unsigned int> n_nz(m_l);
    std::vector<unsigned int> n_oz(m_l);
   for(uint i=0; i<m_l; i++) {
     uint len = sparsity_pattern[ml_start+i].size()-1;  //this is the real number of nonzero dofs in each row
     n_oz[i]  = sparsity_pattern[ml_start+i][len];      //in the last position we put the number of offset dofs
     n_nz[i]  = len - n_oz[i];
   }
//      for(uint i=0;i<m;i++) n_nz[i]=0;
    
//     for(uint i=0;i<m;i++) n_oz[i]=0;   
//    const unsigned int n_l = this->_dof_map->n_dofs_on_processor(proc_id);
//   const unsigned int m_l = n_l;
//   const std::vector<unsigned int>& n_nz = this->_dof_map->get_n_nz();
//   const std::vector<unsigned int>& n_oz = this->_dof_map->get_n_oz();
//
//   // Make sure the sparsity pattern isn't empty unless the matrix is 0x0
//   libmesh_assert (n_nz.size() == n_l);  libmesh_assert (n_oz.size() == n_l);
//
//   // We allow 0x0 matrices now
//   //if (m==0)
//   //  return;
//
  int ierr     = 0;
  int m_global = static_cast<int>(m);
  int n_global = static_cast<int>(n);
  int m_local  = static_cast<int>(m_l);
  int n_local  = static_cast<int>(n_l);

  int numprocs; MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
  
  if (numprocs == 1)    {
      assert ((m_l == m) && (n_l == n));
      if (n_nz.empty())
        ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global,
			        PETSC_NULL, (int*) PETSC_NULL, &_mat);
      else
        ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global,
			        PETSC_NULL, (int*) &n_nz[0], &_mat);
             CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = MatSetFromOptions (_mat);
             CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  else    {
      parallel_onlyM();
      if (n_nz.empty())
        ierr = MatCreateAIJ (MPI_COMM_WORLD,
			        m_local, n_local,
			        m_global, n_global,
			        PETSC_NULL, (int*) PETSC_NULL,
			        PETSC_NULL, (int*) PETSC_NULL, &_mat);
      else
        ierr = MatCreateAIJ (MPI_COMM_WORLD,
			        m_local, n_local,
			        m_global, n_global,
			        PETSC_NULL, (int*) &n_nz[0],
			        PETSC_NULL, (int*) &n_oz[0], &_mat);
             CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatSetFromOptions (_mat);
             CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    
//   int rank,size,Istart,Iend;
//   ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);//CHKERRQ(ierr);
//   ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);//CHKERRQ(ierr);
//   ierr = MatGetOwnershipRange(_mat,&Istart,&Iend);//CHKERRQ(ierr);
//   printf(" Rostartend %d %d %d \n ",Istart,Iend,rank); 
   this->zero();
   
}

// ===================================
void PetscMatrixM::zero () {
    assert (this->initialized());
    semiparallel_onlyM();
    int ierr=0;
    ierr = MatZeroEntries(_mat);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
}
 
// =======================================================
void PetscMatrixM::zero_rows (std::vector<int> & rows, Real diag_value) {
    assert (this->initialized());
    semiparallel_onlyM();
    int ierr=0;
    if (!rows.empty())
        ierr = MatZeroRows(_mat, rows.size(), &rows[0], diag_value, PETSC_NULL, PETSC_NULL);
    else
        ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value, PETSC_NULL, PETSC_NULL);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// =================================================
void PetscMatrixM::clear () {
    int ierr=0;
    if ((this->initialized()) && (this->_destroy_mat_on_exit)) {
        semiparallel_onlyM();
        ierr = MatDestroy (&_mat);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
        this->_is_initialized = false;
    }
}

// ============================================
Real PetscMatrixM::l1_norm () const {
    assert (this->initialized());

    semiparallel_onlyM();
    int ierr=0;
    PetscReal petsc_value;
    Real value;
    assert (this->closed());
    ierr = MatNorm(_mat, NORM_1, &petsc_value);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    value = static_cast<Real>(petsc_value);
    return value;
}

// ==========================================
Real PetscMatrixM::linfty_norm () const {
    assert (this->initialized());
    semiparallel_onlyM();
    int ierr=0;
    PetscReal petsc_value;
    Real value;
    assert (this->closed());
    ierr = MatNorm(_mat, NORM_INFINITY, &petsc_value);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    value = static_cast<Real>(petsc_value);
    return value;
}



// =====================================================
void PetscMatrixM::print_hdf5(const std::string /*name*/) const {
    assert (this->initialized());
    semiparallel_onlyM();
    // libmesh_assert (this->closed());
    this->close();
    int ierr=0;
    PetscViewer petsc_viewer;
    ierr = PetscViewerCreate (MPI_COMM_WORLD, &petsc_viewer);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

   

    // Otherwise the matrix will be dumped to the screen.
    
        ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                                     PETSC_VIEWER_ASCII_MATLAB);
        CHKERRABORT(MPI_COMM_WORLD,ierr);

        ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
    
    // Destroy the viewer.
    ierr = PetscViewerDestroy (&petsc_viewer);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// ============================================================
void PetscMatrixM::add_matrix(const DenseMatrix& dm,
                              const std::vector<unsigned int>& rows,
                              const std::vector<unsigned int>& cols) {
    assert (this->initialized());
    const unsigned int m = dm.m();
    const unsigned int n = dm.n();
    assert (rows.size() == m);
    assert (cols.size() == n);
    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValues(_mat,
                        m, (int*) &rows[0],
                        n, (int*) &cols[0],
                        (PetscScalar*) &dm.get_values()[0],
                        ADD_VALUES);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
}


// ===========================================================
void PetscMatrixM::_get_submatrix(SparseMatrixM& submatrix,
                                  const std::vector<unsigned int> &rows,
                                  const std::vector<unsigned int> &cols,
                                  const bool reuse_submatrix) const {
    // Can only extract submatrices from closed matrices
    this->close();

    // Make sure the SparseMatrix passed in is really a PetscMatrix
    PetscMatrixM* petsc_submatrix = libmeshM_cast_ptr<PetscMatrixM*>(&submatrix);

    // If we're not reusing submatrix and submatrix is already initialized
    // then we need to clear it, otherwise we get a memory leak.
    if ( !reuse_submatrix && submatrix.initialized() )  submatrix.clear();
    // Construct row and column index sets.
    int ierr=0;
    IS isrow, iscol;
    ierr = ISCreateGeneral(MPI_COMM_WORLD,
                           rows.size(),
                           (int*) &rows[0],
			   PETSC_USE_POINTER,
                           &isrow);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = ISCreateGeneral(MPI_COMM_WORLD,
                           cols.size(),
                           (int*) &cols[0],
			   PETSC_USE_POINTER,
                           &iscol);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

#if !PETSC_VERSION_LESS_THAN(3,0,1) || !PETSC_VERSION_RELEASE
    ierr = MatGetSubMatrix(_mat,
                           isrow,
                           iscol,
//                            PETSC_DECIDE,    //PETSC version 3.1 
                           (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                           &(petsc_submatrix->_mat)); CHKERRABORT(MPI_COMM_WORLD,ierr);
#else
  ierr = MatGetSubMatrix(_mat,
			 isrow,
			 iscol,
			 PETSC_DECIDE,
			 (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
			 &(petsc_submatrix->_mat));  CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif


    // Specify that the new submatrix is initialized and close it.
    petsc_submatrix->_is_initialized = true;
    petsc_submatrix->close();

    // Clean up PETSc data structures
    ierr = ISDestroy(&isrow);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = ISDestroy(&iscol);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

 
    

}


// ======================================================
void PetscMatrixM::get_diagonal (NumericVectorM& dest) const {
    // Make sure the NumericVector passed in is really a PetscVector
    PetscVectorM& petsc_dest = libmeshM_cast_ref<PetscVectorM&>(dest);
    // Call PETSc function.
// #if PETSC_VERSION_LESS_THAN(2,3,1)
// 
//     std::cout << "This method has been developed with PETSc 2.3.1.  "
//               << "No one has made it backwards compatible with older "
//               << "versions of PETSc so far; however, it might work "
//               << "without any change with some older version." << std::endl;
//     libmesh_error();
// #else
    // Needs a const_cast since PETSc does not work with const.
    int ierr =  MatGetDiagonal(const_cast<PetscMatrixM*>(this)->mat(),petsc_dest.vec());
    CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
}



// ===========================================================
void PetscMatrixM::get_transpose (SparseMatrixM& dest) const {
    // Make sure the SparseMatrix passed in is really a PetscMatrix
    PetscMatrixM& petsc_dest = libmeshM_cast_ref<PetscMatrixM&>(dest);

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
