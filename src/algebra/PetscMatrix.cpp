/*=========================================================================

 Program: FEMUS
 Module: PetscMatrix
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

// Local includes
#include "Parallel.hpp"
#include "PetscMatrix.hpp"
#include "PetscVector.hpp"
#include "DenseMatrix.hpp"
#include <mpi.h>
#include <hdf5.h>
#include <sstream>


namespace femus {


  using namespace std;

//-----------------------------------------------------------------------
// PetscMatrix members
  void PetscMatrix::init (const  int m,
                          const  int n,
                          const  int m_l,
                          const  int n_l,
                          const  int nnz,
                          const  int noz) {
    // Set matrix dimension
    _m = m;
    _n = n;
    _m_l = m_l;
    _n_l = n_l;

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
    MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);
    int numprocs;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    int ierr     = 0;
    int m_global = static_cast<int> (m);
    int n_global = static_cast<int> (n);
    int m_local  = static_cast<int> (m_l);
    int n_local  = static_cast<int> (n_l);
    int n_nz     = static_cast<int> (nnz);
    int n_oz     = static_cast<int> (noz);

// create a sequential matrix on one processor
    if (numprocs == 1) {
      assert ( (m_l == m) && (n_l == n));

      // Create matrix.  Revisit later to do preallocation and make more efficient
      ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global,
                              n_nz, PETSC_NULL, &_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
      ierr = MatSetFromOptions (_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }

    else {
      parallel_only();

      ierr = MatCreate (MPI_COMM_WORLD, &_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);

      ierr = MatSetSizes (_mat, m_l, n_l, m, n);
      CHKERRABORT (MPI_COMM_WORLD, ierr);

      ierr = MatSetType (_mat, MATMPIAIJ); // Automatically chooses seqaij or mpiaij
      CHKERRABORT (MPI_COMM_WORLD, ierr);

      ierr = MatMPIAIJSetPreallocation (_mat, nnz, PETSC_NULL, noz, PETSC_NULL);
      CHKERRABORT (MPI_COMM_WORLD, ierr);

    }

    this->zero();
  }
  
  
  void PetscMatrix::init (const int nr, const int nc, const std::vector < SparseMatrix*> &P){
    Mat KK[3];
    unsigned dim = P.size();
    for (unsigned k = 0; k < dim; k++) {
      KK[k] = (static_cast<PetscMatrix*> (P[k]))->mat();
    }
    MatCreateNest (MPI_COMM_WORLD, dim, NULL, 1, NULL, KK, &_mat);
    this->_is_initialized = true;
    MatGetSize(_mat, &_n, &_m);
    MatGetLocalSize(_mat, &_n_l, &_m_l);
    _destroy_mat_on_exit = true;
    
    std::cout << _n <<" " <<  _m << " " << _n_l << " " << _m_l << std::endl;
  }
  

  void PetscMatrix::init (const  int m, const  int n, const  int m_l, const  int n_l,
                          const std::vector< int > & n_nz, const std::vector< int > & n_oz) {
    // Set matrix dimension
    _m = m;
    _n = n;
    _m_l = m_l;
    _n_l = n_l;

    // Clear initialized matrices
    if (this->initialized())
      this->clear();

    this->_is_initialized = true;

    // processor info
    int n_procs;
    MPI_Comm_size (MPI_COMM_WORLD, &n_procs);

    int ierr = 0;

// create a sequential matrix on one processor
    if (n_procs == 1) {
      assert (n_nz.size() == _m_l);
      ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, _m, _n, 0, &n_nz[0], &_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
      ierr = MatSetFromOptions (_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
    else {
      parallel_only();
      assert ( (n_nz.size() == _m_l) && (n_oz.size() == _m_l));
      ierr = MatCreate (MPI_COMM_WORLD, &_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
      ierr = MatSetSizes (_mat, _m_l, _n_l, _m, _n);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
      ierr = MatSetType (_mat, MATMPIAIJ);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
      ierr = MatMPIAIJSetPreallocation (_mat, 1, &n_nz[0], 100, &n_oz[0]);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
    this->zero();
  }

// =====================================0
  void PetscMatrix::update_sparsity_pattern (
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
    MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);
    int numprocs;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    int ierr     = 0;

    // create a sequential matrix on one processor -------------------
    if (numprocs == 1)    {
      assert ( (m_local == m_global) && (n_local == n_global));
      if (n_nz.empty())
        ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global,
                                PETSC_DEFAULT, (int*) PETSC_NULL, &_mat);
      else
        ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global,
                                PETSC_DEFAULT, (int*) &n_nz[0], &_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);

      ierr = MatSetFromOptions (_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
    else    {   // multi processors ----------------------------
      parallel_only();             //TODO
      if (n_nz.empty()) {

//old piece of code with Petsc lib version 3.1
//       ierr = MatCreateMPIAIJ(MPI_COMM_WORLD,
//                              m_local, n_local,
//                              m_global, n_global,
//                              PETSC_NULL, (int*) PETSC_NULL,
//                              PETSC_NULL, (int*) PETSC_NULL,
//                              &_mat
//                             );

        ierr = MatCreate (MPI_COMM_WORLD, &_mat);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = MatSetSizes (_mat, m_local, n_local, m_global, n_global);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = MatSetType (_mat, MATMPIAIJ);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = MatMPIAIJSetPreallocation (_mat, 0, 0, 0, 0);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

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

        ierr = MatCreate (MPI_COMM_WORLD, &_mat);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = MatSetSizes (_mat, m_local, n_local, m_global, n_global);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = MatSetType (_mat, MATMPIAIJ); // Automatically chooses seqaij or mpiaij
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = MatMPIAIJSetPreallocation (_mat, 0, (int*) &n_nz[0], 0, (int*) &n_oz[0]);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

      }

      ierr = MatSetFromOptions (_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    } // ----------------------------------------------------

    this->zero();
    return;
  }

// =====================================0
  void PetscMatrix::update_sparsity_pattern (const Graph &sparsity_pattern) {

    // dimension
    const  int m   =  sparsity_pattern._m; // global rows;
    const  int n   =  sparsity_pattern._n; // global columns
    const  int n_l =  sparsity_pattern._nl;  // local rows
    const  int m_l =  sparsity_pattern._ml;  // local columns
    const  int ml_start = sparsity_pattern._ml_start;

    // vectors n_nz (diagonal) and n_oz (offset) -------------------------------
    std::vector< int> n_nz (m_l);
    std::vector< int> n_oz (m_l);
    for (int i = 0; i < m_l; i++) {
      int len = sparsity_pattern[ml_start + i].size() - 1;
      n_oz[i] = sparsity_pattern[ml_start + i][len];
      n_nz[i] = len - n_oz[i];
    }
    update_sparsity_pattern (m, n, m_l, n_l, n_oz, n_nz);
    return;
  }

  void PetscMatrix::update_sparsity_pattern_old (const Graph & sparsity_pattern) {

    // Clear initialized matrices
    if (this->initialized())    this->clear();
    this->_is_initialized = true;
    int proc_id = 0;
    MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);

    const unsigned int m   = sparsity_pattern._m;  //this->_dof_map->n_dofs();
    const unsigned int n   = sparsity_pattern._n;
    const unsigned int n_l = sparsity_pattern._nl;
    const unsigned int m_l = sparsity_pattern._ml;
    const unsigned int ml_start = sparsity_pattern._ml_start;

    //     std::vector<unsigned int> n_nz(m_l); //TODO eugenio
    //     std::vector<unsigned int> n_oz(m_l); //TODO eugenio

    std::vector<int> n_nz (m_l); //TODO eugenio
    std::vector<int> n_oz (m_l); //TODO eugenio

    for (uint i = 0; i < m_l; i++) {
      uint len = sparsity_pattern[ml_start + i].size() - 1; //this is the real number of nonzero dofs in each row
      n_oz[i]  = sparsity_pattern[ml_start + i][len];    //in the last position we put the number of offset dofs
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
    int m_global = static_cast<int> (m);
    int n_global = static_cast<int> (n);
    int m_local  = static_cast<int> (m_l);
    int n_local  = static_cast<int> (n_l);

    int numprocs;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    if (numprocs == 1)    {
      assert ( (m_l == m) && (n_l == n));
      if (n_nz.empty())
//         ierr = MatCreateSeqAIJ(MPI_COMM_WORLD, m_global, n_global, //TODO eugenio
//                                PETSC_NULL, (int*) PETSC_NULL, &_mat);
        ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global, PETSC_DEFAULT, PETSC_NULL, &_mat); //TODO eugenio
      else
//         ierr = MatCreateSeqAIJ(MPI_COMM_WORLD, m_global, n_global, //TODO eugenio
//                                PETSC_NULL, (int*) &n_nz[0], &_mat);
        ierr = MatCreateSeqAIJ (MPI_COMM_WORLD, m_global, n_global, PETSC_DEFAULT, &n_nz[0], &_mat); //TODO eugenio
      CHKERRABORT (MPI_COMM_WORLD, ierr);

      ierr = MatSetFromOptions (_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
    else    {
      parallel_only();
      if (n_nz.empty())
//         ierr = MatCreateAIJ(MPI_COMM_WORLD, //TODO eugenio
//                             m_local, n_local,
//                             m_global, n_global,
//                             PETSC_NULL, (int*) PETSC_NULL,
//                             PETSC_NULL, (int*) PETSC_NULL, &_mat);

        ierr = MatCreateAIJ (MPI_COMM_WORLD, //TODO eugenio
                             m_local, n_local,
                             m_global, n_global,
                             PETSC_DEFAULT, PETSC_NULL,
                             PETSC_DEFAULT, PETSC_NULL, &_mat);
      else
//         ierr = MatCreateAIJ(MPI_COMM_WORLD, //TODO eugenio
//                             m_local, n_local,
//                             m_global, n_global,
//                             PETSC_NULL, (int*) &n_nz[0],
//                             PETSC_NULL, (int*) &n_oz[0], &_mat);

        ierr = MatCreateAIJ (MPI_COMM_WORLD, //TODO eugenio
                             m_local, n_local,
                             m_global, n_global,
                             PETSC_DEFAULT, &n_nz[0],
                             PETSC_DEFAULT, &n_oz[0], &_mat);

      CHKERRABORT (MPI_COMM_WORLD, ierr);
      ierr = MatSetFromOptions (_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }

//   int rank,size,Istart,Iend;
//   ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);//CHKERRQ(ierr);
//   ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);//CHKERRQ(ierr);
//   ierr = MatGetOwnershipRange(_mat,&Istart,&Iend);//CHKERRQ(ierr);
//   printf(" Rostartend %d %d %d \n ",Istart,Iend,rank);
    this->zero();

  }

// ==========================================
/// This function sets all the entries to 0
  void PetscMatrix::zero() {
    assert (this->initialized());
    semiparallel_only();
    int ierr = 0;
    ierr = MatZeroEntries (_mat);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

// =========================================================================
/// This function sets all row entries to 0 then puts diag_value in the diagonal entry
  void PetscMatrix::zero_rows (std::vector<int> & rows, double diag_value) {
    assert (this->initialized());
    semiparallel_only();
    int ierr = 0;
    if (!rows.empty()) ierr = MatZeroRows (_mat, rows.size(), &rows[0], diag_value, 0, 0); // add,0,0,)    !!!!
    else   ierr = MatZeroRows (_mat, 0, PETSC_NULL, diag_value, 0, 0);                    // add,0,0,)    !!!!
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

// =================================================
  void PetscMatrix::clear() {
    int ierr = 0;
    if ( (this->initialized()) && (this->_destroy_mat_on_exit)) {
      semiparallel_only();
      ierr = MatDestroy (&_mat);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
      this->_is_initialized = false;
    }
  }

// ============================================
  double PetscMatrix::l1_norm() const {
    assert (this->initialized());
    semiparallel_only();
    int ierr = 0;
    PetscReal petsc_value;
    double value;
    assert (this->closed());
    ierr = MatNorm (_mat, NORM_1, &petsc_value);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    value = static_cast<double> (petsc_value);
    return value;
  }

// ==========================================
  double PetscMatrix::linfty_norm() const {
    assert (this->initialized());
    semiparallel_only();
    int ierr = 0;
    PetscReal petsc_value;
    double value;
    assert (this->closed());
    ierr = MatNorm (_mat, NORM_INFINITY, &petsc_value);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    value = static_cast<double> (petsc_value);
    return value;
  }

// ===================================================
  void PetscMatrix::print_personal (
    std::ostream& os // pointer stream
  ) const {// =========================================
    assert (this->initialized());
    std::cout << "\n PetscMatrix::print_personal \n";
// #ifndef NDEBUG
//     if(os != std::cout)
//       std::cerr << "Warning! PETSc can only print to std::cout!" << std::endl;
// #endif
    int ierr = 0;
    ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }




  void PetscMatrix::print_matlab (const std::string& name, const std::string& format) const {

    this->close();

    PetscErrorCode ierr = 0;
    PetscViewer petsc_viewer;
    ierr = PetscViewerCreate (MPI_COMM_WORLD,
                              &petsc_viewer);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

    /**
     * Create a binary file containing the matrix
     * if a filename was provided.
     */
    if (name != "") {

      if (format == "binary") {

        ierr = PetscViewerBinaryOpen (MPI_COMM_WORLD,
                                      name.c_str(),
                                      FILE_MODE_WRITE,
                                      &petsc_viewer);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = MatView (_mat, petsc_viewer);
        CHKERRABORT (MPI_COMM_WORLD, ierr);
      }
      else if (format == "ascii") {
        ierr = PetscViewerASCIIOpen (MPI_COMM_WORLD,
                                     name.c_str(),
                                     &petsc_viewer);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = PetscViewerPushFormat (petsc_viewer,
                                      PETSC_VIEWER_ASCII_MATLAB);
        CHKERRABORT (MPI_COMM_WORLD, ierr);

        ierr = MatView (_mat, petsc_viewer);
        CHKERRABORT (MPI_COMM_WORLD, ierr);
      }
      else {
        std::cout << "Provide either \"ascii\" or \"binary\" for the second argument" << std::endl;
        abort();
      }

    }

    /**
     * Otherwise the matrix will be dumped to the screen, regardless of the format you provide
     */
    else {
      ierr = PetscViewerPushFormat (PETSC_VIEWER_STDOUT_WORLD,
                                    PETSC_VIEWER_ASCII_MATLAB);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
      ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }

    /**
     * Destroy the viewer.
     */
    ierr = PetscViewerDestroy (&petsc_viewer);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

  }




// =====================================================
  void PetscMatrix::print_hdf5 (const std::string /*name*/) const {




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
  void PetscMatrix::add_matrix (const DenseMatrix& dm,
                                const std::vector<unsigned int>& rows,
                                const std::vector<unsigned int>& cols) {
    assert (this->initialized());
    const  int m = dm.m();
    assert ( (int) rows.size() == m);
    const  int n = dm.n();
    assert ( (int) cols.size() == n);

    int ierr = 0;
    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValues (_mat,
                         m, (int*) &rows[0],
                         n, (int*) &cols[0],
                         (PetscScalar*) &dm.get_values() [0],
                         ADD_VALUES);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

// ============================================================
/// This Petsc adds a dense matrix to a sparse matrix
  void PetscMatrix::add_matrix_blocked (
    const std::vector< double > &mat_values,  // blocked matrix stored as an m*n array
    const std::vector< int>& rows, // row vector indexes
    const std::vector< int>& cols) { // column vector indices
    // ==================================================
    assert (this->initialized());

    const int m = rows.size();
    const int n = cols.size();
    assert (m * n == mat_values.size());

    MatSetValuesBlocked (_mat, m, &rows[0], n, &cols[0], &mat_values[0], ADD_VALUES);

    return;
  }

  void PetscMatrix::add_matrix_blocked (
    const std::vector< double > &mat_values,  // blocked matrix stored as an m*n array
    const std::vector< unsigned >& rows, // row vector indexes
    const std::vector< unsigned >& cols) { // column vector indices
    // ==================================================
    assert (this->initialized());

    const int m = rows.size();
    const int n = cols.size();
    assert (m * n == mat_values.size());

    MatSetValuesBlocked (_mat, m, (int*) &rows[0], n, (int*) &cols[0], &mat_values[0], ADD_VALUES);

    return;
  }

// // ============================================================

  void PetscMatrix::matrix_PtAP (const SparseMatrix &mat_P, const SparseMatrix &mat_A, const bool &mat_reuse) {

    const PetscMatrix* A = static_cast<const PetscMatrix*> (&mat_A);
    A->close();

    const PetscMatrix* P = static_cast<const PetscMatrix*> (&mat_P);
    P->close();

    int ierr = 0;
    if (mat_reuse) {
      ierr = MatPtAP (const_cast<PetscMatrix*> (A)->mat(), const_cast<PetscMatrix*> (P)->mat(), MAT_REUSE_MATRIX, 1.0, &_mat);
    }
    else {
      this->clear();
      ierr = MatPtAP (const_cast<PetscMatrix*> (A)->mat(), const_cast<PetscMatrix*> (P)->mat(), MAT_INITIAL_MATRIX , 1.0, &_mat);
      this->_is_initialized = true;
    }
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

// // ============================================================

  void PetscMatrix::matrix_ABC (const SparseMatrix &mat_A, const SparseMatrix &mat_B, const SparseMatrix &mat_C, const bool &mat_reuse) {

    const PetscMatrix* A = static_cast<const PetscMatrix*> (&mat_A);
    A->close();

    const PetscMatrix* B = static_cast<const PetscMatrix*> (&mat_B);
    B->close();

    const PetscMatrix* C = static_cast<const PetscMatrix*> (&mat_C);
    C->close();

    int ierr = 0;
    if (mat_reuse) {
      ierr = MatMatMatMult (const_cast<PetscMatrix*> (A)->mat(), const_cast<PetscMatrix*> (B)->mat(),
                            const_cast<PetscMatrix*> (C)->mat(), MAT_REUSE_MATRIX, 1.0, &_mat);
    }
    else {
      this->clear();
      ierr = MatMatMatMult (const_cast<PetscMatrix*> (A)->mat(), const_cast<PetscMatrix*> (B)->mat(),
                            const_cast<PetscMatrix*> (C)->mat(), MAT_INITIAL_MATRIX, 1.0, &_mat);
      this->_is_initialized = true;
    }
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

  void PetscMatrix::matrix_RightMatMult (const SparseMatrix &mat_A) {

    PetscMatrix* A = const_cast< PetscMatrix* > (static_cast < const PetscMatrix* > (&mat_A));
    A->close();
    Mat matCopy;

    int ierr = MatConvert (_mat, MATSAME, MAT_INITIAL_MATRIX, &matCopy);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

    this->clear();

    MatMatMult (matCopy, A->mat(), MAT_INITIAL_MATRIX, 1.0, &_mat);

    ierr = MatDestroy (&matCopy);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

    this->_is_initialized = true;

  }

  void PetscMatrix::matrix_LeftMatMult (const SparseMatrix &mat_A) {

    PetscMatrix* A = const_cast< PetscMatrix* > (static_cast < const PetscMatrix* > (&mat_A));
    A->close();
    Mat matCopy;

    int ierr = MatConvert (_mat, MATSAME, MAT_INITIAL_MATRIX, &matCopy);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

    this->clear();

    MatMatMult (A->mat(), matCopy, MAT_INITIAL_MATRIX, 1.0, &_mat);

    ierr = MatDestroy (&matCopy);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

    this->_is_initialized = true;

  }


// ===========================================================

  void PetscMatrix::matrix_get_diagonal_values (const std::vector< int > &index, std::vector<double> &value) const {
    assert (index.size() == value.size());
    for (int i = 0; i < index.size(); i++) {
      int ierr = MatGetValues (_mat, 1, &index[i], 1, &index[i], &value[i]);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
  }

  // ===========================================================
  
  void PetscMatrix::matrix_set_diagonal_values (NumericVector& D) {
    MatDiagonalSet(_mat, (static_cast< PetscVector& > (D)).vec(),INSERT_VALUES);
  }
  
  
// ===========================================================

  void PetscMatrix::matrix_set_diagonal_values (const std::vector< int > &index, const double &value) {
    for (int i = 0; i < index.size(); i++) {
      int ierr = MatSetValuesBlocked (_mat, 1, &index[i], 1, &index[i], &value, INSERT_VALUES);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
  }

// ===========================================================

  void PetscMatrix::matrix_set_diagonal_values (const std::vector< int > &index, const std::vector<double> &value) {
    assert (index.size() == value.size());
    for (int i = 0; i < index.size(); i++) {
      int ierr = MatSetValuesBlocked (_mat, 1, &index[i], 1, &index[i], &value[i], INSERT_VALUES);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
  }

// ===========================================================
  void PetscMatrix::matrix_set_off_diagonal_values_blocked (const std::vector< int > &index_rows, const std::vector< int > &index_cols, const double &value) {
    assert (index_rows.size() == index_cols.size());
    for (int i = 0; i < index_rows.size(); i++) {
      int ierr = MatSetValuesBlocked (_mat, 1, &index_rows[i], 1, &index_cols[i], &value, INSERT_VALUES);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
  }

// ===========================================================
  void PetscMatrix::matrix_set_off_diagonal_values_blocked (const std::vector< int > &index_rows, const std::vector< int > &index_cols, const std::vector<double> &value) {
    assert (index_rows.size() == index_cols.size());
    assert (index_rows.size() == value.size());
    for (int i = 0; i < index_rows.size(); i++) {
      int ierr = MatSetValuesBlocked (_mat, 1, &index_rows[i], 1, &index_cols[i], &value[i], INSERT_VALUES);
      CHKERRABORT (MPI_COMM_WORLD, ierr);
    }
  }


// ===========================================================
/// This function either creates or re-initializes a matrix called "submatrix".
  void PetscMatrix::_get_submatrix (SparseMatrix& submatrix,
                                    const std::vector< int> &rows,
                                    const std::vector< int> &cols,
                                    const bool reuse_submatrix) const {

    // Can only extract submatrices from closed matrices
    this->close();

    // Make sure the SparseMatrix passed in is really a PetscMatrix
    PetscMatrix* petsc_submatrix = libmeshM_cast_ptr<PetscMatrix*> (&submatrix);

    // If we're not reusing submatrix and submatrix is already initialized
    // then we need to clear it, otherwise we get a memory leak.
    if (!reuse_submatrix && submatrix.initialized())  submatrix.clear();
    // Construct row and column index sets.
    int ierr = 0;
    IS isrow, iscol;

    ierr = ISCreateGeneral (MPI_COMM_WORLD, rows.size(), (int*) &rows[0], PETSC_COPY_VALUES, &isrow); // PETSC_COPY_VALUES is my first choice; see also PETSC_OWN_POINTER, PETSC_USE_POINTER
    CHKERRABORT (MPI_COMM_WORLD, ierr);

    ierr = ISCreateGeneral (MPI_COMM_WORLD, cols.size(), (int*) &cols[0], PETSC_COPY_VALUES, &iscol);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

//---

#if !PETSC_VERSION_LESS_THAN(3,0,1) || !PETSC_VERSION_RELEASE
//     ierr = MatGetSubMatrix(_mat,
//                            isrow,
//                            iscol,
//                            (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
//                            &(petsc_submatrix->_mat));
    ierr = MatCreateSubMatrix (_mat,
                               isrow,
                               iscol,
                               (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                               & (petsc_submatrix->_mat));
    CHKERRABORT (MPI_COMM_WORLD, ierr);
#else
    ierr = MatGetSubMatrix (_mat,
                            isrow,
                            iscol,
                            PETSC_DECIDE,
                            (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                            & (petsc_submatrix->_mat));
    CHKERRABORT (MPI_COMM_WORLD, ierr);
#endif

    // Specify that the new submatrix is initialized and close it.
    petsc_submatrix->_is_initialized = true;
    petsc_submatrix->close();

    // Clean up PETSc data structures
    ierr = ISDestroy (&isrow);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    ierr = ISDestroy (&iscol);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

    return;
  }

// ======================================================
  void PetscMatrix::get_diagonal (NumericVector& dest) const {
    // Make sure the NumericVector passed in is really a PetscVector
    PetscVector& petsc_dest = libmeshM_cast_ref<PetscVector&> (dest);
    // Call PETSc function.
    // Needs a const_cast since PETSc does not work with const.
    int ierr =  MatGetDiagonal (const_cast<PetscMatrix*> (this)->mat(), petsc_dest.vec());
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    return;
  }

/// This function copies the transpose of the matrix
// ===========================================================
  void PetscMatrix::get_transpose (
    SparseMatrix& dest  // output transpose matrix
  ) const {

    // Make sure the SparseMatrix passed in is really a PetscMatrix
    PetscMatrix& petsc_dest = libmeshM_cast_ref<PetscMatrix&> (dest);

    // If we aren't reusing the matrix then need to clear dest,
    // otherwise we get a memory leak
    if (&petsc_dest != this)    dest.clear();

    int ierr;
#if PETSC_VERSION_LESS_THAN(3,0,0)
    if (&petsc_dest == this)
      ierr = MatTranspose (_mat, PETSC_NULL);
    else
      ierr = MatTranspose (_mat, &petsc_dest._mat);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
#else
    // FIXME - we can probably use MAT_REUSE_MATRIX in more situations
    if (&petsc_dest == this)
      //ierr = MatTranspose(_mat, MAT_REUSE_MATRIX, &petsc_dest._mat);
      ierr = MatTranspose (_mat, MAT_INPLACE_MATRIX, &petsc_dest._mat);
    else
      ierr = MatTranspose (_mat, MAT_INITIAL_MATRIX, &petsc_dest._mat);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
#endif

    int temp = _m;
    petsc_dest._m = _n;
    petsc_dest._n = temp;

    temp = _m_l;
    petsc_dest._m_l = _n_l;
    petsc_dest._n_l = temp;

    // Specify that the transposed matrix is initialized and close it.
    petsc_dest._is_initialized = true;
    petsc_dest.close();
  }


  void PetscMatrix::mat_zero_rows (const std::vector <int> &index, const double &diagonal_value) const {
    MatSetOption (_mat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);
    MatSetOption (_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    MatZeroRows (_mat, index.size(), &index[0], diagonal_value, 0, 0);
  }



} //end namespace femus



#endif // #ifdef LIBMESH_HAVE_PETSC
