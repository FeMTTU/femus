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

#ifndef __femus_algebra_PetscMatrix_hpp__
#define __femus_algebra_PetscMatrix_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

// ======================================
#ifdef HAVE_PETSC
// ======================================

// This class
#include "SparseMatrix.hpp"
#include "PetscMacro.hpp"
// Petsc include files.
EXTERN_C_FOR_PETSC_BEGIN
# include <petscmat.h>
EXTERN_C_FOR_PETSC_END

// C++ includes -------------------
#include <algorithm> // equal_range

// Local includes
#include "Parallel.hpp"        // parallel commands

// Macro to identify and debug functions which should be called in
// parallel on parallel matrices but which may be called in serial on
// serial matrices.  This macro will only be valid inside non-static
// PetscMatrix methods
#undef semiparallel_only
#ifndef NDEBUG
#include <cstring>

#define semiparallel_only() do { if (this->initialized()) { const char *mytype; \
      MatGetType(_mat,&mytype); \
      if (!strcmp(mytype, MATSEQAIJ)) \
        parallel_only(); } } while (0)
#else
#define semiparallel_only()
#endif



namespace femus {





// =======================================================
// Petsc matrix
// =======================================================

  class PetscMatrix : public SparseMatrix {

    private:
      // data ------------------------------------
      Mat _mat;                 ///< Petsc matrix pointer
      bool _destroy_mat_on_exit;///< Boolean value (false)

    public:
      // Constructor ---------------------------------------------------------
      /// Constructor I;  initialize the matrix before usage with \p init(...).
      PetscMatrix();
      /// Constructor II.
      PetscMatrix (Mat m);

      /// Initialize a Petsc matrix
      void init (const int m, const int n, const int m_l, const int n_l,
                 const int nnz = 0, const int noz = 0);
      void init (const  int m, const  int n, const  int m_l, const  int n_l,
                 const std::vector< int > & n_nz, const std::vector< int > & n_oz);

      void init (const int nr, const int nc, const std::vector < SparseMatrix*> &P);
            
      void init (const int m,  const int n) {
        _m = m;
        _n = n;
      }
      void init () {};
      // Destructors ----------------------------
      /// Destructor
      ~PetscMatrix();

      void clear(); /// Release all memory
      void zero();///< set to zero
      void zero_rows (std::vector<int> & rows, double diag_value = 0.0); ///< set  rows to zero
      void close() const;///< close


      // Returns -------------------------------------------
      /// PETSc matrix context pointer
      Mat mat() {
        assert (_mat != NULL);
        return _mat;
      }

      // matrix dimensions
      int m() const;          ///< row-dimension
      int n() const;          ///< column dimension
      int row_start() const;  ///< row-start
      int row_stop() const;   ///< row-stop

      ///    Return the value of the entry  (i,j). Do not use
      double operator() (const int i, const int j) const;
      /// Return the row values
      int MatGetRowM (const int i_val, int cols[] = PETSC_NULL , double vals[] = PETSC_NULL);

      // Setting -------------------------------------
      /** @deprecated */
      void update_sparsity_pattern_old (const Graph & sparsity_pattern);

      void update_sparsity_pattern (const Graph &sparsity_pattern); ///<   sparsity patter update (Graph)

      void update_sparsity_pattern (int m, int n, int m_l, int n_l, ///<   sparsity patter update (petsc)
                                    const std::vector<int>  n_oz, const std::vector<int>  n_nz);
      // set values
      void set (const int i, const int j, const double value); ///< Set the value.
      void add (const int i, const int j, const double value); ///<  add  value
      // add
      /// Add the full matrix to the Petsc matrix.
      void add_matrix (const DenseMatrix &dm,
                       const std::vector<unsigned int> &rows,
                       const std::vector<unsigned int> &cols);
      /// Add the full matrix to the Petsc matrix.
      void add_matrix (const DenseMatrix &dm,
                       const std::vector<unsigned int> &dof_indices);

      /// Add a Sparse matrix
      void add (const double a, SparseMatrix &X);
      // Add a row to a Sparse matrix
      void insert_row (const int row, const int ncols,
                       const std::vector<int>& cols, double* values);

      void add_matrix_blocked (const std::vector< double > &mat_value,
                               const std::vector< int> &rows,
                               const std::vector< int> &cols);

      void add_matrix_blocked (const std::vector< double > &mat_values,  // blocked matrix stored as an m*n array
                               const std::vector< unsigned >& rows, // row vector indexes
                               const std::vector< unsigned >& cols);

      void matrix_add (const double a_in, SparseMatrix &X_in, const char pattern []);

      void matrix_PtAP (const SparseMatrix &mat_P, const SparseMatrix &mat_A, const bool &reuse);
      void matrix_ABC (const SparseMatrix &mat_A, const SparseMatrix &mat_B, const SparseMatrix &mat_C, const bool &reuse);

      void matrix_RightMatMult (const SparseMatrix &mat_A);
      void matrix_LeftMatMult (const SparseMatrix &mat_A);

      void matrix_get_diagonal_values (const std::vector< int > &index, std::vector<double> &value) const ;
      void matrix_set_diagonal_values (NumericVector& D);
      void matrix_set_diagonal_values (const std::vector< int > &index, const double &value);
      void matrix_set_diagonal_values (const std::vector< int > &index, const std::vector<double> &value);
      void matrix_set_off_diagonal_values_blocked (const std::vector< int > &index_rows, const std::vector< int > &index_cols, const double &value);
      void matrix_set_off_diagonal_values_blocked (const std::vector< int > &index_rows, const std::vector< int > &index_cols, const std::vector<double> &value);

      // functions ------------------------------------
      /// Return the l1-norm of the matrix
      double l1_norm() const;
      /// Return the linfty-norm of the matrix
      double linfty_norm() const;

      /// Petsc matrix has been closed and fully assembled
      bool closed() const;


      // print -----------------------------------------------------
      void print_personal (std::ostream& os = std::cout) const; ///< print personal
      void print_personal (const std::string name = "NULL") const; ///< print
      void print_hdf5 (const std::string name = "NULL") const;  ///< print hdf5
      /** Print  to Matlab format */
      void print_matlab (const std::string& name, const std::string& format) const;

      // functions ---------------------------------------------------
      /// Copies the diagonal part of the matrix
      virtual void get_diagonal (NumericVector& dest) const;
      /// Transpose Matrix.
      virtual void get_transpose (SparseMatrix& dest) const;
      /// Swaps the raw PETSc matrix context pointers.
      void mat_zero_rows (const std::vector <int> &index, const double &diagonal_value) const;

      void swap (PetscMatrix &);


// ========================================================


    protected:
      ///  get "submatrix" from sparse matrix.
      virtual void _get_submatrix (SparseMatrix& submatrix,
                                   const std::vector<int>& rows,
                                   const std::vector<int>& cols,
                                   const bool reuse_submatrix) const;


  };

// ===============================================
// PetscMatrix inline members
// ===============================================

// ===============================================
  inline PetscMatrix::PetscMatrix()  : _destroy_mat_on_exit (true) {}

// =================================================================
  inline PetscMatrix::PetscMatrix (Mat m) : _destroy_mat_on_exit (false) {
    this->_mat = m;
    this->_is_initialized = true;
  }

// ===============================================
  inline PetscMatrix::~PetscMatrix() {
    this->clear();
  }

// ==========================================
/// This function checks if the matrix is closed
  inline void PetscMatrix::close() const {
    parallel_only();
    int ierr = 0;
    ierr = MatAssemblyBegin (_mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    ierr = MatAssemblyEnd (_mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

// ==================================================
/// This function returns the row number
  inline int PetscMatrix::m() const {
    assert (this->initialized());
    int petsc_m = 0, petsc_n = 0, ierr = 0;
    ierr = MatGetSize (_mat, &petsc_m, &petsc_n);
    return static_cast<int> (petsc_m);
  }

// ============================================
/// This function returns the column number
  inline int PetscMatrix::n() const {
    assert (this->initialized());
    int petsc_m = 0, petsc_n = 0, ierr = 0;
    ierr = MatGetSize (_mat, &petsc_m, &petsc_n);
    return static_cast<int> (petsc_n);
  }

// ==========================================
/// This function returns the start row location
  inline int PetscMatrix::row_start() const {
    assert (this->initialized());
    int start = 0, stop = 0, ierr = 0;
    ierr = MatGetOwnershipRange (_mat, &start, &stop);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    return static_cast<int> (start);
  }

// =========================================
/// This function returns the stop row location
  inline int PetscMatrix::row_stop() const {
    assert (this->initialized());
    int start = 0, stop = 0, ierr = 0;
    ierr = MatGetOwnershipRange (_mat, &start, &stop);
    CHKERRABORT (MPI_COMM_WORLD, ierr);

    return static_cast<int> (stop);
  }

// ======================================================
/// This function sets the value in the (i,j) pos
  inline void PetscMatrix::set (const int i,
                                const int j,
                                const double value) {
    assert (this->initialized());
    int ierr = 0, i_val = i, j_val = j;
    PetscScalar petsc_value = static_cast<PetscScalar> (value);
    ierr = MatSetValues (_mat, 1, &i_val, 1, &j_val,
                         &petsc_value, INSERT_VALUES);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

// =================================================
/// This function add the  value to the element (i,j).
  inline void PetscMatrix::add (const int i,   // index i
                                const int j, // index j
                                const double value      // value
                               ) {
    assert (this->initialized());
    int ierr = 0, i_val = i, j_val = j;

    PetscScalar petsc_value = static_cast<PetscScalar> (value);
    ierr = MatSetValues (_mat, 1, &i_val, 1, &j_val,
                         &petsc_value, ADD_VALUES);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

// =====================================================
/// This function adds a square dense matrix to the main sparse matrix
  inline void PetscMatrix::add_matrix (
    const DenseMatrix& dm,                       // dense matrix
    const std::vector<unsigned int>& dof_indices  // index vector
  ) {
    this->add_matrix (dm, dof_indices, dof_indices);
  }

// =========================================================
/// This function performs \f$\texttt{this} = a*X + \texttt{this} \f$.
  inline void PetscMatrix::add (const double a_in,     // double constant
                                SparseMatrix &X_in  // sparse matrix
                               ) {
    assert (this->initialized());

    // crash due to incompatible sparsity structure...
    assert (this->m() == X_in.m());
    assert (this->n() == X_in.n());

    PetscScalar     a = static_cast<PetscScalar> (a_in);
    PetscMatrix*   X = dynamic_cast<PetscMatrix*> (&X_in);

    // the matrix from which we copy the values has to be assembled/closed
    // X->close ();
    assert (X != NULL);
    assert (X->closed());
    semiparallel_only();
    int ierr = 0;

#if PETSC_VERSION_LESS_THAN(2,3,0)  // 2.2.x & earlier style  
    ierr = MatAXPY (&a,  X->_mat, _mat, SAME_NONZERO_PATTERN);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
#else    // 2.3.x & newer 
    ierr = MatAXPY (_mat, a, X->_mat, DIFFERENT_NONZERO_PATTERN);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
#endif
  }

// ========================================================
  inline double PetscMatrix::operator() (const int i, // index i
                                         const int j  // index j
                                        ) const {

    // do no use this function !!!!!!!!
    assert (this->initialized());
#if PETSC_VERSION_LESS_THAN(2,2,1)  // PETSc 2.2.0 & older
    PetscScalar *petsc_row;
    int* petsc_cols;
#else   // PETSc 2.2.1 & newer
    const PetscScalar *petsc_row;
    const PetscInt    *petsc_cols;
#endif

    // If the entry is not in the sparse matrix, it is 0.
    double value = 0.;
    int  ierr = 0, ncols = 0,  i_val = static_cast<int> (i),    j_val = static_cast<int> (j);

    assert (this->closed()); // the matrix needs to be closed for this to work
    ierr = MatGetRow (_mat, i_val, &ncols, &petsc_cols, &petsc_row);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    // Perform a binary search to find the contiguous index in
    // petsc_cols (resp. petsc_row) corresponding to global index j_val
    std::pair<const int*, const int*> p =
      std::equal_range (&petsc_cols[0], &petsc_cols[0] + ncols, j_val);
    // Found an entry for j_val
    if (p.first != p.second)    {
      // The entry in the contiguous row corresponding to the j_val column of interest
      const int j = std::distance (const_cast<int*> (&petsc_cols[0]), const_cast<int*> (p.first));

      assert (j < ncols);
      assert (petsc_cols[j] == j_val);
      value = static_cast<double> (petsc_row[j]);
    }
    ierr  = MatRestoreRow (_mat, i_val,
                           &ncols, &petsc_cols, &petsc_row);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    return value;
  }


// ========================================================
  inline int PetscMatrix::MatGetRowM (const int i_val, // index i
                                      int cols[],
                                      double vals[]
                                     ) {
    // Set up
    assert (this->initialized());
    const PetscScalar *petsc_row;
    const PetscInt    *petsc_cols;
    int  ierr = 0, ncols = 0;
    // Get row
    assert (this->closed()); // the matrix needs to be closed for this to work
    ierr = MatGetRow (_mat, i_val, &ncols, &petsc_cols, &petsc_row);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    // close row
    ierr  = MatRestoreRow (_mat, i_val, &ncols, &petsc_cols, &petsc_row);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    // print row
    if (&cols[0] != PETSC_NULL)   for (int j = 0; j < ncols; j++)  {
        cols[j] = petsc_cols[j];
        vals[j] = petsc_row[j];
      }
    // return # of columns
    return ncols;
  }


// ================================================
  inline bool PetscMatrix::closed() const {
    assert (this->initialized());
    int ierr = 0;
    PetscBool assembled;
    ierr = MatAssembled (_mat, &assembled);
    CHKERRABORT (MPI_COMM_WORLD, ierr);
    return (assembled == PETSC_TRUE);
  }

// =========================================
  inline void PetscMatrix::swap (
    PetscMatrix &m // pointer to swap
  ) {// =========================================
    std::swap (_mat, m._mat);
    std::swap (_destroy_mat_on_exit, m._destroy_mat_on_exit);
  }

// =========================================================

  inline void PetscMatrix::matrix_add (const double a_in, SparseMatrix &X_in, const char pattern_type []) {
    assert (this->initialized());

    // sanity check. but this cannot avoid
    // crash due to incompatible sparsity structure...
    assert (this->m() == X_in.m());
    assert (this->n() == X_in.n());

    PetscScalar     a = static_cast<PetscScalar> (a_in);
    PetscMatrix* X = dynamic_cast<PetscMatrix*> (&X_in);

    assert (X != NULL);

    int ierr = 0;

    // the matrix from which we copy the values has to be assembled/closed
    // X->close ();
    assert (X->closed());

    semiparallel_only();

    if (!strcmp (pattern_type, "different_nonzero_pattern")) {
      ierr = MatAXPY (_mat, a, X->_mat, DIFFERENT_NONZERO_PATTERN);
    }
    else if (!strcmp (pattern_type, "subset_nonzero_pattern")) {
      ierr = MatAXPY (_mat, a, X->_mat, SUBSET_NONZERO_PATTERN);
    }
    else if (!strcmp (pattern_type, "same_nonzero_pattern")) {
      ierr = MatAXPY (_mat, a, X->_mat, SAME_NONZERO_PATTERN);
    }
    else {
      std::cout << "Error in pattern_type in function PetscMatrix::matrix_add " << std::endl;
      exit (0);
    }
    CHKERRABORT (MPI_COMM_WORLD, ierr);
  }

// =================================================
  inline void PetscMatrix::insert_row (const int row, const int ncols,
                                       const std::vector<int>& cols, double* values) {

    assert (this->initialized());
    int ierr = 0;

    ierr = MatSetValues (_mat, 1, (PetscInt*) &row, (PetscInt) ncols, (PetscInt*) &cols[0],
                         values, INSERT_VALUES);

    CHKERRABORT (MPI_COMM_WORLD, ierr);

  }





} //end namespace femus


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_matrix_h__
