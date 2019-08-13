/*=========================================================================

 Program: FEMuS
 Module: SparseMatrix
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_SparseMatrix_hpp__
#define __femus_algebra_SparseMatrix_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <memory>
#include <iostream>
#include <vector>

// configure files ----------------
#include "FemusConfig.hpp"
#include "SolverPackageEnum.hpp"
#include "Graph.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
  class SparseMatrix;
  class DenseMatrix;
  class NumericVector;

  using std::vector;

  /**
   *             Generic sparse matrix.
  */

  class SparseMatrix {

    public:

      /** Constructor;  before usage call init(...). */
      SparseMatrix() : _is_initialized (false) {}

      /** Destructor */
      virtual ~SparseMatrix () {}

      /** Release all memory and return */
      virtual void clear () = 0;

      /** Builds a \p SparseMatrix using the linear solver package specified by \p solver_package */
      static std::unique_ptr<SparseMatrix>  build (const SolverPackage solver_package = LSOLVER);

      /** Initialize */
      virtual void init (const int  m,  const int  n, const int  m_l, const int  n_l,
                         const int  nnz = 30, const int  noz = 10) {
        _m = m;
        _n = n;
        _m = m_l;
        _n = n_l;
      }
      /** To be Added */
      virtual void init (const  int m, const  int n, const  int m_l, const  int n_l,
                         const std::vector< int > & n_nz, const std::vector< int > & n_oz) = 0;
      /** To be Added */
      virtual void init (const int  m,  const int  n) {
        _m = m; ///< Initialize  matrix  with dims
        _n = n;
      }
      
      virtual void init (const int nr, const int nc, const std::vector < SparseMatrix*> &P) = 0;

      /** To be Added */
      virtual void init () {};

      // Setting data -------------------------------------
      // setting single values

      /** Set a value */
      virtual void set (const  int i, const  int j, const double value) = 0;

      /** Add a value */
      virtual void add (const  int i, const  int j, const double value) = 0;

      /**  Set all entries to 0 */
      virtual void zero () = 0;

      /** Set all row entries to 0 and the diagonal entry */
      virtual void zero_rows (std::vector<int> & rows, double diag_value = 0.0);

      /** set close flag */
      virtual void close () const = 0;

      // Return data -------------------------------------------------
      // matrix values

      /** Return the value */
      virtual double operator () (const  int i, const  int j) const = 0;

      /** Return the row values */
      virtual int MatGetRowM (const int i_val, int *cols = NULL, double *vals = NULL) = 0;

      // flag
      /** Initialization flag */
      virtual bool initialized() const {
        return _is_initialized;
      }

      /**  Full sparsity flag (laspack) */
      virtual bool need_full_sparsity_pattern() const {
        return false;
      }

      /** Close flag */
      virtual bool closed() const = 0;

      /** @deprecated */
      virtual void update_sparsity_pattern_old (const Graph &) = 0;

      /** Full sparsity update */
      virtual void update_sparsity_pattern (const Graph &) = 0;

      /** Partial sparsity update */
      virtual void update_sparsity_pattern (int m, int n, int m_l, int n_l, const std::vector<int >  n_oz, const std::vector<int >  n_nz) = 0;

      // dimensions
      /** Row-dimension  m return */
      virtual  int m () const = 0;

      /** Column-dimension n return */
      virtual  int n () const = 0;

      /** return row_start, the index of the first matrix row */
      virtual  int row_start () const = 0;

      /** return row_stop, the index of the last matrix row (+1) */
      virtual  int row_stop () const = 0;

      // Operations ---------------------------------
      // Add
      /** Add the full matrix to the Sparse matrix. */
      virtual void add_matrix (const DenseMatrix &dm,
                               const std::vector<unsigned int> &rows,
                               const std::vector<unsigned int> &cols) = 0;

      /** Same, but assumes the row and column maps are the same. */
      virtual void add_matrix (const DenseMatrix &dm, const std::vector<unsigned int> &dof_indices) = 0;

      /** Add a row to a Sparse matrix */
      virtual void insert_row (const int row, const int ncols, const std::vector<int>& cols, double* values) = 0;

      /** To be Addded */
      virtual void add_matrix_blocked (const std::vector< double > &mat_value,
                                       const std::vector< int> &rows,
                                       const std::vector< int> &cols) = 0;

      virtual void add_matrix_blocked (const std::vector< double > &mat_value,
                                       const std::vector< unsigned> &rows,
                                       const std::vector< unsigned> &cols) = 0;

      /** To be Added */
      virtual void matrix_set_off_diagonal_values_blocked (const std::vector< int > &index_rows, const std::vector< int > &index_cols, const double &value) = 0;

      /** To be Added */
      virtual void matrix_set_off_diagonal_values_blocked (const std::vector< int > &index_rows, const std::vector< int > &index_cols, const std::vector<double> &value) = 0;

      /** To be Addded */
      virtual void matrix_add (const double a_in, SparseMatrix &X_in, const char pattern []) = 0;

      /** To be Addded */
      virtual void matrix_PtAP (const SparseMatrix &mat_P, const SparseMatrix &mat_A, const bool &reuse) = 0;

      /** To be Addded */
      virtual void matrix_ABC (const SparseMatrix &mat_A, const SparseMatrix &mat_B, const SparseMatrix &mat_C, const bool &reuse) = 0;

      /** To be Addded */
      virtual void matrix_RightMatMult (const SparseMatrix &mat_A) = 0;

      virtual void matrix_LeftMatMult (const SparseMatrix &mat_A) = 0;

      
      /** To be Addded */
      virtual void matrix_get_diagonal_values (const std::vector< int > &index, std::vector<double> &value) const = 0;

      /** To be Addded */
      virtual void matrix_set_diagonal_values (NumericVector& D) = 0;
      
      /** To be Addded */
      virtual void matrix_set_diagonal_values (const std::vector< int > &index, const double &value) = 0;

      /** To be Addded */
      virtual void matrix_set_diagonal_values (const std::vector< int > &index, const std::vector<double> &value) = 0;

      /** Add a Sparse matrix \p _X, scaled with \p _a, to \p  A += cB */
      virtual void add (const double /*c*/, SparseMatrix & /*B*/) = 0;

      // norm of the matrix
      /**  Return the l1-norm */
      virtual double l1_norm () const = 0;

      /** Return the linfty-norm */
      virtual double linfty_norm () const = 0;

      /** Multiplies the matrix with \p arg and stores the result in \p dest. */
      void vector_mult (NumericVector& dest, const NumericVector& arg) const;

      /** Multiplies the matrix with \p arg and adds the result to \p dest. */
      void vector_mult_add (NumericVector& dest, const NumericVector& arg) const;

      // form new matrices
      /** Diagonal */
      virtual void get_diagonal (NumericVector& dest) const = 0;

      /** Matrix transpose */
      virtual void get_transpose (SparseMatrix& dest) const = 0;

      virtual void mat_zero_rows (const std::vector <int> &index, const double &diagonal_value) const = 0;

      // Read - Print ------------------------------
      // print
      /** Print  to file */
      virtual void print (const std::string& name) const;

      /** Print  to Matlab format */
      virtual void print_matlab (const std::string& name, const std::string& format) const {};

      /** Print  to the screen */
      virtual  void print (std::ostream& os = std::cout) const {
        print_personal (os);
      }

      /** Print  to the text file */
      virtual void print_personal (std::ostream& os = std::cout) const = 0;

      /** @deprecated print to hdf5 files */
      virtual void print_hdf5 (const std::string name = "NULL") const = 0;

      // Read
      /** Read */
      void read (const std::string& name);

      /** Read from hdf5 files */
      /** @deprecated Read lengths */
      virtual void read_len_hdf5 (const std::string namefile, const int mode,
                                  int len_row[], int leng_off_row[]);

      /** @deprecated Read pos */
      virtual void read_pos_hdf5 (const std::string namefile, const int mode, int pos[]);

      /** @deprecated Read dimensions */
      virtual void read_dim_hdf5 (const std::string namefile, int dim[]);

      // Submatrices ------------------------------------------------------
      /** SubMatrix creation - false means DO NOT REUSE submatrix */
      virtual void create_submatrix (SparseMatrix& submatrix,
                                     const std::vector< int>& rows,
                                     const std::vector< int>& cols) const {
        this->_get_submatrix (submatrix, rows, cols, false);
      }

      /** SubMatrix reinitialization (quick init  mode same size) -  true means REUSE submatrix */
      virtual void reinit_submatrix (SparseMatrix& submatrix,
                                     const std::vector< int>& rows,
                                     const std::vector< int>& cols) const {
        this->_get_submatrix (submatrix, rows, cols, true);
      }

      /** Get submatrix from main matrix */
      virtual void _get_submatrix (SparseMatrix& ,
                                   const std::vector< int>& ,
                                   const std::vector< int>& ,
                                   const bool) const {
        std::cerr << " SparseMatrix::_get_submatrix not yet implemented \n";
      }

      /** Print  to the screen */
      friend std::ostream& operator << (std::ostream& os, const SparseMatrix& m);

    protected:

      /** matrix dimensions */
      int _m;         ///< number of global rows
      int _n;         ///< number of global columns
      /** parallel matrix dimensions */
      int _m_l;       ///< number of local rows
      int _n_l;       ///< number of local columns
      int _ml_start;  ///< processor start dof

      /** Flag indicating whether or not the matrix has been initialized. */
      bool _is_initialized;

  };

  /**
   * SparseMatrix inline members
   */
  inline std::ostream& operator << (std::ostream& os, const SparseMatrix& m) {
    m.print (os);
    return os;
  }


} //end namespace femus



#endif

