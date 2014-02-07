#ifndef __sparse_rectangular_matrix_h__
#define __sparse_rectangular_matrix_h__

// C++ includes -------------------
#include <memory>
#include <iostream>
#include <vector>

// configure files ----------------
#include "FEMTTUConfig.h"
#include "SolverPackageEnum.hpp"

// forward declarations ----------------------
class SparseMatrix;     // sparse matrix
class DenseMatrix;      // dense matrix
class NumericVector;    // vector

using std::vector;

// ==========================================
/// Rectangular sparse matrix.
// =========================================

class SparseRectangularMatrix {
protected:
  // Data -----------------------------------
  /// Flag indicating whether or not the matrix has been initialized.
  bool _is_initialized;
public:
  // matrix dimensions
  int _m;  ///< number of global rows
  int _n;  ///< number of global columns
  // parallel matrix dimensions
  int _m_l;       ///< number of local rows
  int _n_l;       ///< number of local columns
  int _ml_start;  ///< processor start dof

  // Constructor - Destructor -------------------------------
  /// Constructor;
  SparseRectangularMatrix():_is_initialized(false) {};
  /// Builds a \p SparseRectangularMatrix using the linear solver package specified by \p solver_package
  static std::auto_ptr<SparseRectangularMatrix>  build(const SolverPackage solver_package = LSOLVER);
  // Init functions
  /// Initialize parallel matrix with dims
//     virtual void init(const int  m,  const int  n, const int  m_l,const int  n_l) { _m=m; _n=n; _m=m_l; _n=n_l; }
  virtual void init(const int  m,  const int  n, const int  m_l,const int  n_l,
                    const int  nnz=30,const int  noz=10) {
    _m=m;
    _n=n;
    _m=m_l;
    _n=n_l;
  }
  /// Initialize  matrix  with dims
  virtual void init(const int  m,  const int  n) {
    _m=m;
    _n=n;
  }
  ///Initialize  matrix
  virtual void init() = 0;
  // Destructor
  virtual ~SparseRectangularMatrix() {}; ///< Destructor. Free all memory
  virtual void clear() = 0;  ///< Release all memory and return

  // Setting data -------------------------------------
  // setting single values
  virtual void set(const int i, const int j,const double value) = 0;   ///< Set a value
  virtual void add(const int i, const int j,const double value) = 0;   ///< Add a value
  // setting zero
  virtual void zero() = 0; /// Set all entries to 0.
  virtual void zero_rows(std::vector<int> & rows, double diag_value = 0.0); ///< Set all row entries to 0 and the diagonal entry

  // setting flag
  virtual void close() const = 0;  /// set close flag

// Return data -------------------------------------------------
  // matrix values
  virtual double operator()(const int i,const int j) const = 0;   /// Return the value
//   virtual int MatGetRowM(const int i_val,int *cols,double *vals)=0;      ///< Return the values of a row

  // flag
  virtual bool initialized() const {
    return _is_initialized;  ///< Initialization flag
  }
  virtual bool need_full_sparsity_pattern() const {
    return false;  ///< Full sparsity flag (laspack)
  }
  virtual bool closed() const = 0;                                  ///< Close flag

  // Updates the matrix sparsity pattern
//   virtual void update_sparsity_pattern(const Graph &) =0;                  ///< Full sparsity update
  virtual void update_sparsity_pattern(int m,int n, int m_l,int n_l,
                                       const std::vector<int >  n_oz,
                                       const std::vector<int >  n_nz) =0;  ///< Partial sparsity update

  // dimensions
  virtual int m() const = 0; ///  Row-dimension  m return
  virtual int n() const = 0; ///  Column-dimension n return

  virtual int row_start() const = 0; /// return row_start, the index of the first matrix row
  virtual int row_stop() const = 0;  /// return row_stop, the index of the last matrix row (+1)


  // Operations ---------------------------------
  // Add
  /// Add the full matrix to the Sparse matrix.
  virtual void add_matrix(const DenseMatrix &dm,
                          const std::vector<int> &rows,
                          const std::vector<int> &cols) = 0;
  /// Same, but assumes the row and column maps are the same.
  virtual void add_matrix(const DenseMatrix &dm,const std::vector<int> &dof_indices) = 0;
  /// Add the full matrix by blocks
  virtual void add_matrix_blocked(const std::vector< double > &mat_value,
                                  const std::vector< int> &rows,
                                  const std::vector< int> &cols) = 0;
  
  /// Add a Sparse matrix \p _X, scaled with \p _a, to \p  A += cB
  virtual void add(const double /*c*/, SparseRectangularMatrix & /*B*/) = 0;
  
  /// Add a row to a Sparse matrix
  virtual void insert_row(const int row, const int ncols, const std::vector<int>& cols, const double& values) = 0;
  
  virtual void matrix_add (const double a_in, SparseRectangularMatrix &X_in, const char pattern []) = 0;       
  virtual void matrix_PtAP(const SparseRectangularMatrix &mat_P, const SparseRectangularMatrix &mat_A, const bool &reuse) = 0;	       
  
  // norm of the matrix
  virtual double l1_norm() const = 0;      ///< Return the l1-norm
  virtual double linfty_norm() const = 0;  ///< Return the linfty-norm

  /// Multiplies the matrix with  arg and stores the result in dest.
  void vector_mult(NumericVector& dest,const NumericVector& arg) const;
  /// Multiplies the matrix with  arg and adds the result to  dest.
  void vector_mult_add(NumericVector& dest,const NumericVector& arg) const;

  // form new matrices
  virtual void get_diagonal(NumericVector& dest) const = 0;  ///<  Diagonal
  virtual void get_transpose(SparseRectangularMatrix& dest) const = 0;  ///< Matrix transpose

  // Print - Read  ---------------------------------------------------
  // print
  void print(std::ostream& /*os=std::cout*/) const {};                          ///< Print inot an hdf5 file
  friend std::ostream& operator << (std::ostream& os, const SparseRectangularMatrix& m); ///< Print  to the screen
  virtual void print_personal(std::ostream& os=std::cout) const = 0;            ///< Print  to file
  // print into an hdf5 file
  virtual  void print_hdf5(const std::string name) const=0;                     ///< Print  to hdf5 file
  // read
  void read(const std::string& name);                                 ///< Read
  // Read from an hdf5 file
  virtual void read_len_hdf5(const std::string namefile,const int mode,
                             int len_row[],int leng_off_row[]);      ///< Read lengths
  virtual void read_pos_hdf5(const std::string namefile,const int mode,
                             int pos[],double val[]);                ///< Read pos
  virtual void read_dim_hdf5(const std::string namefile,int dim[]);   ///< Read dimensions


// Submatrices ------------------------------------------------------
  /// SubMatrix creation
  virtual void create_submatrix(SparseRectangularMatrix& submatrix,
                                const std::vector<int>& rows,
                                const std::vector<int>& cols) const {
    this->_get_submatrix(submatrix,rows,cols,false); // false means DO NOT REUSE submatrix
  }

  /// SubMatrix reinitialization (quick init  mode same size)
  virtual void reinit_submatrix(SparseRectangularMatrix& submatrix,
                                const std::vector<int>& rows,
                                const std::vector<int>& cols) const {
    this->_get_submatrix(submatrix,rows,cols,true); // true means REUSE submatrix
  }

protected:

  /// Get submatrix from main matrix
  virtual void _get_submatrix(SparseRectangularMatrix& ,
                              const std::vector<int>& ,
                              const std::vector<int>& ,
                              const bool) const {
    std::cerr << " SparseRectangularMatrix::_get_submatrix not yet implemented \n";
  }
};


// ==============================
/// SparseRectangularMatrix inline members
// ==============================
inline std::ostream& operator << (std::ostream& os, const SparseRectangularMatrix& m) {
  m.print(os);
  return os;
}


#endif
