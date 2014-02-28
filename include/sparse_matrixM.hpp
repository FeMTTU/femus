#ifndef __sparse_matrixM_h__
#define __sparse_matrixM_h__

#include "FemusExtLib_conf.hpp"

#include "Typedefs_conf.hpp"
#include "Graph.hpp"
#include "SolverPackageEnum.hpp"

// #include "linear_solverM.h"

// C++ includes
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>


// Local includes
// #include "libmesh_common.h"
// #include "reference_counted_object.h"

// forward declarations
class DenseMatrixM;
//  inline std::ostream& operator << (std::ostream& os, const SparseMatrixM& m);
// class DofMap;
// namespace SparsityPattern { class Graph; }
class NumericVectorM;


// ==========================================
 // Generic sparse matrix. This class contains
 //  pure virtual members that must be overloaded
 // in derived classes.  Using a derived class
 // allows for uniform access to sparse matrices
 // from various different solver packages in
 // different formats.
 // ===========================================
class SparseMatrixM
#ifdef LM_REFCOUNT
: public ReferenceCountedObject<SparseMatrixM >
#endif

{
protected:
  // ===========================
  // DATA
  // ===========================
   /// The \p DofMap object associated with this object
   //   DofMap const *_dof_map;
  /// Flag indicating whether or not the matrix has been initialized.
  bool _is_initialized;  
    
public:
//   unsigned int _dim;

  // =======================================
  // CONSTRUCTOR
  // =====================================
  /// Constructor; initializes the matrix to be empty, without any structure, i.e.
  /// the matrix is not usable at all. This constructor is therefore only useful
  /// for matrices which are members of a class. All other matrices should be created at a point in the data flow
  /// where all necessary information is available.
  ///You have to initialize the matrix before usage with \p init(...).
  SparseMatrixM();

  /// Destructor. Free all memory
  virtual ~SparseMatrixM ();
  /// Release all memory and return
  virtual void clear () = 0;
  
  /// Builds a \p SparseMatrixM using the linear solver package specified by \p solver_package
  static std::auto_ptr<SparseMatrixM>  build(const SolverPackage solver_package = LSOLVER);
  
  /**
   * Initialize a Sparse matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 10).
   */
  virtual void init (const unsigned int m,
		     const unsigned int n,
		     const unsigned int m_l,
		     const unsigned int n_l,
		     const unsigned int nnz=30,
		     const unsigned int noz=10) = 0;

  /// Initialize using sparsity structure computed by \p dof_map.   
  virtual void init () = 0;
  
  // ===========================
  // SETTING DATA
  // ===========================
  /// Get a pointer to the \p DofMap to use.
//   void attach_dof_map (const DofMap& dof_map) { _dof_map = &dof_map; }
  
  /// Set the element \p (i,j) to \p value.
  virtual void set (const unsigned int i, const unsigned int j,const Real value) = 0;  
  /// Add \p value to the element \p (i,j). 
  virtual void add (const unsigned int i, const unsigned int j,const Real value) = 0;
  
  /// Set all entries to 0.
  virtual void zero () = 0;
  /// Set all row entries to 0 then puts diag_value in the diagonal entry
  virtual void zero_rows (std::vector<int> & rows, Real diag_value = 0.0);
  
  // ===========================
  // RETURN DATA
  // ===========================
 /// @returns true if the matrix has been initialized,
  virtual bool initialized() const { return _is_initialized; }
  
  /// \p returns true if this sparse matrix format needs to be fed the 
  /// graph of the sparse matrix.  This is true in the case of the \p LaspackMatrix, 
  /// but not for the \p PetscMatrix.  In the case where the full graph is not
  /// required we can efficiently approximate it to provide a good estimate of the 
  /// required size of the sparse matrix. 
  virtual bool need_full_sparsity_pattern() const { return false; }  

  /// Updates the matrix sparsity pattern. When your \p SparseMatrixM
  /// implementation does not need this data simply do not overload this method.
  virtual void update_sparsity_pattern (const Graph &) =0;
  
  /// Call the Sparse assemble routines.
  virtual void close () const = 0;
  
  /// @returns \p m, the row-dimension  
  virtual unsigned int m () const = 0;
  /// @returns \p n, the column-dimension   
  virtual unsigned int n () const = 0;

  /// return row_start, the index of the first matrix row 
  virtual unsigned int row_start () const = 0;
  /// return row_stop, the index of the last matrix row (+1) 
  virtual unsigned int row_stop () const = 0;

    /// Return the value of the entry \p (i,j) do not use
  virtual Real operator () (const unsigned int i,const unsigned int j) const = 0;
  
  // ===========================
  // OPERATIONS ON DATA
  // ===========================

  /// Add the full matrix to the Sparse matrix.  
  virtual void add_matrix (const DenseMatrixM &dm,
			   const std::vector<unsigned int> &rows,
			   const std::vector<unsigned int> &cols) = 0;
  
  /// Same, but assumes the row and column maps are the same.
  virtual void add_matrix (const DenseMatrixM &dm,const std::vector<unsigned int> &dof_indices) = 0;
      
  /// Add a Sparse matrix \p _X, scaled with \p _a, to \p  A += cB
  virtual void add (const Real /*c*/, SparseMatrixM & /*B*/) = 0;

  /// Return the l1-norm of the matrix, that is
  virtual Real l1_norm () const = 0;
  /// Return the linfty-norm of the matrix
  virtual Real linfty_norm () const = 0;

  /// see if Sparse matrix has been closed and fully assembled yet
  virtual bool closed() const = 0;
  
  /// Multiplies the matrix with \p arg and stores the result in \p dest.
  void vector_mult (NumericVectorM& dest,
		    const NumericVectorM& arg) const;
  ///Multiplies the matrix with \p arg and adds the result to \p dest.
  void vector_mult_add (NumericVectorM& dest,
			const NumericVectorM& arg) const;

  /// Copies the diagonal part of the matrix into \p dest.
  virtual void get_diagonal (NumericVectorM& dest) const = 0;
  /// Copies the transpose of the matrix into \p dest, which may be this.
  virtual void get_transpose (SparseMatrixM& dest) const = 0;

  // ===========================
  // PRINT 
  // ===========================
  /// Read
  void read(const std::string& name);
   virtual void print(const std::string& name)const{
     std::ofstream out(name.c_str());   print(out);}
  /// Print the contents of the matrix to the screen
  virtual  void print(std::ostream& os=std::cout) const{print_personal(os);}
  ///Same as the print method above
  friend std::ostream& operator << (std::ostream& os, const SparseMatrixM& m);
  /// Print the contents of the matrix to the screen
  virtual void print_personal(std::ostream& os=std::cout) const = 0;
  /// Print the contents of the matrix in Matlab's
  virtual void print_hdf5(const std::string name="NULL") const=0;

  virtual void print_graphic(const bool myvalue) const = 0;

// =======================================
// SUBMATRICES
// ========================================

 /// This function creates a matrix called "submatrix" which is defined
 /// by the row and column indices given in the "rows" and "cols" entries.
 /// Currently this operation is only defined for the PetscMatrix type.
  virtual void create_submatrix(SparseMatrixM& submatrix,
				const std::vector<unsigned int>& rows,
				const std::vector<unsigned int>& cols) const
  {
    this->_get_submatrix(submatrix,rows,cols,false); // false means DO NOT REUSE submatrix
  }

  /// This function is similar to the one above, but it allows you to reuse
  /// the existing sparsity pattern of "submatrix" instead of reallocating
  /// it again.  This should hopefully be more efficient if you are frequently
  /// extracting submatrices of the same size.
  virtual void reinit_submatrix(SparseMatrixM& submatrix,
				const std::vector<unsigned int>& rows,
				const std::vector<unsigned int>& cols) const
  {
    this->_get_submatrix(submatrix,rows,cols,true); // true means REUSE submatrix
  }
  
protected:

  /**
   * Protected implementation of the create_submatrix and reinit_submatrix
   * routines.  Note that this function must be redefined in derived classes
   * for it to work properly!
   */
  virtual void _get_submatrix(SparseMatrixM& ,
			      const std::vector<unsigned int>& ,
			      const std::vector<unsigned int>& ,
			      const bool) const
  {
    std::cerr << "Error! This function is not yet implemented in the base class!"
	      << std::endl;
    abort();
  }
  
 
};


// ==============================
// SparseMatrixM inline members
// ==============================
inline SparseMatrixM::SparseMatrixM():/*_dof_map(NULL),*/_is_initialized(false){}

// ==============================
inline SparseMatrixM::~SparseMatrixM (){}

// ==============================
// For SGI MIPSpro this implementation must occur after
// the partial specialization of the print() member.
inline std::ostream& operator << (std::ostream& os, const SparseMatrixM& m){m.print(os);return os;}


#endif // #ifndef __sparse_matrix_h__
