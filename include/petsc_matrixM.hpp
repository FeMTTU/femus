#ifndef __petsc_matrixMn_h__
#define __petsc_matrixMn_h__

#include "FemusExtLib_conf.hpp"

#ifdef FEMUS_HAVE_PETSC

#include "Typedefs_conf.hpp"
#include "Graph.hpp"

// C++ includes
#include <algorithm>

// Local includes
#include "sparse_matrixM.hpp"
#include "parallelM.hpp"
#include "petsc_macroM.hpp"

// Macro to identify and debug functions which should be called in
// parallel on parallel matrices but which may be called in serial on
// serial matrices.  This macro will only be valid inside non-static
// PetscMatrixM methods
#undef semiparallel_onlyM
#ifndef NDEBUG
  #include <cstring>

  #define semiparallel_onlyM() do { if (this->initialized()) { const char *mytype; \
    MatGetType(_mat,&mytype); \
    if (!strcmp(mytype, MATSEQAIJ)) \
      parallel_onlyM(); } } while (0)
#else
  #define semiparallel_onlyM()
#endif

#include "Casts.hpp"



// Forward Declarations
class DenseMatrixM;

/// Petsc include files.
EXTERN_C_FOR_PETSC_BEGIN
# include <petscmat.h>
EXTERN_C_FOR_PETSC_END



/**
 * Petsc matrix. Provides a nice interface to the
 * Petsc C-based data structures for parallel,
 * sparse matrices.
 */
class PetscMatrixM : public SparseMatrixM
{
public:
  /// Constructor; initializes the matrix to be empty, without any structure, i.e.
  /// the matrix is not usable at all. This  constructor is therefore only useful
  /// for matrices which are members of a
  /// class. All other matrices should be
  /// created at a point in the data flow
  /// where all necessary information is available.
  /// You have to initialize the matrix before usage with \p init(...).
  PetscMatrixM();

  /**
   * Constructor.  Creates a PetscMatrixM assuming you already
   * have a valid Mat object.  In this case, m is NOT destroyed
   * by the PetscMatrixM destructor when this object goes out of scope.
   * This allows ownership of m to remain with the original creator,
   * and to simply provide additional functionality with the PetscMatrixM.
   */
  PetscMatrixM(Mat m);
  
  /// Destructor. Free all memory, but do not release the memory of the sparsity
  ~PetscMatrixM();

  /**
   * Initialize a Petsc matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 30).
   */
  void init (const unsigned int m,
	     const unsigned int n,
	     const unsigned int m_l,
	     const unsigned int n_l,
	     const unsigned int nnz=0,
	     const unsigned int noz=0);

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */   
  void init ();
  
  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. 
   */
  void clear ();

  /**
   * Set all entries to 0. This method retains 
   * sparsity structure.
   */
  void zero ();

  /**
   * Set all row entries to 0 then puts diag_value in the diagonal entry
   */
  void zero_rows (std::vector<int> & rows, Real diag_value = 0.0);
  
  /**
   * Call the Petsc assemble routines.
   * sends necessary messages to other
   * processors
   */
  void close () const;
  
  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  unsigned int m () const;

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  unsigned int n () const;

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  unsigned int row_start () const;

  /**
   * return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  unsigned int row_stop () const;

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  void set (const unsigned int i,
	    const unsigned int j,
	    const Real value);
    
  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  void add (const unsigned int i,
	    const unsigned int j,
	    const Real value);

  /**
   * Add the full matrix to the
   * Petsc matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
    
  void add_matrix (const DenseMatrixM &dm,
		   const std::vector<unsigned int> &rows,
		   const std::vector<unsigned int> &cols);	     

  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  void add_matrix (const DenseMatrixM &dm,
		   const std::vector<unsigned int> &dof_indices);	     
      
  /**
   * Add a Sparse matrix \p X, scaled with \p a, to \p this,
   * stores the result in \p this: 
   * \f$\texttt{this} = a*X + \texttt{this} \f$.
   * Use this with caution, the sparse matrices need to have the
   * same nonzero pattern, otherwise \p PETSc will crash! 
   * It is advisable to not only allocate appropriate memory with 
   * \p init() , but also explicitly zero the terms of \p this
   * whenever you add a non-zero value to \p X.  Note: \p X will 
   * be closed, if not already done, before performing any work.
   */
  void add (const Real a, SparseMatrixM &X);
    
  /**
   * Return the value of the entry
   * \p (i,j).  This may be an
   * expensive operation and you
   * should always take care where
   * to call this function.  In
   * order to avoid abuse, this
   * function throws an exception
   * if the required element does
   * not exist in the matrix.
   *
   * In case you want a function
   * that returns zero instead (for
   * entries that are not in the
   * sparsity pattern of the
   * matrix), use the \p el
   * function.
   */
  Real operator () (const unsigned int i,
		 const unsigned int j) const;

  /**
   * Return the l1-norm of the matrix, that is
   * \f$|M|_1=max_{all columns j}\sum_{all 
   * rows i} |M_ij|\f$,
   * (max. sum of columns).
   * This is the
   * natural matrix norm that is compatible
   * to the l1-norm for vectors, i.e.
   * \f$|Mv|_1\leq |M|_1 |v|_1\f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  Real l1_norm () const;

  /**
   * Return the linfty-norm of the
   * matrix, that is
   * \f$|M|_infty=max_{all rows i}\sum_{all 
   * columns j} |M_ij|\f$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the linfty-norm of vectors, i.e.
   * \f$|Mv|_infty \leq |M|_infty |v|_infty\f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  Real linfty_norm () const;

  /**
   * see if Petsc matrix has been closed
   * and fully assembled yet
   */
  bool closed() const;
  
  ///   sparsity patter update
  void update_sparsity_pattern (const Graph &sparsity_pattern);

 
  void print_graphic(const bool myvalue) const;

  /**
   * Print the contents of the matrix to the screen
   * with the PETSc viewer.  This function only allows
   * printing to standard out, this is because we have
   * limited ourselves to one PETSc implementation for
   * writing.
   */
  void print_personal(std::ostream& os=std::cout) const;

  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  void print_personal(const std::string name="NULL") const;
  
  void print_hdf5(const std::string name="NULL") const;
  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (NumericVectorM& dest) const;

  /**
   * Copies the transpose of the matrix into \p dest, which may be
   * *this.
   */
  virtual void get_transpose (SparseMatrixM& dest) const;

  /**
   * Swaps the raw PETSc matrix context pointers.
   */
  void swap (PetscMatrixM &);

  /**
   * Returns the raw PETSc matrix context pointer.  Note this is generally
   * not required in user-level code. Just don't do anything crazy like
   * calling MatDestroy()!
   */
  Mat mat () { assert (_mat != NULL); return _mat; }

protected:

  /**
   * This function either creates or re-initializes
   * a matrix called "submatrix" which is defined
   * by the row and column indices given in the "rows" and "cols" entries.
   * This function is implemented in terms of the MatGetSubMatrix()
   * routine of PETSc.  The boolean reuse_submatrix parameter determines
   * whether or not PETSc will treat "submatrix" as one which has already
   * been used (had memory allocated) or as a new matrix.
   */
  virtual void _get_submatrix(SparseMatrixM& submatrix,
			      const std::vector<unsigned int>& rows,
			      const std::vector<unsigned int>& cols,
			      const bool reuse_submatrix) const;

private:
  
  /**
   * Petsc matrix datatype to store values
   */				      
  Mat _mat;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a PETSc Mat object. 
   */
  bool _destroy_mat_on_exit;
};




//-----------------------------------------------------------------------
// PetscMatrixM inline members
inline PetscMatrixM::PetscMatrixM()  : _destroy_mat_on_exit(true){}
// =================================================================
inline PetscMatrixM::PetscMatrixM(Mat m): _destroy_mat_on_exit(false){
  this->_mat = m;  this->_is_initialized = true;
}
// ===============================================
inline PetscMatrixM::~PetscMatrixM(){  this->clear();}
// ==========================================
inline void PetscMatrixM::close () const{
  parallel_onlyM();

  // BSK - 1/19/2004
  // strictly this check should be OK, but it seems to
  // fail on matrix-free matrices.  Do they falsely
  // state they are assembled?  Check with the developers...
//   if (this->closed())
//     return;
  
  int ierr=0;
 
  ierr = MatAssemblyBegin (_mat, MAT_FINAL_ASSEMBLY);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = MatAssemblyEnd   (_mat, MAT_FINAL_ASSEMBLY);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
}
// ==================================================
inline unsigned int PetscMatrixM::m () const{
  assert (this->initialized());
  int petsc_m=0, petsc_n=0, ierr=0;
  ierr = MatGetSize (_mat, &petsc_m, &petsc_n);  CHKERRABORT(MPI_COMM_WORLD,ierr);

  return static_cast<unsigned int>(petsc_m);
}
// ============================================
inline unsigned int PetscMatrixM::n () const{
  assert (this->initialized());  
  int petsc_m=0, petsc_n=0, ierr=0;
  ierr = MatGetSize (_mat, &petsc_m, &petsc_n);  CHKERRABORT(MPI_COMM_WORLD,ierr);

  return static_cast<unsigned int>(petsc_n);
}
// ==========================================
inline unsigned int PetscMatrixM::row_start () const{
  assert (this->initialized());
  int start=0, stop=0, ierr=0;
  ierr = MatGetOwnershipRange(_mat, &start, &stop);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<unsigned int>(start);
}
// =========================================
inline unsigned int PetscMatrixM::row_stop () const{
  assert (this->initialized());
  int start=0, stop=0, ierr=0;
  ierr = MatGetOwnershipRange(_mat, &start, &stop);
         CHKERRABORT(MPI_COMM_WORLD,ierr);

  return static_cast<unsigned int>(stop);
}



inline void PetscMatrixM::set (const unsigned int i,
			  const unsigned int j,
			  const Real value){  
  assert (this->initialized());
  int ierr=0, i_val=i, j_val=j;
  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val,
		      &petsc_value, INSERT_VALUES);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
}



// =================================================
inline
void PetscMatrixM::add (const unsigned int i,
			  const unsigned int j,
			  const Real value){
  assert (this->initialized());
  
  int ierr=0, i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val,
		      &petsc_value, ADD_VALUES);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
}



// =====================================================
inline void PetscMatrixM::add_matrix(const DenseMatrixM& dm,
				const std::vector<unsigned int>& dof_indices){
  this->add_matrix (dm, dof_indices, dof_indices);
}
// =========================================================
inline void PetscMatrixM::add (const Real a_in, SparseMatrixM &X_in){
  assert (this->initialized());

  // sanity check. but this cannot avoid 
  // crash due to incompatible sparsity structure...
  assert (this->m() == X_in.m());
  assert (this->n() == X_in.n());

  PetscScalar     a = static_cast<PetscScalar>      (a_in);
  PetscMatrixM* X = dynamic_cast<PetscMatrixM*> (&X_in);

  assert (X != NULL);
  
  int ierr=0;

  // the matrix from which we copy the values has to be assembled/closed
  // X->close ();
  assert(X->closed());

  semiparallel_onlyM();

// 2.2.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,3,0)  
  
  ierr = MatAXPY(&a,  X->_mat, _mat, SAME_NONZERO_PATTERN);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
	 
// 2.3.x & newer
#else
  
  ierr = MatAXPY(_mat, a, X->_mat, DIFFERENT_NONZERO_PATTERN);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
	 
#endif
}




// ========================================================
inline Real PetscMatrixM::operator () (const unsigned int i,
			       const unsigned int j) const{
  assert (this->initialized());
#if PETSC_VERSION_LESS_THAN(2,2,1)  // PETSc 2.2.0 & older
  PetscScalar *petsc_row;
  int* petsc_cols;
#else   // PETSc 2.2.1 & newer
  const PetscScalar *petsc_row;
  const PetscInt    *petsc_cols;
#endif
  // If the entry is not in the sparse matrix, it is 0.  
  Real value=0.;  
  int  ierr=0, ncols=0,  i_val=static_cast<int>(i),    j_val=static_cast<int>(j);
  // the matrix needs to be closed for this to work
  // this->close();
  // but closing it is a semiparallel operation; we want operator()
  // to run on one processor.
  assert(this->closed());
  ierr = MatGetRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Perform a binary search to find the contiguous index in
  // petsc_cols (resp. petsc_row) corresponding to global index j_val
  std::pair<const int*, const int*> p =
    std::equal_range (&petsc_cols[0], &petsc_cols[0] + ncols, j_val);
  // Found an entry for j_val
  if (p.first != p.second)    {
      // The entry in the contiguous row corresponding
      // to the j_val column of interest
      const int j = std::distance (const_cast<int*>(&petsc_cols[0]),
				   const_cast<int*>(p.first));
      
      assert (j < ncols);
      assert (petsc_cols[j] == j_val);
      value = static_cast<Real> (petsc_row[j]);
    } 
  ierr  = MatRestoreRow(_mat, i_val,
			&ncols, &petsc_cols, &petsc_row);
          CHKERRABORT(MPI_COMM_WORLD,ierr);
  return value;
}
// ================================================
inline bool PetscMatrixM::closed() const{
  assert (this->initialized());
  int ierr=0;
  PetscBool assembled;
  ierr = MatAssembled(_mat, &assembled);
         CHKERRABORT(MPI_COMM_WORLD,ierr);
  return (assembled == PETSC_TRUE);
}
// =========================================
inline void PetscMatrixM::swap(PetscMatrixM &m){
  std::swap(_mat, m._mat);
  std::swap(_destroy_mat_on_exit, m._destroy_mat_on_exit);
}
// ===================================================
inline void PetscMatrixM::print_personal(std::ostream& os) const{
  assert (this->initialized());
#ifndef NDEBUG
  if (os != std::cout)
    std::cerr << "Warning! PETSc can only print to std::cout!" << std::endl;
#endif
  int ierr=0;
  ierr = MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
          CHKERRABORT(MPI_COMM_WORLD,ierr);
}

inline void PetscMatrixM::print_graphic(const bool myvalue) const {

  if (myvalue) {
      PetscViewer viewer;
      PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);
      MatView(_mat,viewer);
      double ff;
      std::cin>>ff;
  }
      
      return;
  }

#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_matrix_h__
