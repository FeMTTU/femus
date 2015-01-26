#ifndef __petsc_linear_solverMaa_h__
#define __petsc_linear_solverMaa_h__

#include "FEMTTUConfig.h"

#ifdef HAVE_PETSC

#include "Mesh.hpp"

#include "Typedefs.hpp"

#include "LinearEquationSolver.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include "PetscMacro.hpp"

// Petsc include files. 
EXTERN_C_FOR_PETSC_BEGIN
// #if PETSC_VERSION_LESS_THAN(2,2,0)
// #  include <petscsles.h>
// #else
#  include <petscksp.h>
// #endif
EXTERN_C_FOR_PETSC_END

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed for preconditioning
// 
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
  // Older versions of PETSc do not have the different int typedefs.
  // On 64-bit machines, PetscInt may actually be a long long int.
  // This change occurred in Petsc-2.2.1.
// #if PETSC_VERSION_LESS_THAN(2,2,1)
//   typedef int PetscErrorCode;
//   typedef int PetscInt;
// #endif
#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
  /// This function is called by PETSc to initialize the preconditioner ctx
  PetscErrorCode __libmesh_petsc_preconditioner_setup (void * ctx);
  /// This function is called by PETSc to acctually apply the preconditioner ctx 
  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y);
#else
  PetscErrorCode __libmesh_petsc_preconditioner_setup (PC);
  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC, Vec x, Vec y);
#endif
} // end extern "C"



namespace femus {


// ==========================================
/// This class provides an interface to PETSc iterative solvers 
class PetscLinearSolverM : public LinearEquationSolver
{// ============================================

   private:
     // data ---------------------------------
// #if PETSC_VERSION_LESS_THAN(2,2,0)  // SLES removed from >= PETSc 2.2.0
//   SLES _sles;///< Linear solver context   
// #endif
  PC _pc;  ///< Preconditioner context  
  KSP _ksp;///< Krylov subspace context 
  
public:
  // Constructor --------------------------------------
  ///  Constructor. Initializes Petsc data structures 
  PetscLinearSolverM (const unsigned &igrid, Mesh *other_mesh);
  /// Destructor.  
  ~PetscLinearSolverM ();
  /// Release all memory and clear data structures.
  void clear ();
  /// Initialize data structures if not done so already.
  void init ();
  /// Initialize data structures if not done so already plus much more
  void init (PetscMatrix* matrix);
  
  void set_tolerances(const double &rtol, const double &atol,
                                const double &divtol, const unsigned &maxits);

  void solve(const vector <unsigned> &VankaIndex, const bool &ksp_clean)
  { std::cout << "To be implemented" << std::endl; abort(); }

  // Solvers ------------------------------------------------------
  /// Call the Petsc solver.  This function calls the method below, using the
  /// same matrix for the system and preconditioner matrices.    
  std::pair<unsigned int, double>   solve (
    SparseMatrix  &matrix_in, NumericVector &solution_in, NumericVector &rhs_in,
    const double tol,const unsigned int m_its)  {
    return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its);
  }

  /// This method allows you to call a linear solver while specifying
  /// the matrix to use as the (left) preconditioning matrix.  Note
  /// that the linear solver will not compute a preconditioner in this
  /// case, and will instead premultiply by the matrix you provide.
  /// In PETSc, this is accomplished by calling  PCSetType(_pc, PCMAT);
  /// before invoking KSPSolve().  Note: this functionality is not implemented
  /// in the LinearSolver class since there is not a built-in analog to this method for LasPack 
  std::pair<unsigned int, double>   solve (
         SparseMatrix  &/*matrix*/,SparseMatrix  &/*preconditioner*/,
	 NumericVector &/*solution*/,NumericVector &/*rhs*/,
	 const double /*tol*/,const unsigned int /*Niter*/);  

// This function solves a system whose matrix is a shell matrix.
//   std::pair<unsigned int, double>
//     solve (const ShellMatrix<T>& shell_matrix,
// 	   NumericVector& solution_in,
// 	   NumericVector& rhs_in,
// 	   const double tol,
// 	   const unsigned int m_its);
  
//  /** This function solves a system whose matrix is a shell matrix, but
//   * a sparse matrix is used as preconditioning matrix, this allowing
//   * other preconditioners than JACOBI.
//   */
//   virtual std::pair<unsigned int, double>
//     solve (const ShellMatrix<T>& shell_matrix,
// 	   const SparseMatrix& precond_matrix,
// 	   NumericVector& solution_in,
// 	   NumericVector& rhs_in,
// 	   const double tol,
// 	   const unsigned int m_its);
// void MGSolve(std::vector<SparseMatrix *>  &/*matrix_in*/,   // System Matrix
//	              std::vector<NumericVector *> &/*solution_in*/,// Solution vector
//		      std::vector<NumericVector *> &/*rhs_in*/,     // RHS vector
//		      Matrix */*P*/, Matrix */*R*/,  
//		      const double /*tol*/,const unsigned int /*m_its*/){
// std::cout << "not implemented"; abort(); 
//}	    
// Returns -----------------------------------	      
  /// Returns the raw PETSc preconditioner context pointer.  This allows
  /// you to specify the PCShellSetApply() and PCShellSetSetUp() functions if you desire.
  PC pc() { this->init(); return _pc; }
  /// Returns the raw PETSc ksp context pointer.  This is useful if
  /// you are for example setting a custom convergence test with KSPSetConvergenceTest().
  KSP ksp() { this->init(); return _ksp; }
  /// Fills the input vector with residual norms from the latest iterative solve.
  void get_residual_history(std::vector<double>& hist);
  /// Returns just the initial residual for the solve just
  /// completed with this interface.  Use this method instead
  /// of the one above if you just want the starting residual and not the entire history.
  double get_initial_residual();

  // Print ----------------------------------------
  /// Prints a useful message about why the latest linear solve con(di)verged.
  virtual void print_converged_reason();
  
  // Setting --------------------------------------------
private:
  ///  Set the user-specified solver stored in \p _solver_type
  void set_petsc_solver_type ();
 // /// Internal function if shell matrix mode is used.
 // static PetscErrorCode _petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest);
 // /// Internal function if shell matrix mode is used.
 // static PetscErrorCode _petsc_shell_matrix_get_diagonal(Mat mat, Vec dest);
  
    PetscReal _rtol;
    PetscReal _abstol;
    PetscReal _dtol;
    PetscInt  _maxits;  
  
};


/*----------------------- functions ----------------------------------*/
// =================================================
inline PetscLinearSolverM::PetscLinearSolverM (const unsigned &igrid, Mesh *other_mesh) :
  LinearEquationSolver(igrid,other_mesh) {
  
  int i; MPI_Comm_size(MPI_COMM_WORLD,&i);  //TODO
  if (i == 1) this->_preconditioner_type = ILU_PRECOND;
  else  this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
}
// =============================================
inline PetscLinearSolverM::~PetscLinearSolverM (){this->clear ();}



} //end namespace femus


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
