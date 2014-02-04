#ifndef __petsc_linear_solverMaa_h__
#define __petsc_linear_solverMaa_h__

#include "FEMTTUConfig.h"
//#include "SolverlibConf.hpp"

#if HAVE_PETSC == 1

#if HAVE_MPI == 1
  #include <mpi.h> 
#endif

// Local includes
#include "LinearSolver.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include "PetscMacro.hpp"

// Petsc include files.
EXTERN_C_FOR_PETSC_BEGIN
#include <petscksp.h>
EXTERN_C_FOR_PETSC_END

// ==========================================
/// This class provides an interface to PETSc iterative solvers
class PetscLinearSolver : public LinearSolverM {
  // ============================================

private:
  // data ---------------------------------
  PC _pc;      ///< Preconditioner context
  KSP _ksp;    ///< Krylov subspace context
  KSP _ksp2;   ///< Krylov subspace context
  unsigned num_elem_vanka_block_;

  PetscReal _rtol;
  PetscReal _abstol;
  PetscReal _dtol;
  PetscInt  _maxits;

  PetscReal _rtol_;
  PetscReal _abstol_;
  PetscReal _dtol_;
  PetscInt  _maxits_;
  IS isPA;
  IS isPB;
  IS isB;
  IS isCl;   //local indices for the Schur Complement
  IS isDl; 
  Vec Pw, Pr, Ps;     
  Mat PA;  
  Mat PB;
  Vec f,w;
  Vec z,g;
  Mat A,B,Bt,D,C;
  const PetscInt *ind;
  const PetscInt *ind2;
  PetscInt PAmBsize;
  PetscInt PBsize;
  PetscInt PCsize;
  PetscInt PDsize;
  PetscInt Asize;
  PetscInt counterb;
  PetscInt Csize;
  PetscInt Dsize;
  
  clock_t SearchTime, AssemblyTime, SolveTime0, SolveTime1, SolveTime2, UpdateTime, start_time, end_time;
  PetscErrorCode ierr;
  int its_A, its_C, its;
    
  
public:
  // Constructor --------------------------------------
  ///  Constructor. Initializes Petsc data structures
//   PetscLinearSolver (const char infile[], vector < vector < double> > &vt, const double Lref);
// 
//   PetscLinearSolver (const unsigned &igrid,elem *elc);

  PetscLinearSolver (const unsigned &igrid, mesh *other_mesh);
  
  /// Destructor.
  ~PetscLinearSolver ();
  /// Release all memory and clear data structures.
  void clear ();
  /// Initialize data structures if not done so already.
  void init ();
  /// Initialize data structures if not done so already plus much more
  void init (Mat& matrix);

  void init(Mat& matrix, const bool pc_flag, const bool Schur);

  void init_schur(Mat& matrix);


  void set_tolerances(const double rtol, const double atol,
                      const double divtol, const unsigned maxits);

  void set_schur_tolerances(const double rtol, const double atol,
                            const double divtol, const unsigned maxits);

  // Solvers ------------------------------------------------------
  // ========================================================
  /// Call the Vanka(Schur) smoother-solver using the PetscLibrary.
  int Vanka_Smoother(const vector <unsigned> &_SolPdeIndex,const vector <unsigned> &VankaIndex,
                     const short unsigned &NSchurVar,const bool &Schur);


  /// Call the Vanka smoother-solver using the PetscLibrary.
  int Vanka_Smoother(const vector <unsigned> &_SolPdeIndex, const vector <unsigned> &VankaIndex);

  /// Call the Gmres smoother-solver
  std::pair< int, double> solve();

//   /// Call the Petsc solver.  This function calls the method below, using the
//   /// same matrix for the system and preconditioner matrices.
//   std::pair< int, double>   solve (
//     SparseMatrix  &matrix_in, NumericVector &solution_in, NumericVector &rhs_in,
//     const double tol,const  int m_its)  {
//     return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its);
//   }

  /// This method allows you to call a linear solver while specifying
  /// the matrix to use as the (left) preconditioning matrix.  Note
  /// that the linear solver will not compute a preconditioner in this
  /// case, and will instead premultiply by the matrix you provide.
  /// In PETSc, this is accomplished by calling  PCSetType(_pc, PCMAT);
  /// before invoking KSPSolve().  Note: this functionality is not implemented
  /// in the LinearSolver class since there is not a built-in analog to this method for LasPack
  std::pair< int, double>   solve (
    SparseMatrix  &/*matrix*/,SparseMatrix  &/*preconditioner*/,
    NumericVector &/*solution*/,NumericVector &/*rhs*/,
    const double /*tol*/,const  int /*Niter*/);

// Returns -----------------------------------
  /// Returns the raw PETSc preconditioner context pointer.  This allows
  /// you to specify the PCShellSetApply() and PCShellSetSetUp() functions if you desire.
  PC pc() {
    this->init();
    return _pc;
  }
  /// Returns the raw PETSc ksp context pointer.  This is useful if
  /// you are for example setting a custom convergence test with KSPSetConvergenceTest().
  KSP ksp() {
    this->init();
    return _ksp;
  }
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
  void set_petsc_solver_type2 ();

};


/*----------------------- functions ----------------------------------*/
// =================================================
// inline PetscLinearSolver::PetscLinearSolver (const char infile[],
//     vector < vector < double> > &vt, const double Lref)
//   : LinearSolverM(infile,vt,Lref) {
//    
//     this->_preconditioner_type = MLU_PRECOND;
// 
//   _rtol   = 1.e-8;
//   _abstol = 1.e-40;
//   _dtol   = 1.e+50;
//   _maxits = 10;
// 
//   _rtol_   = 1.e-8;
//   _abstol_ = 1.e-40;
//   _dtol_   = 1.e+50;
//   _maxits_ = 10;
// }
// 
// inline PetscLinearSolver::PetscLinearSolver (const unsigned &igrid,elem *elc)
//   : LinearSolverM(igrid,elc) {
//    
//   std::cout << " PetscLin " << std::endl;  
//   if(_msh->_nprocs==1) {  
//    this->_preconditioner_type = ILU_PRECOND;
//   } else {
//     this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
//   }
// 
//   _rtol   = 1.e-8;
//   _abstol = 1.e-40;
//   _dtol   = 1.e+50;
//   _maxits = 10;
// 
//   _rtol_   = 1.e-8;
//   _abstol_ = 1.e-40;
//   _dtol_   = 1.e+50;
//   _maxits_ = 10;
// 
// 
// }

inline PetscLinearSolver::PetscLinearSolver (const unsigned &igrid, mesh* other_msh)
  : LinearSolverM(igrid, other_msh) {
  
  if(igrid==0){
    this->_preconditioner_type = MLU_PRECOND;  
  }
  else{
    std::cout << " PetscLin " << std::endl;  
    if(_msh->_nprocs==1) {  
      this->_preconditioner_type = ILU_PRECOND;
    } 
    else {
      this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
  }

  _rtol   = 1.e-8;
  _abstol = 1.e-40;
  _dtol   = 1.e+50;
  _maxits = 10;

  _rtol_   = 1.e-8;
  _abstol_ = 1.e-40;
  _dtol_   = 1.e+50;
  _maxits_ = 10;


}

// =============================================
inline PetscLinearSolver::~PetscLinearSolver () {
  this->clear ();
}

#endif // #ifdef FEMTTU_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
