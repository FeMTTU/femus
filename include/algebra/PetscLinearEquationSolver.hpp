/*=========================================================================

 Program: FEMUS
 Module: PetscLinearEquationSolver
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __PetscLinearEquationSolver_hpp__
#define __PetscLinearEquationSolver_hpp__

#include "FEMTTUConfig.h"

#ifdef HAVE_PETSC

#ifdef HAVE_MPI
  #include <mpi.h> 
#endif

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "LinearEquationSolver.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include "PetscMacro.hpp"


namespace femus {

/**
* This class inherits the abstract class LinearEquationSolver. In this class the solver is implemented using the PETSc package
*/

class PetscLinearEquationSolver : public LinearEquationSolver {

public:
  
  // Constructor --------------------------------------
  /**  Constructor. Initializes Petsc data structures */
  PetscLinearEquationSolver (const unsigned &igrid, mesh *other_mesh);
  
  /// Destructor.
  ~PetscLinearEquationSolver ();

  /// Release all memory and clear data structures.
  void clear ();
  /// Initialize data structures if not done so already.
  void init ();
  /// Initialize data structures if not done so already plus much more
  void init (Mat& Amat, Mat &Pmat);

  void init(Mat& matrix, const bool pc_flag, const bool Schur);

  void init_schur(Mat& matrix);

  void set_tolerances(const double &rtol, const double &atol,
                      const double &divtol, const unsigned &maxits,const unsigned &index);

  void set_num_elem_vanka_block(const unsigned num_elem_vanka_block);
  
  // Solvers ------------------------------------------------------
  // ========================================================
  /// Call the Vanka(Schur) smoother-solver using the PetscLibrary.
  std::pair< int, double> solve(const vector <unsigned> &VankaIndex,
				const short unsigned &NSchurVar,const bool &Schur,
       			       const bool &ksp_clean);
  /// Call the Gmres smoother-solver
  std::pair< int, double> solve( const bool &clean = true);

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
    return _ksp[0];
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
  
private:
  // data ---------------------------------
  PC _pc;      ///< Preconditioner context
  vector <KSP> _ksp;    ///< Krylov subspace context
  vector< PetscReal > _rtol;
  vector< PetscReal > _abstol;
  vector< PetscReal > _dtol;
  vector< PetscInt >  _maxits;
  unsigned _num_elem_vanka_block;
  
  vector< vector <PetscInt> > _indexai;
  
  vector< vector <unsigned> > _Psize;
  bool _indexai_init;
  
  vector <IS> _isA;
  vector <IS> _isB;
   
  vector< vector <PetscInt> > _is_ovl_idx;
  vector< vector <PetscInt> > _is_loc_idx;
  
  vector <IS> _is_ovl;
  vector <IS> _is_loc;
   
  KSP       *_subksp;  							
  PetscInt  _nlocal,_first;   						
  PC        _subpc;   		
    
  Mat _Pmat;
  bool _Pmat_is_initialized;  

  // Setting --------------------------------------------
  ///  Set the user-specified solver stored in \p _solver_type
  void set_petsc_solver_type ();
  void set_petsc_solver_type2 ();

  clock_t BuildIndex(const vector <unsigned> &VankaIndex, const short unsigned &NSchurVar);
  clock_t BuildIndex();
  clock_t BuildAMSIndex(const vector <unsigned> &VankaIndex,const short unsigned &NSchurVar);
  
};

inline PetscLinearEquationSolver::PetscLinearEquationSolver (const unsigned &igrid, mesh* other_msh)
  : LinearEquationSolver(igrid, other_msh) {
        
  if(igrid==0){
    this->_preconditioner_type = MLU_PRECOND;
    this->_solver_type         = PREONLY;
    _num_elem_vanka_block = _msh->el->GetElementNumber();  
  }
  else{
    if(_msh->_nprocs==1) {  
      this->_preconditioner_type = ILU_PRECOND;
    } 
    else {
      this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
    unsigned dim = _msh->GetDimension();
    unsigned base = pow(2,dim);
    unsigned exponent = 5 - dim;
    _num_elem_vanka_block = pow(base,exponent);
  }

  _ksp.resize(2);  
  
  _rtol.resize(2,1.e-8);
  _abstol.resize(2, 1.e-40);
  _dtol.resize(2, 1.e+50);
  _maxits.resize(2,10);
  
  _indexai_init=0;
  
  _Pmat_is_initialized = false;

}

// =============================================
inline PetscLinearEquationSolver::~PetscLinearEquationSolver () {
  this->clear ();
  
  for(unsigned i=0;i<_isA.size();i++){
    ISDestroy(&_isA[i]); 	
  }
  
  for(unsigned i=0;i<_isB.size();i++){
    ISDestroy(&_isB[i]); 	
  }
  
  for(unsigned i=0;i<_is_loc.size();i++){
    ISDestroy(&_is_loc[i]); 	
  }
  
  for(unsigned i=0;i<_is_ovl.size();i++){
    ISDestroy(&_is_ovl[i]); 	
  }
  
}



} //end namespace femus


#endif 
#endif 
