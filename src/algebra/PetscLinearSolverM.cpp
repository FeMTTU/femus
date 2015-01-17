#include "FEMTTUConfig.h"

#ifdef HAVE_PETSC

#include "Typedefs.hpp"

// Local Includes
#include "PetscMacro.hpp"
#include "PetscLinearSolverM.hpp"
#include "PetscPreconditioner.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"


namespace femus {





// ==========================================================
extern "C" {
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif

#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
  // ------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_setup(void * ctx) {
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    preconditioner->init();
    return 0;
  }
// ------------------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y)  {
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    PetscVector x_vec(x);
    PetscVector y_vec(y);
    preconditioner->apply(x_vec,y_vec);
    return 0;
  }
#else
// se if it compiles...
PetscErrorCode __libmesh_petsc_preconditioner_setup(PC pc) {
  void *ctx;
  PetscErrorCode ierr = PCShellGetContext(pc,&ctx);
  CHKERRQ(ierr);
  Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
  preconditioner->init();
  return 0;
}
// --------------------------------------------------------------------------
PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y) {
  void *ctx;
  PetscErrorCode ierr = PCShellGetContext(pc,&ctx);
  CHKERRQ(ierr);
  Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
  PetscVector x_vec(x);
  PetscVector y_vec(y);
  preconditioner->apply(x_vec,y_vec);
  return 0;
}
#endif
} // end extern "C
// ================================================


// /*extern "C"
// {
// #if PETSC_VERSION_LESS_THAN(2,2,1)
//   typedef int PetscErrorCode;
//   typedef int PetscInt;
// #endif
// 
// 
// #if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
//   PetscErrorCode __libmesh_petsc_preconditioner_setup (void * ctx)
//   {
//     Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
//     preconditioner->init();
// 
//     return 0;
//   }
//   
// 
//   PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y)
//   {
//     Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
// 
//     PetscVector<Number> x_vec(x);
//     PetscVector<Number> y_vec(y);
// 
//     preconditioner->apply(x_vec,y_vec);
// 
//     return 0;
//   }
// #else
//   PetscErrorCode __libmesh_petsc_preconditioner_setup (PC pc)
//   {
//     void *ctx;
//     PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
//     Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
//     preconditioner->init();
// 
//     return 0;
//   }
// 
//   PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y)
//   {
//     void *ctx;
//     PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
//     Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
// 
//     PetscVector<Number> x_vec(x);
//     PetscVector<Number> y_vec(y);
// 
//     preconditioner->apply(x_vec,y_vec);
// 
//     return 0;
//   }
// #endif
// } // end extern "C"*/




  void PetscLinearSolverM::set_tolerances(const double &rtol, const double &atol,
						    const double &divtol, const unsigned &maxits) {

    _rtol   = static_cast<PetscReal>(rtol);
    _abstol = static_cast<PetscReal>(atol);
    _dtol   = static_cast<PetscReal>(divtol);
    _maxits = static_cast<PetscInt>(maxits);
						      
  }




// ==============================================
// ----------------------- functions ------
// ==============================================
void PetscLinearSolverM::clear() {
  if (this->initialized()) {
    this->_is_initialized = false;
    int ierr=0;
// #if PETSC_VERSION_LESS_THAN(2,2,0) // 2.1.x & earlier style 
//     ierr = SLESDestroy(_sles); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else    // 2.2.0 & newer style
    ierr = KSPDestroy(&_ksp);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
    // Mimic PETSc default solver and preconditioner
    this->_solver_type  = GMRES;
    if (!this->_preconditioner)      {
      int i; MPI_Comm_size(MPI_COMM_WORLD,&i);
      if (i == 1) this->_preconditioner_type = ILU_PRECOND;
      else this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
  }
}

// ==============================================================
void PetscLinearSolverM::init() {
  // Initialize the data structures if not done so already.
  if (!this->initialized()) {
    this->_is_initialized = true;
    int ierr=0;
// 2.1.x & earlier style
// #if PETSC_VERSION_LESS_THAN(2,2,0)
//     // Create the linear solver context
//     ierr = SLESCreate(MPI_COMM_WORLD, &_sles);  CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Create the Krylov subspace & preconditioner contexts
//     ierr = SLESGetKSP(_sles, &_ksp);   CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = SLESGetPC(_sles, &_pc);    CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Have the Krylov subspace method use our good initial guess rather than 0
//     ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);  CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Set user-specified  solver and preconditioner types
//     this->set_petsc_solver_type();
//
//     // Set the options from user-input
//     // Set runtime options, e.g.,
//     //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//     //  These options will override those specified above as long as
//     //  SLESSetFromOptions() is called _after_ any other customization
//     //  routines.
//
//     ierr = SLESSetFromOptions(_sles);   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
// // 2.2.0 & newer style
// #else
    // Create the linear solver context
    ierr = KSPCreate(MPI_COMM_WORLD, &_ksp);   CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Create the preconditioner context
    ierr = KSPGetPC(_ksp, &_pc); CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Have the Krylov subspace method use our good initial guess rather than 0
    ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set user-specified  solver and preconditioner types
    this->set_petsc_solver_type();
    // Set the options from user-input
    // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    //  These options will override those specified above as long as
    //  KSPSetFromOptions() is called _after_ any other customization  routines.
    ierr = KSPSetFromOptions(_ksp);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
    // NOT NECESSARY!!!!
    //ierr = PCSetFromOptions (_pc);
    //CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif

    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
    ierr = KSPSetResidualHistory(_ksp,
                                 PETSC_NULL,   // pointer to the array which holds the history
                                 PETSC_DECIDE, // size of the array holding the history
                                 PETSC_TRUE);  // Whether or not to reset the history for each solve.
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);

    //If there is a preconditioner object we need to set the internal setup and apply routines
    if (this->_preconditioner) {
      PCShellSetContext(_pc,(void*)this->_preconditioner);
      PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
      PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
    }
  }
}

// ========================================================
void PetscLinearSolverM::init(PetscMatrix* matrix) {
  // Initialize the data structures if not done so already.
  if (!this->initialized())    {
    this->_is_initialized = true;   int ierr=0;
// #if PETSC_VERSION_LESS_THAN(2,2,0)  // 2.1.x & earlier style
//     // Create the linear solver context
//     ierr = SLESCreate(MPI_COMM_WORLD, &_sles);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Create the Krylov subspace & preconditioner contexts
//     ierr = SLESGetKSP(_sles, &_ksp);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = SLESGetPC(_sles, &_pc);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Have the Krylov subspace method use our good initial guess rather than 0
//     ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Set user-specified  solver and preconditioner types
//     this->set_petsc_solver_type();
//     // Set the options from user-input
//     // Set runtime options, e.g.,
//     //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//     //  These options will override those specified above as long as
//     //  SLESSetFromOptions() is called _after_ any other customization
//     //  routines.
//     ierr = SLESSetFromOptions(_sles);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//
// #else // 2.2.0 & newer style
    // Create the linear solver context
    ierr = KSPCreate(MPI_COMM_WORLD, &_ksp);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    //ierr = PCCreate (MPI_COMM_WORLD, &_pc); CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Create the preconditioner context
    ierr = KSPGetPC(_ksp, &_pc);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set operators. The input matrix works as the preconditioning matrix
    //ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat(),SAME_NONZERO_PATTERN);
    ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat()); //PETSC3p5
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Have the Krylov subspace method use our good initial guess rather than 0
    ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set user-specified  solver and preconditioner types
    this->set_petsc_solver_type();
    // Set the options from user-input
    // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    //  These options will override those specified above as long as
    //  KSPSetFromOptions() is called _after_ any other customization  routines.
    ierr = KSPSetFromOptions(_ksp);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
    //ierr = PCSetFromOptions (_pc);CHKERRABORT(MPI_COMM_WORLD,ierr);

// #endif

    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
    ierr = KSPSetResidualHistory(_ksp,
                                 PETSC_NULL,   // pointer to the array which holds the history
                                 PETSC_DECIDE, // size of the array holding the history
                                 PETSC_TRUE);  // Whether or not to reset the history for each solve.
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
    if (this->_preconditioner) {
      this->_preconditioner->set_matrix(*matrix);
      PCShellSetContext(_pc,(void*)this->_preconditioner);
      PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
      PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
    }
  }
}



// /*void PetscLinearSolverM::init ()
// {
//   // Initialize the data structures if not done so already.
//   if (!this->initialized())
//     {
//       this->_is_initialized = true;
//       
//       int ierr=0;
//   
// // 2.1.x & earlier style      
// #if PETSC_VERSION_LESS_THAN(2,2,0)
//       
//       // Create the linear solver context
//       ierr = SLESCreate (MPI_COMM_WORLD, &_sles);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       // Create the Krylov subspace & preconditioner contexts
//       ierr = SLESGetKSP       (_sles, &_ksp);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       ierr = SLESGetPC        (_sles, &_pc);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       // Set user-specified  solver and preconditioner types
//       this->set_petsc_solver_type();
//       
//       // Set the options from user-input
//       // Set runtime options, e.g.,
//       //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//       //  These options will override those specified above as long as
//       //  SLESSetFromOptions() is called _after_ any other customization
//       //  routines.
//       
//       ierr = SLESSetFromOptions (_sles);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
// // 2.2.0 & newer style
// #else
//       
//       // Create the linear solver context
//       ierr = KSPCreate (MPI_COMM_WORLD, &_ksp);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       // Create the preconditioner context
//       ierr = KSPGetPC        (_ksp, &_pc);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//       // Set user-specified  solver and preconditioner types
//       this->set_petsc_solver_type();
// 
//       // Set the options from user-input
//       // Set runtime options, e.g.,
//       //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//       //  These options will override those specified above as long as
//       //  KSPSetFromOptions() is called _after_ any other customization
//       //  routines.
//       
//       ierr = KSPSetFromOptions (_ksp);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//       // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
//       // NOT NECESSARY!!!!
//       //ierr = PCSetFromOptions (_pc);
//       //CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
// 	       
// #endif
// 	     
//       // Have the Krylov subspace method use our good initial guess
//       // rather than 0, unless the user requested a KSPType of
//       // preonly, which complains if asked to use initial guesses.
// #if PETSC_VERSION_LESS_THAN(3,0,0)
//       KSPType ksp_type;
// #else
//       const KSPType ksp_type;
// #endif
// 
//       ierr = KSPGetType (_ksp, &ksp_type);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//       if (strcmp(ksp_type, "preonly"))
//         {
//           ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
//                  CHKERRABORT(MPI_COMM_WORLD,ierr);
//         }
//       
//       // Notify PETSc of location to store residual history.
//       // This needs to be called before any solves, since
//       // it sets the residual history length to zero.  The default
//       // behavior is for PETSc to allocate (internally) an array
//       // of size 1000 to hold the residual norm history.
//       ierr = KSPSetResidualHistory(_ksp,
// 				   PETSC_NULL,   // pointer to the array which holds the history
// 				   PETSC_DECIDE, // size of the array holding the history
// 				   PETSC_TRUE);  // Whether or not to reset the history for each solve. 
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//       PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
// 
//       //If there is a preconditioner object we need to set the internal setup and apply routines
//       if(this->_preconditioner)
//       {
//         PCShellSetContext(_pc,(void*)this->_preconditioner);
//         PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
//         PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
//       }
//     }
// }
// 
// 
// template <typename T>
// void PetscLinearSolver<T>::init ( PetscMatrix<T>* matrix )
// {
//   // Initialize the data structures if not done so already.
//   if (!this->initialized())
//     {
//       this->_is_initialized = true;
//       
//       int ierr=0;
// 
// // 2.1.x & earlier style      
// #if PETSC_VERSION_LESS_THAN(2,2,0)
//       
//       // Create the linear solver context
//       ierr = SLESCreate (MPI_COMM_WORLD, &_sles);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       // Create the Krylov subspace & preconditioner contexts
//       ierr = SLESGetKSP       (_sles, &_ksp);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       ierr = SLESGetPC        (_sles, &_pc);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       // Set user-specified  solver and preconditioner types
//       this->set_petsc_solver_type();
//       
//       // Set the options from user-input
//       // Set runtime options, e.g.,
//       //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//       //  These options will override those specified above as long as
//       //  SLESSetFromOptions() is called _after_ any other customization
//       //  routines.
//       
//       ierr = SLESSetFromOptions (_sles);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
// // 2.2.0 & newer style
// #else
//       
//       // Create the linear solver context
//       ierr = KSPCreate (MPI_COMM_WORLD, &_ksp);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//     
//       //ierr = PCCreate (MPI_COMM_WORLD, &_pc);
//       //     CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       // Create the preconditioner context
//       ierr = KSPGetPC        (_ksp, &_pc);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       // Set operators. The input matrix works as the preconditioning matrix
//       ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat(),SAME_NONZERO_PATTERN);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);       
// 
//       // Set user-specified  solver and preconditioner types
//       this->set_petsc_solver_type();
// 
//       // Set the options from user-input
//       // Set runtime options, e.g.,
//       //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//       //  These options will override those specified above as long as
//       //  KSPSetFromOptions() is called _after_ any other customization
//       //  routines.
//       
//       ierr = KSPSetFromOptions (_ksp);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//       // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
//       // NOT NECESSARY!!!!
//       //ierr = PCSetFromOptions (_pc);
//       //CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
// 	       
// #endif
// 	     
//       // Have the Krylov subspace method use our good initial guess
//       // rather than 0, unless the user requested a KSPType of
//       // preonly, which complains if asked to use initial guesses.
// #if PETSC_VERSION_LESS_THAN(3,0,0)
//       KSPType ksp_type;
// #else
//       const KSPType ksp_type;
// #endif
// 
//       ierr = KSPGetType (_ksp, &ksp_type);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//       if (strcmp(ksp_type, "preonly"))
//         {
//           ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);
//                  CHKERRABORT(MPI_COMM_WORLD,ierr);
//         }
//       
//       // Notify PETSc of location to store residual history.
//       // This needs to be called before any solves, since
//       // it sets the residual history length to zero.  The default
//       // behavior is for PETSc to allocate (internally) an array
//       // of size 1000 to hold the residual norm history.
//       ierr = KSPSetResidualHistory(_ksp,
// 				   PETSC_NULL,   // pointer to the array which holds the history
// 				   PETSC_DECIDE, // size of the array holding the history
// 				   PETSC_TRUE);  // Whether or not to reset the history for each solve. 
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//       PetscPreconditioner<T>::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);      
//       if(this->_preconditioner)
//       {
//         this->_preconditioner->set_matrix(*matrix);
//         PCShellSetContext(_pc,(void*)this->_preconditioner);
//         PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
//         PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
//       }
//     }
// }*/











// ========================================================
std::pair<unsigned int, double> PetscLinearSolverM::solve(SparseMatrix&  matrix_in,
    SparseMatrix&  precond_in,  NumericVector& solution_in,  NumericVector& rhs_in,
    const double tol,   const unsigned int m_its) {

//   START_LOG("solve()", "PetscLinearSolverM");
  // Make sure the data passed in are really of Petsc types
  PetscMatrix* matrix   = libmeshM_cast_ptr<PetscMatrix*>(&matrix_in);
  PetscMatrix* precond  = libmeshM_cast_ptr<PetscMatrix*>(&precond_in);
  PetscVector* solution = libmeshM_cast_ptr<PetscVector*>(&solution_in);
  PetscVector* rhs      = libmeshM_cast_ptr<PetscVector*>(&rhs_in);
  this->init(matrix);

  int ierr=0;  int its=0, max_its = static_cast<int>(m_its);
  PetscReal final_resid=0.;
  // Close the matrices and vectors in case this wasn't already done.
  matrix->close();  precond->close();  solution->close(); rhs->close();
//   // If matrix != precond, then this means we have specified a
//   // special preconditioner, so reset preconditioner type to PCMAT.
//   if (matrix != precond)
//     {
//       this->_preconditioner_type = USER_PRECOND;
//       this->set_petsc_preconditioner_type ();
//     }
  if (this->_preconditioner) this->_preconditioner->set_matrix(matrix_in);
  // 2.2.1 & newer style
  // Set operators. The input matrix works as the preconditioning matrix
  if (!this->same_preconditioner)  {
    //ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),SAME_NONZERO_PATTERN);
    ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat()); //PETSC3p5
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else  {
    //ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),SAME_PRECONDITIONER);
    ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat()); //PETSC3p5
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances(_ksp, tol, PETSC_DEFAULT,PETSC_DEFAULT, max_its);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Solve the linear system
  
//        PetscLogEvent USER_EVENT;
//      PetscLogDouble user_event_flops;
//      PetscLogEventRegister("User event",0,&USER_EVENT);
//      PetscLogEventBegin(USER_EVENT,0,0,0,0);
       
 
  ierr = KSPSolve(_ksp, rhs->vec(), solution->vec()); CHKERRABORT(MPI_COMM_WORLD,ierr);
//         PetscLogFlops(user_event_flops);
//      PetscLogEventEnd(USER_EVENT,0,0,0,0);

     // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber(_ksp, &its); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm(_ksp, &final_resid); CHKERRABORT(MPI_COMM_WORLD,ierr);

//   STOP_LOG("solve()", "PetscLinearSolverM");
  return std::make_pair(its, final_resid);
}

// std::pair<unsigned int, double>
// PetscLinearSolverM::solve (const ShellMatrix<T>& shell_matrix,
// 			     NumericVector<T>& solution_in,
// 			     NumericVector<T>& rhs_in,
// 			     const double tol,
// 			     const unsigned int m_its)
// {
//
// #if PETSC_VERSION_LESS_THAN(2,3,1)
//   // FIXME[JWP]: There will be a bunch of unused variable warnings
//   // for older PETScs here.
//   std::cout << "This method has been developed with PETSc 2.3.1.  "
// 	    << "No one has made it backwards compatible with older "
// 	    << "versions of PETSc so far; however, it might work "
// 	    << "without any change with some older version." << std::endl;
//   libmesh_error();
//   return std::make_pair(0,0.0);
//
// #else
//
//   START_LOG("solve()", "PetscLinearSolver");
//
//   // Make sure the data passed in are really of Petsc types
//   PetscVector* solution = libmesh_cast_ptr<PetscVector*>(&solution_in);
//   PetscVector* rhs      = libmesh_cast_ptr<PetscVector*>(&rhs_in);
//
//   this->init ();
//
//   int ierr=0;
//   int its=0, max_its = static_cast<int>(m_its);
//   PetscReal final_resid=0.;
//
//   // Close the matrices and vectors in case this wasn't already done.
//   solution->close ();
//   rhs->close ();
//
//   // Prepare the matrix.
//   Mat mat;
//   ierr = MatCreateShell(MPI_COMM_WORLD,
// 			rhs_in.local_size(),
// 			solution_in.local_size(),
// 			rhs_in.size(),
// 			solution_in.size(),
// 			const_cast<void*>(static_cast<const void*>(&shell_matrix)),
// 			&mat);
//   /* Note that the const_cast above is only necessary because PETSc
//      does not accept a const void*.  Inside the member function
//      _petsc_shell_matrix() below, the pointer is casted back to a
//      const ShellMatrix<T>*.  */
//
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
//   ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Set operators. The input matrix works as the preconditioning matrix
//   ierr = KSPSetOperators(_ksp, mat, mat,
// 			 SAME_NONZERO_PATTERN);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Set the tolerances for the iterative solver.  Use the user-supplied
//   // tolerance for the relative residual & leave the others at default values.
//   ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
//  			   PETSC_DEFAULT, max_its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Solve the linear system
//   ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Get the number of iterations required for convergence
//   ierr = KSPGetIterationNumber (_ksp, &its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Get the norm of the final residual to return to the user.
//   ierr = KSPGetResidualNorm (_ksp, &final_resid);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Destroy the matrix.
//   ierr = MatDestroy(mat);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   STOP_LOG("solve()", "PetscLinearSolver");
//   // return the # of its. and the final residual norm.
//   return std::make_pair(its, final_resid);
//
// #endif
//
// }


// std::pair<unsigned int, double>
// PetscLinearSolverM::solve (const ShellMatrix<T>& shell_matrix,
// 			     const SparseMatrix<T>& precond_matrix,
// 			     NumericVector<T> &solution_in,
// 			     NumericVector<T> &rhs_in,
// 			     const double tol,
// 			     const unsigned int m_its)
// {
//
// #if PETSC_VERSION_LESS_THAN(2,3,1)
//   // FIXME[JWP]: There will be a bunch of unused variable warnings
//   // for older PETScs here.
//   std::cout << "This method has been developed with PETSc 2.3.1.  "
// 	    << "No one has made it backwards compatible with older "
// 	    << "versions of PETSc so far; however, it might work "
// 	    << "without any change with some older version." << std::endl;
//   libmesh_error();
//   return std::make_pair(0,0.0);
//
// #else
//
//   START_LOG("solve()", "PetscLinearSolver");
//
//   // Make sure the data passed in are really of Petsc types
//   const PetscMatrix* precond  = libmesh_cast_ptr<const PetscMatrix*>(&precond_matrix);
//   PetscVector* solution = libmesh_cast_ptr<PetscVector*>(&solution_in);
//   PetscVector* rhs      = libmesh_cast_ptr<PetscVector*>(&rhs_in);
//
//   this->init ();
//
//   int ierr=0;
//   int its=0, max_its = static_cast<int>(m_its);
//   PetscReal final_resid=0.;
//
//   // Close the matrices and vectors in case this wasn't already done.
//   solution->close ();
//   rhs->close ();
//
//   // Prepare the matrix.
//   Mat mat;
//   ierr = MatCreateShell(MPI_COMM_WORLD,
// 			rhs_in.local_size(),
// 			solution_in.local_size(),
// 			rhs_in.size(),
// 			solution_in.size(),
// 			const_cast<void*>(static_cast<const void*>(&shell_matrix)),
// 			&mat);
//   /* Note that the const_cast above is only necessary because PETSc
//      does not accept a const void*.  Inside the member function
//      _petsc_shell_matrix() below, the pointer is casted back to a
//      const ShellMatrix<T>*.  */
//
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
//   ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Set operators. The input matrix works as the preconditioning matrix
//   ierr = KSPSetOperators(_ksp, mat, const_cast<PetscMatrix*>(precond)->mat(),
// 			 DIFFERENT_NONZERO_PATTERN);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   if(this->_preconditioner)
//     this->_preconditioner->set_matrix(const_cast<SparseMatrix<Number>&>(precond_matrix));
//
//   // Set the tolerances for the iterative solver.  Use the user-supplied
//   // tolerance for the relative residual & leave the others at default values.
//   ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
//  			   PETSC_DEFAULT, max_its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Solve the linear system
//   ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Get the number of iterations required for convergence
//   ierr = KSPGetIterationNumber (_ksp, &its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Get the norm of the final residual to return to the user.
//   ierr = KSPGetResidualNorm (_ksp, &final_resid);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Destroy the matrix.
//   ierr = MatDestroy(mat);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   STOP_LOG("solve()", "PetscLinearSolver");
//   // return the # of its. and the final residual norm.
//   return std::make_pair(its, final_resid);
//
// #endif
//
// }


// =========================================================================
void PetscLinearSolverM::get_residual_history(std::vector<double>& hist) {
  int ierr = 0;  int its  = 0;
  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p; ierr = KSPGetResidualHistory(_ksp, &p, &its); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Check for early return
  if (its == 0) return;
  // Create space to store the result
  hist.resize(its);
  // Copy history into the vector provided by the user.
  for (int i=0; i<its; ++i) { hist[i] = *p; p++; }
}

// ======================================================
double PetscLinearSolverM::get_initial_residual() {
  int ierr = 0;  int its  = 0;
  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p; ierr = KSPGetResidualHistory(_ksp, &p, &its); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Check no residual history
  if (its == 0) {std::cerr << "No iterations have been performed, returning 0." << std::endl; return 0.; }
  // Otherwise, return the value pointed to by p.
  return *p;
}

// =================================================
void PetscLinearSolverM::set_petsc_solver_type() {
  int ierr = 0;
  switch (this->_solver_type) {
  case CG:
    ierr = KSPSetType(_ksp, (char*) KSPCG);         CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case CR:
    ierr = KSPSetType(_ksp, (char*) KSPCR);         CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case CGS:
    ierr = KSPSetType(_ksp, (char*) KSPCGS);        CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case BICG:
    ierr = KSPSetType(_ksp, (char*) KSPBICG);       CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case TCQMR:
    ierr = KSPSetType(_ksp, (char*) KSPTCQMR);      CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case TFQMR:
    ierr = KSPSetType(_ksp, (char*) KSPTFQMR);      CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case LSQR:
    ierr = KSPSetType(_ksp, (char*) KSPLSQR);       CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case BICGSTAB:
    ierr = KSPSetType(_ksp, (char*) KSPBCGS);       CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case MINRES:
    ierr = KSPSetType(_ksp, (char*) KSPMINRES);     CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case GMRES:
    ierr = KSPSetType(_ksp, (char*) KSPGMRES);      CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case RICHARDSON:
    ierr = KSPSetType(_ksp, (char*) KSPRICHARDSON); CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case CHEBYSHEV:
    ierr = KSPSetType(_ksp, (char*) KSPCHEBYSHEV);  CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  default:
    std::cerr << "ERROR:  Unsupported PETSC Solver: "
              << this->_solver_type               << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }
}

// =======================================================
void PetscLinearSolverM::print_converged_reason() {

   KSPConvergedReason reason;
  KSPGetConvergedReason(_ksp, &reason);
  std::cout << "Linear solver convergence/divergence reason: " << KSPConvergedReasons[reason] << std::endl;

}

// PetscErrorCode PetscLinearSolverM::_petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest)
// {
//   /* Get the matrix context.  */
//   int ierr=0;
//   void* ctx;
//   ierr = MatShellGetContext(mat,&ctx);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   /* Get user shell matrix object.  */
//   const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);
//
//   /* Make \p NumericVector instances around the vectors.  */
//   PetscVector arg_global(arg);
//   PetscVector dest_global(dest);
//
//   /* Call the user function.  */
//   shell_matrix.vector_mult(dest_global,arg_global);
//
//   return ierr;
// }

// =====================================================
// PetscErrorCode PetscLinearSolverM::_petsc_shell_matrix_get_diagonal(Mat mat, Vec dest){
//   /* Get the matrix context.  */
//   int ierr=0;
//   void* ctx;
//   ierr = MatShellGetContext(mat,&ctx);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   /* Get user shell matrix object.  */
//   const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);
//
//   /* Make \p NumericVector instances around the vector.  */
//   PetscVector dest_global(dest);
//
//   /* Call the user function.  */
//   shell_matrix.get_diagonal(dest_global);
//
//   return ierr;
// }



//------------------------------------------------------------------
// Explicit instantiations
// template class PetscLinearSolver<Number>;




} //end namespace femus



#endif // #ifdef LIBMESH_HAVE_PETSC
