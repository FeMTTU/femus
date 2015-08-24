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

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

// Local Includes
#include "PetscMacro.hpp"
#include "GmresPetscLinearEquationSolver.hpp"
#include "PetscPreconditioner.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include <iomanip>
#include <sstream>

namespace femus {

  using namespace std;

// ==============================================
// ----------------------- functions ------
// ==============================================

  void GmresPetscLinearEquationSolver::set_tolerances(const double& rtol, const double& atol,
      const double& divtol, const unsigned& maxits) {

    _rtol   = static_cast<PetscReal>(rtol);
    _abstol = static_cast<PetscReal>(atol);
    _dtol   = static_cast<PetscReal>(divtol);
    _maxits = static_cast<PetscInt>(maxits);

  };

// ================================================

  clock_t GmresPetscLinearEquationSolver::BuildIndex(const vector <unsigned>& variable_to_be_solved) {

    clock_t SearchTime = 0;
    clock_t start_time = clock();
    _indexai_init = 1;

    unsigned IndexaSize = KKoffset[KKIndex.size() - 1][processor_id()] - KKoffset[0][processor_id()];
    _indexai.resize(2);
    _indexai[0].resize(IndexaSize);
    _indexai[1].resize(IndexaSize);

    vector <bool> ThisSolutionIsIncluded(_SolPdeIndex.size(), false);
    for (unsigned iind = 0; iind < variable_to_be_solved.size(); iind++) {
      unsigned PdeIndexSol = variable_to_be_solved[iind];
      ThisSolutionIsIncluded[PdeIndexSol] = true;
    }

    unsigned count0 = 0;
    unsigned count1 = 0;
    for (int k = 0; k < _SolPdeIndex.size(); k++) {
      unsigned indexSol = _SolPdeIndex[k];
      unsigned soltype = _SolType[indexSol];
      for (unsigned inode_mts = _msh->MetisOffset[soltype][processor_id()];
           inode_mts < _msh->MetisOffset[soltype][processor_id() + 1]; inode_mts++) {
        int local_mts = inode_mts - _msh->MetisOffset[soltype][processor_id()];
        int idof_kk = KKoffset[k][processor_id()] + local_mts;
        if (!ThisSolutionIsIncluded[k] || (* (*_Bdc) [indexSol])(inode_mts) < 1.9) {
          _indexai[0][count0] = idof_kk;
          count0++;
        } else {
          _indexai[1][count1] = idof_kk;
          count1++;
        }
      }
    }
    _indexai[0].resize(count0);
    _indexai[1].resize(count1);

    std::sort(_indexai[0].begin(), _indexai[0].end());
    std::sort(_indexai[1].begin(), _indexai[1].end());



    //BEGIN Generate std::vector<IS> for GMRES solve by elimination ***********
    _isA.resize(1);
    PetscInt Asize = _indexai[1].size();
    int ierr = ISCreateGeneral(MPI_COMM_WORLD, Asize, &_indexai[1][0], PETSC_USE_POINTER , &_isA[0]);
    CHKERRABORT(MPI_COMM_WORLD, ierr);
    //END Generate std::vector<IS> for GMRES solve by elimination ***********


    clock_t end_time = clock();
    SearchTime = (end_time - start_time);

    return SearchTime;
  }

// ================================================

  void GmresPetscLinearEquationSolver::solve(const vector <unsigned>& variable_to_be_solved, const bool& ksp_clean) {

    clock_t SearchTime, AssemblyTime, SolveTime, UpdateTime;
    PetscErrorCode ierr;

    // ***************** NODE/ELEMENT SEARCH *******************
    clock_t start_time = clock();
    if (_indexai_init == 0) BuildIndex(variable_to_be_solved);
    SearchTime = clock() - start_time;
    // ***************** END NODE/ELEMENT SEARCH *******************

    if (_DirichletBCsHandlingMode == 0) {   // By penalty
      PetscVector* EPSCp = static_cast<PetscVector*>(_EPSC);
      Vec EPSC = EPSCp->vec();
      PetscVector* RESp = static_cast<PetscVector*>(_RES);
      Vec RES = RESp->vec();
      PetscMatrix* KKp = static_cast<PetscMatrix*>(_KK);
      Mat KK = KKp->mat();

      // ***************** ASSEMBLE matrix to set Dirichlet BCs by penalty *******************
      start_time = clock();

      if (ksp_clean) {
        this->clear();
        // initialize Pmat wiwth penaly diagonal on the Dirichlet Nodes
        MatDuplicate(KK, MAT_COPY_VALUES, &_Pmat);
        MatSetOption(_Pmat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);
        MatZeroRows(_Pmat, _indexai[0].size(), &_indexai[0][0], 1.e100, 0, 0);
        _Pmat_is_initialized = true;
        this->init(KK, _Pmat);
      }
      AssemblyTime = clock() - start_time;
      // ***************** END ASSEMBLE ******************

      // ***************** SOLVE ******************
      start_time = clock();

      // Solve the linear system
      ierr = KSPSolve(_ksp, RES, EPSC);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      SolveTime = clock() - start_time;
      // ***************** END SOLVE ******************

      // ***************** RES/EPS UPDATE RES ******************
      start_time = clock();

      *_EPS += *_EPSC;

      _RESC->matrix_mult(*_EPSC, *_KK);
      *_RES -= *_RESC;

      UpdateTime = clock() - start_time;
      // ***************** END RES/EPS UPDATE ******************
    } else if (_DirichletBCsHandlingMode == 1) {   // By elimination

      PetscVector* RESp = static_cast<PetscVector*>(_RES);
      Vec RES = RESp->vec();
      PetscMatrix* KKp = static_cast<PetscMatrix*>(_KK);
      Mat KK = KKp->mat();
      PetscVector* EPSCp = static_cast<PetscVector*>(_EPSC);
      Vec EPSC = EPSCp->vec();
      PetscVector* EPSp = static_cast<PetscVector*>(_EPS);
      Vec EPS = EPSp->vec();

      // ***************** ASSEMBLE *******************
      clock_t start_time = clock();

      IS& isA = _isA[0];

      Vec Pr;
      ierr = VecGetSubVector(RES, isA, &Pr);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      // initialize _Pmat,_Ksp,_pc,_Pw,
      if (ksp_clean) {
        this->clear();
        ierr = MatGetSubMatrix(KK, isA, isA, MAT_INITIAL_MATRIX, &_Pmat);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        _Pmat_is_initialized = true;
        this->init(_Pmat, _Pmat);
        ierr = VecDuplicate(Pr, &_Pw);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        _Pw_is_initialized = true;
        ierr = VecScatterCreate(RES, isA, _Pw, NULL, &_scat);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        _scat_is_initialized = true;
      }

      AssemblyTime = clock() - start_time;
      // ***************** END ASSEMBLE ******************

      // ***************** SOLVE ******************
      start_time = clock();

      // Solve the linear system
      ierr = KSPSolve(_ksp, Pr, _Pw);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      ierr = VecRestoreSubVector(RES, isA, &Pr);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
      ierr = VecDestroy(&Pr);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      SolveTime = clock() - start_time;
      // ***************** END SOLVE ******************

      // ***************** RES/EPS UPDATE ******************
      start_time = clock();

      _EPSC->zero();

      ierr = VecScatterBegin(_scat, _Pw, EPSC, INSERT_VALUES, SCATTER_REVERSE);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
      ierr = VecScatterEnd(_scat, _Pw, EPSC, INSERT_VALUES, SCATTER_REVERSE);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      ierr = VecScatterBegin(_scat, _Pw, EPS, ADD_VALUES, SCATTER_REVERSE);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
      ierr = VecScatterEnd(_scat, _Pw, EPS, ADD_VALUES, SCATTER_REVERSE);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      _RESC->matrix_mult(*_EPSC, *_KK);
      *_RES -= *_RESC;

      UpdateTime = clock() - start_time;
      // ***************** END RES/EPS UPDATE ******************
    }

    // *** Computational info ***
#ifndef NDEBUG
    cout << "GMRES Grid: " << _msh->GetLevel() << "      SOLVER TIME:        "  << std::setw(11) << std::setprecision(6) << std::fixed <<
         static_cast<double>(SearchTime + AssemblyTime + SolveTime + UpdateTime) / CLOCKS_PER_SEC <<
         "  ITS: " << _maxits  << "\t ksp_clean = " << ksp_clean << endl;
#endif

  }



  void GmresPetscLinearEquationSolver::MGsetLevels(
    LinearEquationSolver* LinSolver, const unsigned& level, const unsigned& levelMax,
    const vector <unsigned>& variable_to_be_solved, SparseMatrix* PP, SparseMatrix* RR,
    const unsigned& npre, const unsigned& npost) {

    PetscErrorCode ierr;

    // ***************** NODE/ELEMENT SEARCH *******************
    if (_indexai_init == 0) BuildIndex(variable_to_be_solved);
    // ***************** END NODE/ELEMENT SEARCH *******************

    if (_DirichletBCsHandlingMode != 0) {
      std::cout << "Warning MGsolve does not allow BC by ELIMINATION, switched to PENALTY" << std::endl;
      _DirichletBCsHandlingMode = 0;
    }

    KSP* kspMG = LinSolver->GetKSP();
    PC pcMG;
    KSPGetPC(*kspMG, &pcMG);

    KSP subksp;
    KSP subkspUp;
    if (level == 0){
      PCMGGetCoarseSolve(pcMG, &subksp);
    }
    else {
      PCMGGetSmoother(pcMG, level , &subksp);
      KSPSetTolerances(subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npre);
      if (npre != npost) {
        PCMGGetSmootherUp(pcMG, level , &subkspUp);
        KSPSetTolerances(subkspUp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npost);
        this->set_petsc_solver_type(subkspUp);
      }
    }
    this->set_petsc_solver_type(subksp);

    if (_Pmat_is_initialized) MatDestroy(&_Pmat);

    PetscMatrix* KKp = static_cast<PetscMatrix*>(_KK);
    Mat KK = KKp->mat();

    MatDuplicate(KK, MAT_COPY_VALUES, &_Pmat);
    MatSetOption(_Pmat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);
    MatZeroRows(_Pmat, _indexai[0].size(), &_indexai[0][0], 1.e100, 0, 0);
    _Pmat_is_initialized = true;

    std::ostringstream levelName;
    levelName << "level-" << level;

    KSPSetOptionsPrefix(subksp, levelName.str().c_str());
    KSPSetFromOptions(subksp);
    KSPSetOperators(subksp, KK, _Pmat);

    PC subpc;

    ierr = KSPGetPC(subksp, &subpc);

    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpc);
    PetscReal zero = 1.e-16;
    PCFactorSetZeroPivot(subpc, zero);
    PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);

    if (level < levelMax) {
      PetscVector* EPSp = static_cast< PetscVector* >(_EPS);
      Vec EPS = EPSp->vec();
      PetscVector* RESp = static_cast< PetscVector* >(_RES);
      Vec RES = RESp->vec();
      PCMGSetX(pcMG, level, EPS);
      PCMGSetRhs(pcMG, level, RES);
    }
    if (level > 0) {
      PetscVector* RESCp = static_cast<PetscVector*>(_RESC);
      Vec RESC = RESCp->vec();
      PCMGSetR(pcMG, level, RESC);

      PetscMatrix* PPp = static_cast< PetscMatrix* >(PP);
      Mat P = PPp->mat();
      PCMGSetInterpolation(pcMG, level, P);

      PetscMatrix* RRp = static_cast< PetscMatrix* >(RR);
      Mat R = RRp->mat();
      PCMGSetRestriction(pcMG, level, R);

      if (npre != npost) {
        KSPSetPC(subkspUp, subpc);
      }

    }


  }

  void GmresPetscLinearEquationSolver::MGsolve(const bool ksp_clean) {

    if (ksp_clean) {
      PetscMatrix* KKp = static_cast< PetscMatrix* >(_KK);
      Mat KK = KKp->mat();
      KSPSetOperators(_ksp, KK, _Pmat);

      KSPSetTolerances(_ksp, _rtol, _abstol, _dtol, _maxits);
      KSPSetInitialGuessKnoll(_ksp, PETSC_TRUE);

      KSPSetFromOptions(_ksp);
    }

    PetscVector* EPSCp = static_cast< PetscVector* >(_EPSC);
    Vec EPSC = EPSCp->vec();
    PetscVector* RESp = static_cast< PetscVector* >(_RES);
    Vec RES = RESp->vec();

    KSPSolve(_ksp, RES, EPSC);

    _RESC->matrix_mult(*_EPSC, *_KK);
    *_RES -= *_RESC;

    *_EPS += *_EPSC;

#ifndef NDEBUG
    int its;
    KSPGetIterationNumber(_ksp, &its);

    KSPConvergedReason reason;
    KSPGetConvergedReason(_ksp, &reason);

    std::cout << "Number of iterations = " << its << "\t convergence reason = " << reason << std::endl;
    std::cout << _rtol << " " << _abstol << " " << _dtol << " " << _maxits << std::endl;
#endif

  }

// ================================================

  void GmresPetscLinearEquationSolver::clear() {

    int ierr;
    if (_Pmat_is_initialized) {
      _Pmat_is_initialized = false;
      ierr = MatDestroy(&_Pmat);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }
    if (_scat_is_initialized) {
      _scat_is_initialized = false;
      ierr = VecScatterDestroy(&_scat);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }
    if (_Pw_is_initialized) {
      _Pw_is_initialized = false;
      ierr = VecDestroy(&_Pw);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }

    if (this->initialized()) {
      this->_is_initialized = false;
      int ierr = 0;
      ierr = KSPDestroy(&_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }
  }

// ================================================

  void GmresPetscLinearEquationSolver::init(Mat& Amat, Mat& Pmat) {

    // Initialize the data structures if not done so already.
    if (!this->initialized())    {
      this->_is_initialized = true;
      int ierr = 0;
      // Create the linear solver context
      ierr = KSPCreate(MPI_COMM_WORLD, &_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      // Create the preconditioner context
      ierr = KSPGetPC(_ksp, &_pc);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type(_ksp);


//       ierr = KSPSetOperators(_ksp, Amat, Pmat, SAME_PRECONDITIONER);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = KSPSetOperators(_ksp, Amat, Pmat);
      CHKERRABORT(MPI_COMM_WORLD, ierr);    //PETSC3p5

      // Set the tolerances for the iterative solver.  Use the user-supplied
      // tolerance for the relative residual & leave the others at default values.
      ierr = KSPSetTolerances(_ksp, _rtol, _abstol, _dtol, _maxits);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      if (_msh->GetLevel() != 0)
        KSPSetInitialGuessKnoll(_ksp, PETSC_TRUE);

      if (_msh->GetLevel() != 0)
        KSPSetNormType(_ksp, KSP_NORM_NONE);

      // Set the options from user-input
      // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization  routines.
      ierr = KSPSetFromOptions(_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      // Notify PETSc of location to store residual history.
      // This needs to be called before any solves, since
      // it sets the residual history length to zero.  The default
      // behavior is for PETSc to allocate (internally) an array
      // of size 1000 to hold the residual norm history.
//       ierr = KSPSetResidualHistory(_ksp,
// 				   PETSC_NULL,   // pointer to the array which holds the history
// 				   PETSC_DECIDE, // size of the array holding the history
// 				   PETSC_TRUE);  // Whether or not to reset the history for each solve.
//       CHKERRABORT(MPI_COMM_WORLD,ierr);

      //PCSetType(_pc,PCREDISTRIBUTE);
      PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, _pc);
      PetscReal zero = 1.e-16;
      PCFactorSetZeroPivot(_pc, zero);
      PCFactorSetShiftType(_pc, MAT_SHIFT_NONZERO);
    }
  }

// ================================================

  void GmresPetscLinearEquationSolver::set_petsc_solver_type(KSP& ksp) {
    int ierr = 0;
    switch (this->_solver_type) {
      case CG:
        ierr = KSPSetType(ksp, (char*) KSPCG);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case CR:
        ierr = KSPSetType(ksp, (char*) KSPCR);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case CGS:
        ierr = KSPSetType(ksp, (char*) KSPCGS);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case BICG:
        ierr = KSPSetType(ksp, (char*) KSPBICG);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case TCQMR:
        ierr = KSPSetType(ksp, (char*) KSPTCQMR);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case TFQMR:
        ierr = KSPSetType(ksp, (char*) KSPTFQMR);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case LSQR:
        ierr = KSPSetType(ksp, (char*) KSPLSQR);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case BICGSTAB:
        ierr = KSPSetType(ksp, (char*) KSPBCGS);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case MINRES:
        ierr = KSPSetType(ksp, (char*) KSPMINRES);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case GMRES:
        ierr = KSPSetType(ksp, (char*) KSPGMRES);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case RICHARDSON:
        ierr = KSPSetType(ksp, (char*) KSPRICHARDSON);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case CHEBYSHEV:
        ierr = KSPSetType(ksp, (char*) KSPCHEBYSHEV);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      case PREONLY:
        ierr = KSPSetType(ksp, (char*) KSPPREONLY);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        return;
      default:
        std::cerr << "ERROR:  Unsupported PETSC Solver: "
                  << this->_solver_type               << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }
  }



// ========================================================
  std::pair<unsigned int, double> GmresPetscLinearEquationSolver::solve(SparseMatrix&  matrix_in,
      SparseMatrix&  precond_in,  NumericVector& solution_in,  NumericVector& rhs_in,
      const double tol,   const unsigned int m_its) {

//   START_LOG("solve()", "PetscLinearSolverM");
    // Make sure the data passed in are really of Petsc types
    PetscMatrix* matrix   = libmeshM_cast_ptr<PetscMatrix*> (&matrix_in);
    PetscMatrix* precond  = libmeshM_cast_ptr<PetscMatrix*> (&precond_in);
    PetscVector* solution = libmeshM_cast_ptr<PetscVector*> (&solution_in);
    PetscVector* rhs      = libmeshM_cast_ptr<PetscVector*> (&rhs_in);
    this->init(matrix);

    int ierr = 0;
    int its = 0, max_its = static_cast<int>(m_its);
    PetscReal final_resid = 0.;
    // Close the matrices and vectors in case this wasn't already done.
    matrix->close();
    precond->close();
    solution->close();
    rhs->close();
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
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat());    //PETSC3p5
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    } else  {
      //ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),SAME_PRECONDITIONER);
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat());    //PETSC3p5
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }
    // Set the tolerances for the iterative solver.  Use the user-supplied
    // tolerance for the relative residual & leave the others at default values.
    ierr = KSPSetTolerances(_ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_its);
    CHKERRABORT(MPI_COMM_WORLD, ierr);
    // Solve the linear system

//        PetscLogEvent USER_EVENT;
//      PetscLogDouble user_event_flops;
//      PetscLogEventRegister("User event",0,&USER_EVENT);
//      PetscLogEventBegin(USER_EVENT,0,0,0,0);


    ierr = KSPSolve(_ksp, rhs->vec(), solution->vec());
    CHKERRABORT(MPI_COMM_WORLD, ierr);
//         PetscLogFlops(user_event_flops);
//      PetscLogEventEnd(USER_EVENT,0,0,0,0);

    // Get the number of iterations required for convergence
    ierr = KSPGetIterationNumber(_ksp, &its);
    CHKERRABORT(MPI_COMM_WORLD, ierr);
    // Get the norm of the final residual to return to the user.
    ierr = KSPGetResidualNorm(_ksp, &final_resid);
    CHKERRABORT(MPI_COMM_WORLD, ierr);

//   STOP_LOG("solve()", "PetscLinearSolverM");
    return std::make_pair(its, final_resid);
  }


// @deprecated ========================================================
  PetscErrorCode __libmesh_petsc_preconditioner_setup(PC pc) {
    void* ctx;
    PetscErrorCode ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    Preconditioner* preconditioner = static_cast<Preconditioner*>(ctx);
    preconditioner->init();
    return 0;
  }


// @deprecated ========================================================
  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y) {
    void* ctx;
    PetscErrorCode ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    Preconditioner* preconditioner = static_cast<Preconditioner*>(ctx);
    PetscVector x_vec(x);
    PetscVector y_vec(y);
    preconditioner->apply(x_vec, y_vec);
    return 0;
  }


// @deprecated ========================================================
  void GmresPetscLinearEquationSolver::init(SparseMatrix* matrix) {

    PetscMatrix* matrix_two   = libmeshM_cast_ptr<PetscMatrix*> (matrix);

    // Initialize the data structures if not done so already.
    if (!this->initialized())    {
      this->_is_initialized = true;
      int ierr = 0;
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
      ierr = KSPCreate(MPI_COMM_WORLD, &_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
      //ierr = PCCreate (MPI_COMM_WORLD, &_pc); CHKERRABORT(MPI_COMM_WORLD,ierr);
      // Create the preconditioner context
      ierr = KSPGetPC(_ksp, &_pc);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
      // Set operators. The input matrix works as the preconditioning matrix
      //ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat(),SAME_NONZERO_PATTERN);
      ierr = KSPSetOperators(_ksp, matrix_two->mat(), matrix_two->mat());    //PETSC3p5
      CHKERRABORT(MPI_COMM_WORLD, ierr);
      // Have the Krylov subspace method use our good initial guess rather than 0
      ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
      // Set user-specified  solver and preconditioner types
      this->set_petsc_solver_type(_ksp);
      // Set the options from user-input
      // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization  routines.
      ierr = KSPSetFromOptions(_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
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
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, _pc);
      if (this->_preconditioner) {
        this->_preconditioner->set_matrix(*matrix);
        PCShellSetContext(_pc, (void*) this->_preconditioner);
        PCShellSetSetUp(_pc, __libmesh_petsc_preconditioner_setup);
        PCShellSetApply(_pc, __libmesh_petsc_preconditioner_apply);
      }
    }
  }



} //end namespace femus


#endif

