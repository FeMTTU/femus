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
#include "FieldSplitPetscLinearEquationSolver.hpp"
#include "MeshASMPartitioning.hpp"
#include "PetscPreconditioner.hpp"
#include "PetscMatrix.hpp"
#include <iomanip>
#include <sstream>

namespace femus {



  using namespace std;

  // ====================================================
  // ------------------- Class functions ------------
  // ====================================================

  void FieldSplitPetscLinearEquationSolver::set_tolerances(const double& rtol, const double& atol,
      const double& divtol, const unsigned& maxits) {

    _rtol   = static_cast<PetscReal>(rtol);
    _abstol = static_cast<PetscReal>(atol);
    _dtol   = static_cast<PetscReal>(divtol);
    _maxits = static_cast<PetscInt>(maxits);

  }

  // ==============================================
  void FieldSplitPetscLinearEquationSolver::SetElementBlockNumber(const char all[], const unsigned& overlap) {
    _element_block_number[0] = _msh->GetNumberOfElements();
    _element_block_number[1] = _msh->GetNumberOfElements();
    _standard_ASM = 1;
    _overlap = overlap;
  }

  void FieldSplitPetscLinearEquationSolver::SetElementBlockNumber(const unsigned& block_elemet_number) {
    _element_block_number[0] = block_elemet_number;
    _element_block_number[1] = block_elemet_number;
    _indexai_init = 0;
    _standard_ASM = 0;
  }

  void FieldSplitPetscLinearEquationSolver::SetElementBlockNumberSolid(const unsigned& block_elemet_number, const unsigned& overlap) {
    _element_block_number[0] = block_elemet_number;
    _indexai_init = 0;
    _standard_ASM = 0;
    _overlap = overlap;
  }

  void FieldSplitPetscLinearEquationSolver::SetElementBlockNumberFluid(const unsigned& block_elemet_number, const unsigned& overlap) {
    _element_block_number[1] = block_elemet_number;
    _indexai_init = 0;
    _standard_ASM = 0;
    _overlap = overlap;
  }

  // ==============================================
  clock_t FieldSplitPetscLinearEquationSolver::BuildBDCIndex(const vector <unsigned>& variable_to_be_solved) {

    clock_t SearchTime = 0;
    clock_t start_time = clock();

    unsigned IndexaSize = KKoffset[KKIndex.size() - 1][processor_id()] - KKoffset[0][processor_id()];
    _indexai.resize(2);

    _indexai[0].clear();
    _indexai[1].clear();

    _indexai[0].resize(IndexaSize);
    _indexai[1].resize(IndexaSize);

    unsigned count0 = 0;
    unsigned count1 = 0;


    vector <bool> ThisSolutionIsIncluded(_SolPdeIndex.size(), false);
    for (unsigned iind = 0; iind < variable_to_be_solved.size(); iind++) {
      unsigned PdeIndexSol = variable_to_be_solved[iind];
      ThisSolutionIsIncluded[PdeIndexSol] = true;
    }

    for (int k = 0; k < _SolPdeIndex.size(); k++) {
      unsigned indexSol = _SolPdeIndex[k];
      unsigned soltype = _SolType[indexSol];
      for (unsigned inode_mts = _msh->MetisOffset[soltype][processor_id()]; inode_mts < _msh->MetisOffset[soltype][processor_id() + 1]; inode_mts++) {
        int local_mts = inode_mts - _msh->MetisOffset[soltype][processor_id()];
        int idof_kk = KKoffset[k][processor_id()] + local_mts;
        if (!ThisSolutionIsIncluded[k] || (*(*_Bdc)[indexSol])(inode_mts) < 1.9) {
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

    clock_t end_time = clock();
    SearchTime = (end_time - start_time);

    return SearchTime;
  }


  // ==============================================

  clock_t FieldSplitPetscLinearEquationSolver::BuildFieldSplitIndex(const vector <unsigned>& variable_to_be_solved) {
    clock_t SearchTime = 0;
    clock_t start_time = clock();

//     unsigned nVariables = variable_to_be_solved.size();
//     unsigned iproc = processor_id(); 
//     
//     _is_loc_idx.resize(nVariables);
//     _is_loc.resize(nVariables);
//     
//     for( unsigned i = 0; i < nVariables; i++){
//       unsigned offset = KKoffset[i][iproc];
//       unsigned offsetp1 = KKoffset[i+1][iproc];
//       unsigned variableSize = offsetp1 - offset;
//       _is_loc_idx[i].resize(variableSize);
//       for(int j = 0; j < variableSize; j++){
// 	_is_loc_idx[i][j] = offset + j;
//       }
//       PetscErrorCode ierr;
//       ierr = ISCreateGeneral(MPI_COMM_WORLD, _is_loc_idx[i].size(), &_is_loc_idx[i][0], PETSC_USE_POINTER, &_is_loc[i]);
//       CHKERRABORT(MPI_COMM_WORLD, ierr);
//     }
    
    
    unsigned nVariables = 2;
    unsigned iproc = processor_id(); 
    
    _is_loc_idx.resize(nVariables);
    _is_loc.resize(nVariables);
    
    for( unsigned i = 0; i < nVariables; i++){
      unsigned start = (i == 0) ? 0 : 3;
      unsigned end = (i == 1) ? 3:4;
      unsigned offset = KKoffset[start][iproc];
      unsigned offsetp1 = KKoffset[end][iproc];
      unsigned variableSize = offsetp1 - offset;
      _is_loc_idx[i].resize(variableSize);
      for(int j = 0; j < variableSize; j++){
	_is_loc_idx[i][j] = offset + j;
      }
      PetscErrorCode ierr;
      ierr = ISCreateGeneral(MPI_COMM_WORLD, _is_loc_idx[i].size(), &_is_loc_idx[i][0], PETSC_USE_POINTER, &_is_loc[i]);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }
    
    clock_t end_time = clock();
    SearchTime += (end_time - start_time);
    return SearchTime;
  }

  // =================================================

  void FieldSplitPetscLinearEquationSolver::solve(const vector <unsigned>& variable_to_be_solved, const bool& ksp_clean) {
    PetscVector* EPSCp = static_cast<PetscVector*>(_EPSC);
    Vec EPSC = EPSCp->vec();
    PetscVector* RESp = static_cast<PetscVector*>(_RES);
    Vec RES = RESp->vec();
    PetscMatrix* KKp = static_cast<PetscMatrix*>(_KK);
    Mat KK = KKp->mat();

    PetscErrorCode ierr;
    clock_t SearchTime, AssemblyTime, SolveTime, UpdateTime;

    // ***************** NODE/ELEMENT SEARCH *******************
    clock_t start_time = clock();
    if (_indexai_init == 0) {
      _indexai_init = 1;
      if (!_standard_ASM)
        BuildFieldSplitIndex(variable_to_be_solved);
      BuildBDCIndex(variable_to_be_solved);
    }
    SearchTime = start_time - clock();
    // ***************** END NODE/ELEMENT SEARCH *******************

    // ***************** ASSEMBLE matrix to set Dirichlet BCs by penalty *******************
    start_time = clock();
    if (ksp_clean) {
      this->clear();
      // initialize Pmat wiwth penaly diagonal on the Dirichlet Nodes
      MatDuplicate(KK, MAT_COPY_VALUES, &_Pmat);
      MatSetOption(_Pmat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);
      MatZeroRows(_Pmat, _indexai[0].size(), &_indexai[0][0], 1.e100, 0, 0);
      _Pmat_is_initialized = true;



//       PetscViewer    viewer;
//       ierr=PetscViewerDrawOpen(MPI_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);
//       ierr= MatView(_Pmat,viewer);
//       double ff;
//       std::cin>>ff;
//       PetscViewerDestroy(&viewer);
//
      init(KK, _Pmat);

    }

    AssemblyTime = clock() - start_time;

    // ***************** END ASSEMBLE ***********

    // ***************** SOLVE ******************
    start_time = clock();

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

    // **************** END RES/EPS UPDATE RES ***************

    // *** Computational info ***
//#ifndef NDEBUG
    int its;
    KSPGetIterationNumber(_ksp, &its);

    KSPConvergedReason reason;
    KSPGetConvergedReason(_ksp, &reason);

    PetscReal rnorm;
    KSPGetResidualNorm(_ksp, &rnorm);

    std::cout << "Number of iterations = " << its << "\t convergence reason = " << reason << std::endl;
    std::cout << "Residual Norm ="<< rnorm <<std::endl;
    std::cout << _rtol << " " << _abstol << " " << _dtol << " " << _maxits << std::endl;




//     cout << "ASM Grid: " << _msh->GetLevel() << "        SOLVER TIME:        "  << std::setw(11) << std::setprecision(6) << std::fixed <<
//          static_cast<double>(SearchTime + AssemblyTime + SolveTime + UpdateTime) / CLOCKS_PER_SEC <<
//          "  ITS: " << _maxits  << "\t ksp_clean = " << ksp_clean << endl;
//#endif

  }


  void FieldSplitPetscLinearEquationSolver::MGsetLevels(
    LinearEquationSolver* LinSolver, const unsigned& level, const unsigned& levelMax,
    const vector <unsigned>& variable_to_be_solved, SparseMatrix* PP, SparseMatrix* RR,
    const unsigned& npre, const unsigned& npost) {

    // ***************** NODE/ELEMENT SEARCH *******************
    if (_indexai_init == 0) {
      _indexai_init = 1;
      if (!_standard_ASM)
	//BEGIN here
        BuildFieldSplitIndex(variable_to_be_solved);
        //END here
      BuildBDCIndex(variable_to_be_solved);
    }

    // ***************** END NODE/ELEMENT SEARCH *******************
    KSP* kspMG = LinSolver->GetKSP();
    PC pcMG;
    KSPGetPC(*kspMG, &pcMG);

    PetscMatrix* KKp = static_cast< PetscMatrix* >(_KK);
    Mat KK = KKp->mat();

    if (_Pmat_is_initialized) MatDestroy(&_Pmat);
    MatDuplicate(KK, MAT_COPY_VALUES, &_Pmat);
    MatSetOption(_Pmat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);
    MatZeroRows(_Pmat, _indexai[0].size(), &_indexai[0][0], 1.e100, 0, 0);
    _Pmat_is_initialized = true;


    KSP subksp;
    KSP subkspUp;
    if (level == 0)
      PCMGGetCoarseSolve(pcMG, &subksp);
    else {
      PCMGGetSmoother(pcMG, level , &subksp);
      //KSPSetInitialGuessKnoll(subksp, PETSC_TRUE);
      KSPSetTolerances(subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npre);
      if (npre != npost) {
        PCMGGetSmootherUp(pcMG, level , &subkspUp);
        KSPSetTolerances(subkspUp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npost);
        this->set_petsc_solver_type(subkspUp);
      }
    }
    this->set_petsc_solver_type(subksp);

    std::ostringstream levelName;
    levelName << "level-" << level;

    KSPSetOptionsPrefix(subksp, levelName.str().c_str());
    KSPSetFromOptions(subksp);
    KSPSetOperators(subksp, KK, _Pmat);

    PC subpc;
    KSPGetPC(subksp, &subpc);

    //BEGIN from here
    
    PCSetType(subpc, (char*) PCFIELDSPLIT);
    PCFieldSplitSetType(subpc, PC_COMPOSITE_ADDITIVE);
    for(int i=0; i<_is_loc.size(); i++ ){
      PCFieldSplitSetIS( subpc, NULL, _is_loc[i]);
    }
    
    KSPSetUp(subksp);
       
    KSP* subksps;
    PCFieldSplitGetSubKSP(subpc, &_nlocal, &subksps);
// 
    PetscReal epsilon = 1.e-16;
    for (int i = 0; i < _is_loc.size(); i++) {
      KSPSetType(subksps[i], (char*) KSPPREONLY);
      PC subpcs;
      KSPGetPC(subksps[i], &subpcs);
      KSPSetTolerances(subksps[i], _rtol, _abstol, _dtol, 1);
      KSPSetFromOptions(subksps[i]);
      
      PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);
      PCFactorSetZeroPivot(subpcs, epsilon);
      PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
              
    }
      

    //END here
    
    if (level < levelMax) {   //all but finest
      PetscVector* EPSp = static_cast< PetscVector* >(_EPS);
      Vec EPS = EPSp->vec();
      PetscVector* RESp = static_cast< PetscVector* >(_RES);
      Vec RES = RESp->vec();
      PCMGSetX(pcMG, level, EPS);
      PCMGSetRhs(pcMG, level, RES);
    }
    if (level > 0) {   //all but coarsest
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


  void FieldSplitPetscLinearEquationSolver::MGsolve(const bool ksp_clean) {

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

    PetscReal rnorm;
    KSPGetResidualNorm(_ksp, &rnorm);

    std::cout << "Number of iterations = " << its << "\t convergence reason = " << reason << std::endl;
    std::cout << "Residual Norm ="<< rnorm <<std::endl;
    std::cout << _rtol << " " << _abstol << " " << _dtol << " " << _maxits << std::endl;
#endif

  }

// ================================================

  void FieldSplitPetscLinearEquationSolver::clear() {
    int ierr = 0;
    if (_Pmat_is_initialized) {
      _Pmat_is_initialized = false;
      ierr = MatDestroy(&_Pmat);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }

    if (this->initialized()) {
      this->_is_initialized = false;
      ierr = KSPDestroy(&_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }
  }

// ========================================================

  void FieldSplitPetscLinearEquationSolver::init(Mat& Amat, Mat& Pmat) {

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
      // Set operators. The input matrix works as the preconditioning matrix
      this->set_petsc_solver_type(_ksp);

      //ierr = KSPSetOperators(_ksp, Amat, Pmat, SAME_PRECONDITIONER);		CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = KSPSetOperators(_ksp, Amat, Pmat);
      CHKERRABORT(MPI_COMM_WORLD, ierr); //PETSC3p5

      // Set the tolerances for the iterative solver.  Use the user-supplied
      // tolerance for the relative residual & leave the others at default values.
      ierr = KSPSetTolerances(_ksp, _rtol, _abstol, _dtol, _maxits);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      if ( _msh->GetLevel() != 0)
        KSPSetInitialGuessKnoll(_ksp, PETSC_TRUE);

      if (_msh->GetLevel() != 0)
        KSPSetNormType(_ksp, KSP_NORM_NONE);

      ierr = KSPSetFromOptions(_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      PetscPreconditioner::set_petsc_preconditioner_type(ASM_PRECOND, _pc);
      if (!_standard_ASM) {
        ierr = PCASMSetLocalSubdomains(_pc, _is_loc_idx.size(), &_is_ovl[0], &_is_loc[0]);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
      }
      ierr = PCASMSetOverlap(_pc, _overlap);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
      //ierr = PCASMSetLocalType(_pc, PC_COMPOSITE_MULTIPLICATIVE);                   CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = KSPSetUp(_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);


      ierr = PCASMGetSubKSP(_pc, &_nlocal, &_first, &_ksp_asm);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      PetscReal epsilon = 1.e-16;

      if (!_standard_ASM) {
        _pc_asm.resize(2);
        for (int i = 0; i < _block_type_range[0]; i++) {
          ierr = KSPGetPC(_ksp_asm[i], &_pc_asm[0]);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = KSPSetTolerances(_ksp_asm[i], _rtol, _abstol, _dtol, 1);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = KSPSetFromOptions(_ksp_asm[i]);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          PetscPreconditioner::set_petsc_preconditioner_type(MLU_PRECOND, _pc_asm[0]);
          ierr = PCFactorSetZeroPivot(_pc_asm[0], epsilon);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = PCFactorSetShiftType(_pc_asm[0], MAT_SHIFT_NONZERO);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
        }
        for (int i = _block_type_range[0]; i < _block_type_range[1]; i++) {
          ierr = KSPGetPC(_ksp_asm[i], &_pc_asm[1]);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = KSPSetTolerances(_ksp_asm[i], _rtol, _abstol, _dtol, 1);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = KSPSetFromOptions(_ksp_asm[i]);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, _pc_asm[1]);
          ierr = PCFactorSetZeroPivot(_pc_asm[1], epsilon);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = PCFactorSetShiftType(_pc_asm[1], MAT_SHIFT_NONZERO);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
        }
      } else {
        _pc_asm.resize(1);
        for (int i = 0; i < _nlocal; i++) {
          ierr = KSPGetPC(_ksp_asm[i], &_pc_asm[0]);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = KSPSetTolerances(_ksp_asm[i], _rtol, _abstol, _dtol, 1);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = KSPSetFromOptions(_ksp_asm[i]);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, _pc_asm[0]);
          ierr = PCFactorSetZeroPivot(_pc_asm[0], epsilon);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = PCFactorSetShiftType(_pc_asm[0], MAT_SHIFT_NONZERO);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
        }
      }
    }
  }

// =================================================

  void FieldSplitPetscLinearEquationSolver::set_petsc_solver_type(KSP& ksp) {
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
} //end namespace femus

#endif

