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
#include "AsmPetscLinearEquationSolver.hpp"
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

  void AsmPetscLinearEquationSolver::set_tolerances(const double& rtol, const double& atol,
      const double& divtol, const unsigned& maxits, const unsigned& restart) {

    _rtol    = static_cast<PetscReal>(rtol);
    _abstol  = static_cast<PetscReal>(atol);
    _dtol    = static_cast<PetscReal>(divtol);
    _maxits  = static_cast<PetscInt>(maxits);
    _restart = static_cast<PetscInt>(restart);

  }

  // ==============================================
  void AsmPetscLinearEquationSolver::SetElementBlockNumber(const char all[], const unsigned& overlap) {
    _element_block_number[0] = _msh->GetNumberOfElements();
    _element_block_number[1] = _msh->GetNumberOfElements();
    _standard_ASM = 1;
    _overlap = overlap;
  }

  void AsmPetscLinearEquationSolver::SetElementBlockNumber(const unsigned& block_elemet_number) {
    _element_block_number[0] = block_elemet_number;
    _element_block_number[1] = block_elemet_number;
    _indexai_init = 0;
    _standard_ASM = 0;
  }

  void AsmPetscLinearEquationSolver::SetElementBlockNumberSolid(const unsigned& block_elemet_number, const unsigned& overlap) {
    _element_block_number[0] = block_elemet_number;
    _indexai_init = 0;
    _standard_ASM = 0;
    _overlap = overlap;
  }

  void AsmPetscLinearEquationSolver::SetElementBlockNumberFluid(const unsigned& block_elemet_number, const unsigned& overlap) {
    _element_block_number[1] = block_elemet_number;
    _indexai_init = 0;
    _standard_ASM = 0;
    _overlap = overlap;
  }

  // ==============================================
  clock_t AsmPetscLinearEquationSolver::BuildBDCIndex(const vector <unsigned>& variable_to_be_solved) {

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
      for (unsigned inode_mts = _msh->_dofOffset[soltype][processor_id()]; inode_mts < _msh->_dofOffset[soltype][processor_id() + 1]; inode_mts++) {
        int local_mts = inode_mts - _msh->_dofOffset[soltype][processor_id()];
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

  clock_t AsmPetscLinearEquationSolver::BuildAMSIndex(const vector <unsigned>& variable_to_be_solved) {
    clock_t SearchTime = 0;
    clock_t start_time = clock();

    unsigned nel = _msh->GetNumberOfElements();

    bool FastVankaBlock = true;
    if (_NSchurVar == !0) {
      FastVankaBlock = (_SolType[_SolPdeIndex[variable_to_be_solved[variable_to_be_solved.size() - _NSchurVar]]] < 3) ? false : true;
    }

    unsigned iproc = processor_id();

    unsigned DofOffset = KKoffset[0][iproc];
    unsigned DofOffsetSize = KKoffset[KKIndex.size() - 1][iproc] - KKoffset[0][iproc];
    vector < unsigned > indexa(DofOffsetSize, DofOffsetSize);
    vector < unsigned > indexb(DofOffsetSize, DofOffsetSize);
    vector <bool> owned(DofOffsetSize, false);

    map<int, bool> mymap;

    unsigned ElemOffset   = _msh->_dofOffset[3][iproc];
    unsigned ElemOffsetp1 = _msh->_dofOffset[3][iproc + 1];
    unsigned ElemOffsetSize = ElemOffsetp1 - ElemOffset;
    vector <PetscInt> indexci(ElemOffsetSize);
    vector < unsigned > indexc(ElemOffsetSize, ElemOffsetSize);

    vector < vector < unsigned > > block_elements;

    MeshASMPartitioning meshasmpartitioning(*_msh);
    meshasmpartitioning.DoPartition(_element_block_number, block_elements, _block_type_range);

    vector <bool> ThisVaribaleIsNonSchur(_SolPdeIndex.size(), true);
    for (unsigned iind = variable_to_be_solved.size() - _NSchurVar; iind < variable_to_be_solved.size(); iind++) {
      unsigned PdeIndexSol = variable_to_be_solved[iind];
      ThisVaribaleIsNonSchur[PdeIndexSol] = false;
    }

    // *** Start Vanka Block ***

    _is_loc_idx.resize(block_elements.size());
    _is_ovl_idx.resize(block_elements.size());

    for (int vb_index = 0; vb_index < block_elements.size(); vb_index++) {
      _is_loc_idx[vb_index].resize(DofOffsetSize);
      _is_ovl_idx[vb_index].resize(DofOffsetSize);

      PetscInt PAsize = 0;
      PetscInt PBsize = 0;

      PetscInt Csize = 0;

      // ***************** NODE/ELEMENT SERCH *******************
      for (int kel = 0; kel < block_elements[vb_index].size(); kel++) {
        unsigned iel = block_elements[vb_index][kel];

        for (unsigned i = 0; i < _msh->GetElementDofNumber(iel, 0); i++) {
          unsigned inode = _msh->el->GetElementDofIndex(iel, i);
          const std::vector < unsigned > & localElementNearVertexNumber = _msh->el->GetLocalElementNearVertex(inode);
          unsigned nve = (FastVankaBlock) ? 1 : localElementNearVertexNumber.size();
          for (unsigned j = 0; j < nve; j++) {
            unsigned jel = (!FastVankaBlock) ? localElementNearVertexNumber[j] : iel;
            //add elements for velocity to be solved
            if (indexc[jel - ElemOffset] == ElemOffsetSize) {
              indexci[Csize] = jel - ElemOffset;
              indexc[jel - ElemOffset] = Csize++;
              //add non-schur node to be solved
              for (int indexSol = 0; indexSol < _SolPdeIndex.size(); indexSol++) {
                if (ThisVaribaleIsNonSchur[indexSol]) {
                  unsigned SolPdeIndex = _SolPdeIndex[indexSol];
                  unsigned SolType = _SolType[SolPdeIndex];
                  unsigned nvej = _msh->GetElementDofNumber(jel, SolType);
                  for (unsigned jj = 0; jj < nvej; jj++) {
		    unsigned jdof = _msh->GetSolutionDof(jj, jel, SolType);
		    unsigned kkdof = GetSystemDof(SolPdeIndex, indexSol, jj, jel);
                    if (jdof >= _msh->_dofOffset[SolType][iproc] &&
                        jdof <  _msh->_dofOffset[SolType][iproc + 1]) {
                      if (indexa[kkdof - DofOffset] == DofOffsetSize && owned[kkdof - DofOffset] == false) {
                        owned[kkdof - DofOffset] = true;
                        _is_loc_idx[vb_index][PAsize] = kkdof;
                        indexa[kkdof - DofOffset] = PAsize++;
                      }
                      if (indexb[kkdof - DofOffset] == DofOffsetSize) {
                        _is_ovl_idx[vb_index][PBsize] = kkdof;
                        indexb[kkdof - DofOffset] = PBsize++;
                      }
                    } else {
                      mymap[kkdof] = true;
                    }
                  }
                }
              }
            }
          }
        }
        //-----------------------------------------------------------------------------------------
        //Add Schur nodes (generally pressure) to be solved
        {
          for (int indexSol = 0; indexSol < _SolPdeIndex.size(); indexSol++) {
            if (!ThisVaribaleIsNonSchur[indexSol]) {
              unsigned SolPdeIndex = _SolPdeIndex[indexSol];
              unsigned SolType = _SolType[SolPdeIndex];
              unsigned nvei = _msh->GetElementDofNumber(iel, SolType);
              for (unsigned ii = 0; ii < nvei; ii++) {
		unsigned inode_Metis = _msh->GetSolutionDof(ii, iel, SolType);
		unsigned kkdof = GetSystemDof(SolPdeIndex, indexSol, ii, iel);
                if (inode_Metis >= _msh->_dofOffset[SolType][iproc] &&
                    inode_Metis <  _msh->_dofOffset[SolType][iproc + 1]) {
                  if (indexa[kkdof - DofOffset] == DofOffsetSize && owned[kkdof - DofOffset] == false) {
                    owned[kkdof - DofOffset] = true;
                    _is_loc_idx[vb_index][PAsize] = kkdof;
                    indexa[kkdof - DofOffset] = PAsize++;
                  }
                  if (indexb[kkdof - DofOffset] == DofOffsetSize) {
                    _is_ovl_idx[vb_index][PBsize] = kkdof;
                    indexb[kkdof - DofOffset] = PBsize++;
                  }
                } else {
                  mymap[kkdof] = true;
                }
              }
            }
          }
        }
        //-----------------------------------------------------------------------------------------
      }

      // *** re-initialize indeces(a,c,d)
      for (PetscInt i = 0; i < PAsize; i++) {
        indexa[_is_loc_idx[vb_index][i] - DofOffset] = DofOffsetSize;
      }
      for (PetscInt i = 0; i < PBsize; i++) {
        indexb[_is_ovl_idx[vb_index][i] - DofOffset] = DofOffsetSize;
      }
      for (PetscInt i = 0; i < Csize; i++) {
        indexc[indexci[i]] = ElemOffsetSize;
      }

      _is_loc_idx[vb_index].resize(PAsize);

      _is_ovl_idx[vb_index].resize(PBsize + mymap.size());
      int i = 0;
      for (std::map<int, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it, ++i) {
        _is_ovl_idx[vb_index][PBsize + i] = it->first;
      }

      std::sort(_is_loc_idx[vb_index].begin(), _is_loc_idx[vb_index].end());
      std::sort(_is_ovl_idx[vb_index].begin(), _is_ovl_idx[vb_index].end());
    }

    //BEGIN Generate std::vector<IS> for vanka solve ***********
    _is_loc.resize(_is_loc_idx.size());
    _is_ovl.resize(_is_ovl_idx.size());

    for (unsigned vb_index = 0; vb_index < _is_loc_idx.size(); vb_index++) {
      PetscErrorCode ierr;
      ierr = ISCreateGeneral(MPI_COMM_SELF, _is_loc_idx[vb_index].size(), &_is_loc_idx[vb_index][0], PETSC_USE_POINTER, &_is_loc[vb_index]);
      CHKERRABORT(MPI_COMM_SELF, ierr);
      ierr = ISCreateGeneral(MPI_COMM_SELF, _is_ovl_idx[vb_index].size(), &_is_ovl_idx[vb_index][0], PETSC_USE_POINTER, &_is_ovl[vb_index]);
      CHKERRABORT(MPI_COMM_SELF, ierr);
    }
    //END Generate std::vector<IS> for vanka solve ***********

    clock_t end_time = clock();
    SearchTime += (end_time - start_time);
    return SearchTime;
  }

  // =================================================

  void AsmPetscLinearEquationSolver::solve(const vector <unsigned>& variable_to_be_solved, const bool& ksp_clean) {

    PetscLogDouble t1;
    PetscLogDouble t2;
    PetscTime(&t1);

    PetscVector* EPSCp = static_cast<PetscVector*>(_EPSC);
    Vec EPSC = EPSCp->vec();
    PetscVector* RESp = static_cast<PetscVector*>(_RES);
    Vec RES = RESp->vec();
    PetscMatrix* KKp = static_cast<PetscMatrix*>(_KK);
    Mat KK = KKp->mat();

    //BEGIN SEARCH
    if (_indexai_init == 0) {
      _indexai_init = 1;
      if (!_standard_ASM)
        BuildAMSIndex(variable_to_be_solved);
      BuildBDCIndex(variable_to_be_solved);
    }
    //END SEARCH

    //BEGIN ASSEMBLE matrix with Dirichlet penalty BCs by penalty
    if (ksp_clean) {
      this->clear();
      // initialize Pmat wiwth penaly diagonal on the Dirichlet Nodes
      MatDuplicate(KK, MAT_COPY_VALUES, &_Pmat);
      MatSetOption(_Pmat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);
      MatZeroRows(_Pmat, _indexai[0].size(), &_indexai[0][0], 1.e100, 0, 0);
      _Pmat_is_initialized = true;
      init(KK, _Pmat);
    }
    //END ASSEMBLE


    //BEGIN SOLVE and UPDATE
    KSPSolve(_ksp, RES, EPSC);
    *_EPS += *_EPSC;
    _RESC->matrix_mult(*_EPSC, *_KK);
    *_RES -= *_RESC;
    //END SOLVE and UPDATE


    //BEGIN PRINT Computational info
// #ifndef NDEBUG
    int its;
    KSPGetIterationNumber(_ksp, &its);

    KSPConvergedReason reason;
    KSPGetConvergedReason(_ksp, &reason);

    PetscReal rnorm;
    KSPGetResidualNorm(_ksp, &rnorm);

    PetscTime(&t2);
    PetscPrintf(PETSC_COMM_WORLD, " *************** MG linear solver time: %8.3f \n", t2 -t1);
    PetscPrintf(PETSC_COMM_WORLD, " *************** Number of outer ksp solver iterations = %i \n", its);
    PetscPrintf(PETSC_COMM_WORLD, " *************** Convergence reason = %i \n", reason);
    PetscPrintf(PETSC_COMM_WORLD, " *************** Residual norm = %10.8g \n", rnorm);

// #endif
    //END PRINT
  }


  void AsmPetscLinearEquationSolver::MGsetLevels(
    LinearEquationSolver* LinSolver, const unsigned& level, const unsigned& levelMax,
    const vector <unsigned>& variable_to_be_solved, SparseMatrix* PP, SparseMatrix* RR,
    const unsigned& npre, const unsigned& npost) {

    // ***************** NODE/ELEMENT SEARCH *******************
    if (_indexai_init == 0) {
      _indexai_init = 1;
      if (!_standard_ASM)
        BuildAMSIndex(variable_to_be_solved);
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
      KSPSetTolerances(subksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npre);
    }
    this->set_petsc_solver_type(subksp);

    std::ostringstream levelName;
    levelName << "level-" << level;

    KSPSetOptionsPrefix(subksp, levelName.str().c_str());
    KSPSetFromOptions(subksp);
    KSPSetOperators(subksp, KK, _Pmat);

    PC subpc;
    KSPGetPC(subksp, &subpc);

    PetscPreconditioner::set_petsc_preconditioner_type(ASM_PRECOND, subpc);
    if (!_standard_ASM) {
      PCASMSetLocalSubdomains(subpc, _is_loc_idx.size(), &_is_ovl[0], &_is_loc[0]);
    }
    PCASMSetOverlap(subpc, _overlap);
    //PCASMSetLocalType(subpc, PC_COMPOSITE_MULTIPLICATIVE);

    KSPSetUp(subksp);
    KSP* subksps;
    PCASMGetSubKSP(subpc, &_nlocal, PETSC_NULL, &subksps);

    PetscReal epsilon = 1.e-16;
    if (!_standard_ASM) {
      for (int i = 0; i < _block_type_range[0]; i++) {
        PC subpcs;
        KSPGetPC(subksps[i], &subpcs);
        KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
        KSPSetFromOptions(subksps[i]);
        PetscPreconditioner::set_petsc_preconditioner_type(MLU_PRECOND, subpcs);
        PCFactorSetZeroPivot(subpcs, epsilon);
        PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
      }
      for (int i = _block_type_range[0]; i < _block_type_range[1]; i++) {
        PC subpcs;
        KSPGetPC(subksps[i], &subpcs);
        KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
        KSPSetFromOptions(subksps[i]);
        if (this->_preconditioner_type == ILU_PRECOND)
          PCSetType(subpcs, (char*) PCILU);
        else
          PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);
        PCFactorSetZeroPivot(subpcs, epsilon);
        PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
      }
    }
    else {
      for (int i = 0; i < _nlocal; i++) {
        PC subpcs;
        KSPGetPC(subksps[i], &subpcs);
        KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
        KSPSetFromOptions(subksps[i]);
        if (this->_preconditioner_type == ILU_PRECOND)
          PCSetType(subpcs, (char*) PCILU);
        else
          PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);
        PCFactorSetZeroPivot(subpcs, epsilon);
        PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
      }
    }

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
        PCMGGetSmootherUp(pcMG, level , &subkspUp);
        KSPSetTolerances(subkspUp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, npost);
        this->set_petsc_solver_type(subkspUp);
        KSPSetPC(subkspUp, subpc);
        PC subpcUp;
        KSPGetPC(subkspUp, &subpcUp);
        //PCASMSetLocalType(subpcUp, PC_COMPOSITE_MULTIPLICATIVE);
        KSPSetUp(subkspUp);
      }
    }
  }


  void AsmPetscLinearEquationSolver::MGsolve(const bool ksp_clean) {

    PetscLogDouble t1; PetscLogDouble t2;
    PetscTime(&t1);

    if (ksp_clean) {
      PetscMatrix* KKp = static_cast< PetscMatrix* >(_KK);
      Mat KK = KKp->mat();
      KSPSetOperators(_ksp, KK, _Pmat);
      KSPSetTolerances(_ksp, _rtol, _abstol, _dtol, _maxits);
      KSPSetInitialGuessKnoll(_ksp, PETSC_TRUE);
      KSPSetFromOptions(_ksp);

      KSPGMRESSetRestart(_ksp, _restart);
      KSPSetUp(_ksp);

//       PetscReal rnorm1;
//       KSPGetResidualNorm(_ksp, &rnorm1);
//       PetscPrintf(PETSC_COMM_WORLD, " Residual norm = %10.8g \n", rnorm1);
    }

    PetscVector* EPSCp = static_cast< PetscVector* >(_EPSC);
    Vec EPSC = EPSCp->vec();
    PetscVector* RESp = static_cast< PetscVector* >(_RES);
    Vec RES = RESp->vec();

    KSPSolve(_ksp, RES, EPSC);

    _RESC->matrix_mult(*_EPSC, *_KK);

    *_RES -= *_RESC;
    *_EPS += *_EPSC;

//#ifndef NDEBUG
    int its;
    KSPGetIterationNumber(_ksp, &its);

    KSPConvergedReason reason;
    KSPGetConvergedReason(_ksp, &reason);

    PetscReal rnorm;
    KSPGetResidualNorm(_ksp, &rnorm);

    PetscTime(&t2);
    PetscPrintf(PETSC_COMM_WORLD, " *************** MG linear solver time: %8.3f \n", t2 -t1);
    PetscPrintf(PETSC_COMM_WORLD, " *************** Number of outer ksp solver iterations = %i \n", its);
    PetscPrintf(PETSC_COMM_WORLD, " *************** Convergence reason = %i \n", reason);
    PetscPrintf(PETSC_COMM_WORLD, " *************** Residual norm = %10.8g \n", rnorm);

//#endif

  }

// ================================================

  void AsmPetscLinearEquationSolver::clear() {
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

  void AsmPetscLinearEquationSolver::init(Mat& Amat, Mat& Pmat) {

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
      ierr = KSPSetTolerances(_ksp, _rtol, _abstol, _dtol, _maxits); CHKERRABORT(MPI_COMM_WORLD, ierr);

      if ( _msh->GetLevel() != 0)
        KSPSetInitialGuessKnoll(_ksp, PETSC_TRUE);

      if (_msh->GetLevel() != 0)
        KSPSetNormType(_ksp, KSP_NORM_NONE);

      ierr = KSPSetFromOptions(_ksp);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      ierr = KSPGMRESSetRestart(_ksp, _restart); CHKERRABORT(MPI_COMM_WORLD, ierr);


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

  void AsmPetscLinearEquationSolver::set_petsc_solver_type(KSP& ksp) {
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
	ierr =  KSPRichardsonSetScale(ksp, 0.7);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        KSPRichardsonSetScale(ksp, 0.7);
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

