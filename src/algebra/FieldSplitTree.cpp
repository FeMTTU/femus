/*=========================================================================

  Program: FEMUS
  Module: FieldSplitPetscLinearEquationSolver
  Authors: Eugenio Aulisa, Guoyi Ke

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

#include <map>
#include <algorithm>

#include "FieldSplitTree.hpp"
#include "FieldSplitPetscLinearEquationSolver.hpp"
#include "MeshASMPartitioning.hpp"

namespace femus {

  FieldSplitTree::FieldSplitTree(const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, const std::vector < unsigned > & solutionType, std::string name) {
    _solutionType = solutionType;

    FieldSplitTreeBuild(solver, preconditioner, fields, name);
  }

  FieldSplitTree::FieldSplitTree(const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, std::string name) {
    FieldSplitTreeBuild(solver, preconditioner, fields, name);
  }


  void FieldSplitTree::FieldSplitTreeBuild(const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, std::string name) {
    _father = NULL;
    _name = name;

    _solver = solver;
    _preconditioner = preconditioner;
    _numberOfSplits = 1;
    _child.resize(0);


    _fieldsSplit.resize(1);
    _fieldsSplit[0] = fields;

    //BEGIN ALL FIELD COLLECTION
    std::map < unsigned, bool > mymap;

    for(unsigned i = 0; i < _numberOfSplits; i++) {
      for(unsigned j = 0; j < _fieldsSplit[i].size(); j++) {
        mymap[_fieldsSplit[i][j]] = true;
      }
    }

    _fieldsAll.resize(mymap.size());

    unsigned j = 0;

    for(std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it) {
      _fieldsAll[j] = it->first;
      j++;
    }

    //END ALL FIELD COLLECTION
    _rtol = 1.e-3;
    _abstol = 1.e-20;
    _dtol = 1.e+50;
    _maxits = 1;
    _schurFactType = SCHUR_FACT_AUTOMATIC;
    _schurPreType = SCHUR_PRE_AUTOMATIC;

    if(preconditioner == ASM_PRECOND && _solutionType.size() != fields.size()) {
      std::cout << "Error! The Solution type with PCFieldSplit - ASM preconditioner has to be specified" << std::endl;
      abort();
    }

    //Only for ASM preconditioner
    _asmLocalIs.reserve(10);
    _asmOverlappingIs.reserve(10);
    _asmLocalIsIndex.reserve(10);
    _asmOverlappingIsIndex.reserve(10);
    _asmBlockSize = 100;
    _asmStandard = true;
    _asmOverlapping = 1;
    _asmSchurVariableNumber = 1;

  };

  //multiple split constructor
  FieldSplitTree::FieldSplitTree(const SolverType& solver, const PreconditionerType& preconditioner, std::vector < FieldSplitTree*> childBranch, std::string name) {

    _father = NULL;

    _name = name;
    _solver = solver;
    _preconditioner = preconditioner;
    _numberOfSplits = childBranch.size();
    _fieldsSplit.resize(_numberOfSplits);
    _child.resize(_numberOfSplits);

    for(unsigned i = 0; i < _numberOfSplits; i++) {
      childBranch[i]->_father = this;
      _child[i] = childBranch[i];
      _fieldsSplit[i].resize(childBranch[i]->_fieldsAll.size());

      for(unsigned j = 0; j < childBranch[i]->_fieldsAll.size(); j++) {
        _fieldsSplit[i][j] = childBranch[i]->_fieldsAll[j];
      }
    }

    //BEGIN ALL FIELD COLLECTION
    std::map < unsigned, bool > mymap;

    for(unsigned i = 0; i < _numberOfSplits; i++) {
      for(unsigned j = 0; j < _fieldsSplit[i].size(); j++) {
        mymap[_fieldsSplit[i][j]] = true;
      }
    }

    _fieldsAll.resize(mymap.size());

    unsigned j = 0;

    for(std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it) {
      _fieldsAll[j] = it->first;
      j++;
    }

    //END ALL FIELD COLLECTION
    _rtol = 1.e-3;
    _abstol = 1.e-20;
    _dtol = 1.e+50;
    _maxits = 1;
    _schurFactType = SCHUR_FACT_AUTOMATIC;
    _schurPreType = SCHUR_PRE_AUTOMATIC;

    //Only for ASM preconditioner
    _asmLocalIs.reserve(10);
    _asmOverlappingIs.reserve(10);
    _asmLocalIsIndex.reserve(10);
    _asmOverlappingIsIndex.reserve(10);
    _asmBlockSize = 100;
    _asmStandard = true;
    _asmOverlapping = 1;
    _asmSchurVariableNumber = 1;
  }


  FieldSplitTree::~FieldSplitTree() {
    if(_numberOfSplits > 1) {
      for(unsigned i = 0; i < _isSplit.size(); i++) {
        for(unsigned j = 0; j < _isSplit[i].size(); j++) {
          ISDestroy(&_isSplit[i][j]);
        }
      }
    }

    for(unsigned i = 0; i < _asmLocalIs.size(); i++) {
      for(unsigned j = 0; j < _asmLocalIs[i].size(); j++) {
        ISDestroy(&_asmLocalIs[i][j]);
        ISDestroy(& _asmOverlappingIs[i][j]);
      }
    }


    for(unsigned i = 0; i < _isSplitIndexPt.size(); i++) {
      delete [] _isSplitIndexPt[i];
    }
  }

  void FieldSplitTree::PrintFieldSplitTree(const unsigned& counter) {

    std::string sub = " ";

    for(int i = 0; i < counter; i++) {
      sub += "sub-";
    }

    std::cout << "Fields in the " << _name << sub << "system:\n";

    for(int j = 0; j < _fieldsAll.size(); j++) std::cout << _fieldsAll[j] << " ";

    std::cout << std::endl;

    if(_father != NULL) {
      std::cout << "My father is " << _father->GetName() << std::endl << std::endl;
    }

    if(GetNumberOfSplits() > 1) {
      for(unsigned i = 0; i < GetNumberOfSplits(); i++) {
        _child[i]->PrintFieldSplitTree(counter + 1);
      }
    }

  }


  void FieldSplitTree::BuildIndexSet(const std::vector< std::vector < unsigned > >& KKoffset, const unsigned& iproc,
                                     const unsigned& nprocs, const unsigned& level, const FieldSplitPetscLinearEquationSolver *solver) {

    if(_MatrixOffset.size() < level) _MatrixOffset.resize(level);

    _MatrixOffset[level - 1] = KKoffset;

    if(GetNumberOfSplits() == 1) {
      if(_preconditioner == ASM_PRECOND && !_asmStandard) {
        BuildASMIndexSet(level, solver);
      }

      return;
    }

    if(_isSplit.size() < level) _isSplit.resize(level);

    _isSplit[level - 1].resize(GetNumberOfSplits());

    for(unsigned i = 0; i < GetNumberOfSplits(); i++) {

      //on the actual structure
      unsigned size = 0;

      for(unsigned k = 0; k < _fieldsSplit[i].size(); k++) {
        unsigned index = _fieldsSplit[i][k];
        unsigned offset = KKoffset[index][iproc];
        unsigned offsetp1 = KKoffset[index + 1][iproc];
        size += offsetp1 - offset;
      }


      PetscInt* isSplitIndex = new PetscInt [size];
      unsigned ptSize = _isSplitIndexPt.size();
      _isSplitIndexPt.resize(ptSize + 1);
      _isSplitIndexPt[ptSize] = isSplitIndex;

      unsigned counter = 0;

      for(int k = 0; k < _fieldsSplit[i].size(); k++) {
        unsigned index = _fieldsSplit[i][k];
        unsigned offset = KKoffset[index][iproc];
        unsigned offsetp1 = KKoffset[index + 1][iproc];

        for(int j = offset; j < offsetp1; j++) {
          isSplitIndex[counter] = j;
          counter++;
        }
      }

      ISCreateGeneral(MPI_COMM_WORLD, size, isSplitIndex, PETSC_USE_POINTER, &_isSplit[level - 1][i]);

      // on the child branches

      std::vector < std::vector < unsigned > > tempFields = _child[i]->_fieldsSplit;

      for(unsigned j = 0; j < _child[i]->_fieldsAll.size(); j++) {
        unsigned index = _child[i]->_fieldsAll[j];
        _child[i]->_fieldsAll[j] = j;

        for(unsigned k = 0; k < _child[i]->_numberOfSplits; k++) {
          for(unsigned l = 0; l < _child[i]->_fieldsSplit[k].size(); l++) {
            if(tempFields[k][l] == index) {
              _child[i]->_fieldsSplit[k][l] = j;

            }
          }
        }
      }

      std::vector< std::vector < unsigned > > fieldsInSplitOffset(_fieldsSplit[i].size() + 1);

      for(unsigned j = 0; j < _fieldsSplit[i].size() + 1; j++) {
        fieldsInSplitOffset[j].resize(nprocs);
      }

      fieldsInSplitOffset[0][0] = 0;

      for(unsigned jproc = 0; jproc < nprocs; jproc++) {
        for(unsigned k = 0; k < _fieldsSplit[i].size(); k++) {
          unsigned index = _fieldsSplit[i][k];
          unsigned ownedDofs =  KKoffset[index + 1][jproc] -  KKoffset[index][jproc];
          fieldsInSplitOffset[k + 1][jproc] = fieldsInSplitOffset[k][jproc] + ownedDofs;
        }

        if(jproc < nprocs - 1) {
          fieldsInSplitOffset[0][jproc + 1] = fieldsInSplitOffset[_fieldsSplit[i].size()][jproc];
        }
      }

      _child[i]->BuildIndexSet(fieldsInSplitOffset, iproc, nprocs, level, solver);
    }

  }

  /*---------adjusted by Guoyi Ke-----------*/
  void FieldSplitTree::SetupKSPTolerances(const double& rtol, const double& abstol, const double& dtol, const unsigned& maxits) {
    _rtol = rtol;
    _abstol = abstol;
    _dtol = dtol;
    _maxits = maxits;
  }
  void FieldSplitTree::SetupSchurFactorizationType(const SchurFactType& schurFactType) {
    _schurFactType = schurFactType;
  }
  void FieldSplitTree::SetupSchurPreType(const SchurPreType& schurPreType) {
    _schurPreType = schurPreType;
  }
  /*---------adjusted by Guoyi Ke-----------*/
  void FieldSplitTree::SetPC(KSP& ksp, const unsigned& level) {

//     std::cout<<level<<std::endl;
//     for(int i = 0; i < _MatrixOffset[level-1].size(); i++){
//       for(int j = 0; j < _MatrixOffset[level-1][i].size(); j++){
// 	std::cout << _MatrixOffset[level-1][i][j] << " ";
//       }
//       std::cout <<std::endl;
//     }
//     std::cout <<std::endl;

    PC pc;
    KSPGetPC(ksp, &pc);

    //BEGIN from here
    if(_preconditioner == FIELDSPLIT_PRECOND) {
      PetscPreconditioner::set_petsc_preconditioner_type(_preconditioner, pc);
      PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);

      //PCFieldSplitSetType( pc, PC_COMPOSITE_MULTIPLICATIVE );
      for(unsigned i = 0; i < _numberOfSplits; i++) {
        PCFieldSplitSetIS(pc, NULL, _isSplit[level - 1][i]);
      }

      PCSetUp(pc);
      KSP* subksp;
      PetscInt nlocal = static_cast < PetscInt >(_numberOfSplits);
      PCFieldSplitGetSubKSP(pc, &nlocal, &subksp);

      for(unsigned i = 0; i < _numberOfSplits; i++) {
        _child[i]->SetPC(subksp[i], level);
      }

      PetscFree(subksp);
    }

    else if(_preconditioner == ASM_PRECOND) {
      SetPetscSolverType(ksp);
      KSPSetTolerances(ksp, _rtol, _abstol, _dtol, _maxits);
      
//       PetscPreconditioner::set_petsc_preconditioner_type(FIELDSPLIT_PRECOND, pc);
//       
//       if(!_asmStandard) {
// 	for(unsigned i=0;i<_asmLocalIsIndex[level - 1].size();i++){
// 	  PCFieldSplitSetIS(pc,NULL,_asmOverlappingIs[level - 1][i]);
// 	}
//       }
// 
//       //PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);
//       KSPSetUp(ksp);
// 
//       KSP* subksps;
//       PetscInt nlocal;
//       PCFieldSplitGetSubKSP(pc, &nlocal, &subksps);
      
      
      PetscPreconditioner::set_petsc_preconditioner_type(ASM_PRECOND, pc);

      if(!_asmStandard) {
        PCASMSetLocalSubdomains(pc, _asmLocalIsIndex[level - 1].size(), &_asmOverlappingIs[level - 1][0], &_asmLocalIs[level - 1][0]);
      }

      PCASMSetOverlap(pc, _asmOverlapping);
      KSPSetUp(ksp);

      //PCASMSetLocalType(pc, PC_COMPOSITE_MULTIPLICATIVE);

      KSP* subksps;
      PetscInt nlocal;
      PCASMGetSubKSP(pc, &nlocal, PETSC_NULL, &subksps);

      PetscReal epsilon = 1.e-16;

      if(!_asmStandard) {
        for(int i = 0; i < _asmBlockMaterialRange[level - 1][0]; i++) {
          PC subpcs;
          KSPGetPC(subksps[i], &subpcs);
          KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
          KSPSetFromOptions(subksps[i]);
          PetscPreconditioner::set_petsc_preconditioner_type(MLU_PRECOND, subpcs);
          PCFactorSetZeroPivot(subpcs, epsilon);
          PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
        }

        for(int i = _asmBlockMaterialRange[level - 1][0]; i < _asmBlockMaterialRange[level - 1][1]; i++) {
          PC subpcs;
          KSPGetPC(subksps[i], &subpcs);
          KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
          KSPSetFromOptions(subksps[i]);
          PCSetType(subpcs, (char*) PCILU);
          PCFactorSetZeroPivot(subpcs, epsilon);
          PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
        }
      }
      else {
        for(int i = 0; i < nlocal; i++) {
          PC subpcs;
          KSPGetPC(subksps[i], &subpcs);
          KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
          KSPSetFromOptions(subksps[i]);
          PCSetType(subpcs, (char*) PCILU);
          PCFactorSetZeroPivot(subpcs, epsilon);
          PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
        }
      }
    }
    else if(_preconditioner == FS_SCHUR_PRECOND) {

      PetscPreconditioner::set_petsc_preconditioner_type(FIELDSPLIT_PRECOND, pc);

      PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);

      SetSchurFactorizationType(pc);
      SetSchurPreType(pc);

      for(int i = 0; i < _numberOfSplits; i++) {
        PCFieldSplitSetIS(pc, NULL, _isSplit[level - 1][i]);
      }

      PCSetUp(pc);

      KSP* subksp;
      PetscInt nlocal = static_cast < PetscInt >(_numberOfSplits);
      PCFieldSplitGetSubKSP(pc, &nlocal, &subksp);

      for(unsigned i = 0; i < _numberOfSplits; i++) {
        _child[i]->SetPC(subksp[i], level);
      }

      PetscFree(subksp);
    }
    else {
      if(_preconditioner == LSC_PRECOND) {
        if(_solver == PREONLY) {
          std::cout << "LSC does not allow to use PREONLY, and ksp switches to GMRES instead!" << std::endl;
          _solver = GMRES;
        }
      }

      SetPetscSolverType(ksp);
      PC pc;
      KSPGetPC(ksp, &pc);
      KSPSetTolerances(ksp, _rtol, _abstol, _dtol, _maxits);
      KSPSetFromOptions(ksp);
      PetscReal epsilon = 1.e-16;

      PetscPreconditioner::set_petsc_preconditioner_type(_preconditioner, pc);
      PCFactorSetZeroPivot(pc, epsilon);
      PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
    }
  }

  FieldSplitTree* FieldSplitTree::GetFather() const {
    if(_father != NULL) {
      return _father;
    }
    else {
      std::cout << "Warning this split has no father" << std::endl;
      abort();
    }
  }

  FieldSplitTree* FieldSplitTree::GetChild(const unsigned& i) {
    if(i < _numberOfSplits) {
      return _child[i];
    }
    else {
      std::cout << "Wrong input (= " << i << ") in function FieldSplitTree::GetChildrenBranchstd" << std::cout;
      std::cout << "Number of Splits = " << _numberOfSplits << std::endl;
      abort();
    }

  }

  void FieldSplitTree::SetPetscSolverType(KSP& ksp) {
    int ierr = 0;

    switch(_solver) {
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
                  << _solver               << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }
  }


  void FieldSplitTree::SetSchurFactorizationType(PC &pc) {

    switch(_schurFactType) {
      case  SCHUR_FACT_UPPER:
        PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER);
        return;

      case  SCHUR_FACT_LOWER:
        PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER);
        return;

      case SCHUR_FACT_FULL:
        PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL);
        return;

      case SCHUR_FACT_AUTOMATIC:
        PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER);

        for(int i = 0; i < _numberOfSplits; i++) {
          if(GetChild(i)->_preconditioner == LSC_PRECOND) { //it goes with pressure LSC
            PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL);
          }
        }

        return;

      default:
        std::cerr << "ERROR:  Unsupported Schur Factorization Type: "
                  << _schurFactType << std::endl;
        abort();
    }
  }

  void FieldSplitTree::SetSchurPreType(PC &pc) {

    switch(_schurPreType) {
      case  SCHUR_PRE_SELF:
        PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELF, NULL);
        return;

      case  SCHUR_PRE_SELFP:
        PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
        return;

      case SCHUR_PRE_USER:
        PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, NULL);
        return;

      case SCHUR_PRE_A11:
        PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_A11, NULL);
        return;

      case SCHUR_PRE_FULL:
        PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_FULL, NULL);
        return;

      case SCHUR_PRE_AUTOMATIC:
        PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);

        for(int i = 0; i < _numberOfSplits; i++) {
          if(GetChild(i)->_preconditioner == LSC_PRECOND) { //it goes with pressure LSC
            PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELF, NULL);
          }
        }

        return;

      default:
        std::cerr << "ERROR:  Unsupported Schur Precondition Type: "
                  << _schurPreType << std::endl;
        abort();
    }
  }





  void FieldSplitTree::BuildASMIndexSet(const unsigned& level, const FieldSplitPetscLinearEquationSolver *solver) {


    Mesh* msh = solver->_msh;

    unsigned nel = msh->GetNumberOfElements();

    bool FastVankaBlock = true;

    if(_asmSchurVariableNumber != 0) {
      FastVankaBlock = (_solutionType[ _solutionType.size() - _asmSchurVariableNumber ] < 3) ? false : true;
    }

    unsigned iproc = solver->processor_id();

    unsigned DofOffset = _MatrixOffset[level - 1][0][iproc];
    unsigned DofOffsetSize = _MatrixOffset[level - 1][_solutionType.size()][iproc] - DofOffset;

    vector < unsigned > indexa(DofOffsetSize, DofOffsetSize);
    vector < unsigned > indexb(DofOffsetSize, DofOffsetSize);

    vector <bool> owned(DofOffsetSize, false);

    std::map<int, bool> mymap;

    unsigned ElemOffset   = msh->_dofOffset[3][iproc];
    unsigned ElemOffsetp1 = msh->_dofOffset[3][iproc + 1];
    unsigned ElemOffsetSize = ElemOffsetp1 - ElemOffset;
    vector <PetscInt> indexci(ElemOffsetSize);
    vector < unsigned > indexc(ElemOffsetSize, ElemOffsetSize);

    vector < vector < unsigned > > block_elements;

    unsigned base = pow(2, msh->GetDimension());
    unsigned elementBlockNumber[2];
    elementBlockNumber[0] = pow(base, _asmBlockSize);
    elementBlockNumber[1] = pow(base, _asmBlockSize);


    elementBlockNumber[0] = 256;
    elementBlockNumber[1] = 256;

    _asmBlockMaterialRange.resize(level);

    MeshASMPartitioning meshasmpartitioning(*msh);

    meshasmpartitioning.DoPartition(elementBlockNumber, block_elements, _asmBlockMaterialRange[level - 1]);

    vector <bool> ThisVaribaleIsNonSchur(_solutionType.size(), true);

    for(unsigned i = _solutionType.size() - _asmSchurVariableNumber; i < _solutionType.size(); i++) {
      ThisVaribaleIsNonSchur[i] = false;
    }

    // *** Start Vanka Block ***

    _asmLocalIsIndex.resize(level);
    _asmOverlappingIsIndex.resize(level);

    _asmLocalIsIndex[level - 1].resize(block_elements.size());
    _asmOverlappingIsIndex[level - 1].resize(block_elements.size());

    for(int vb_index = 0; vb_index < block_elements.size(); vb_index++) { //loop on the vanka-blocks
      _asmLocalIsIndex[level - 1][vb_index].resize(DofOffsetSize);
      _asmOverlappingIsIndex[level - 1][vb_index].resize(DofOffsetSize);

      PetscInt PAsize = 0;
      PetscInt PBsize = 0;

      PetscInt Csize = 0;

      // ***************** NODE/ELEMENT SERCH *******************
      for(int kel = 0; kel < block_elements[vb_index].size(); kel++) { //loop on the vanka-block elements
        unsigned iel = block_elements[vb_index][kel];
	for(unsigned j = 0; j < msh->el->GetElementNearElementSize(iel,!FastVankaBlock);j++){
	  unsigned jel = msh->el->GetElementNearElement(iel,j);
	  if( jel >= ElemOffset && jel<ElemOffsetp1 ){
	    if(indexc[jel - ElemOffset] == ElemOffsetSize) {
              indexci[Csize] = jel - ElemOffset;
              indexc[jel - ElemOffset] = Csize++;

              //add non-schur variables to be solved
              for(int indexSol = 0; indexSol < _solutionType.size(); indexSol++) {
                if(ThisVaribaleIsNonSchur[indexSol]) {
                  unsigned SolType = _solutionType[indexSol];
                  unsigned nvej = msh->GetElementDofNumber(jel, SolType);

                  for(unsigned jj = 0; jj < nvej; jj++) {
                    unsigned jdof = msh->GetSolutionDof(jj, jel, SolType);

                    unsigned kkdof = solver->GetSystemDof(SolType, indexSol, jj, jel, _MatrixOffset[level - 1]);
                  
                    if(jdof >= msh->_dofOffset[SolType][iproc] &&
                        jdof <  msh->_dofOffset[SolType][iproc + 1]) {
                      if(indexa[kkdof - DofOffset] == DofOffsetSize && owned[kkdof - DofOffset] == false) {
                        owned[kkdof - DofOffset] = true;
                        _asmLocalIsIndex[level - 1][vb_index][PAsize] = kkdof;
                        indexa[kkdof - DofOffset] = PAsize++;
                      }

                      if(indexb[kkdof - DofOffset] == DofOffsetSize) {
                        _asmOverlappingIsIndex[level - 1][vb_index][PBsize] = kkdof;
                        indexb[kkdof - DofOffset] = PBsize++;
                      }
                    }
                    else {
                      mymap[kkdof] = true;
                    }
                  }
                }
              }
            }
          }
        }

        //-----------------------------------------------------------------------------------------
        //Add Schur nodes (generally pressure type variables) to be solved
        {
          for(int indexSol = 0; indexSol < _solutionType.size(); indexSol++) {
            if(!ThisVaribaleIsNonSchur[indexSol]) {
              unsigned SolPdeIndex = indexSol;
              unsigned SolType = _solutionType[SolPdeIndex];
              unsigned nvei = msh->GetElementDofNumber(iel, SolType);

              for(unsigned ii = 0; ii < nvei; ii++) {
                unsigned inode_Metis = msh->GetSolutionDof(ii, iel, SolType);

                unsigned kkdof = solver->GetSystemDof(SolType, indexSol, ii, iel, _MatrixOffset[level - 1]);

                if(inode_Metis >= msh->_dofOffset[SolType][iproc] &&
                    inode_Metis <  msh->_dofOffset[SolType][iproc + 1]) {
                  if(indexa[kkdof - DofOffset] == DofOffsetSize && owned[kkdof - DofOffset] == false) {
                    owned[kkdof - DofOffset] = true;
                    _asmLocalIsIndex[level - 1][vb_index][PAsize] = kkdof;
                    indexa[kkdof - DofOffset] = PAsize++;
                  }

                  if(indexb[kkdof - DofOffset] == DofOffsetSize) {
                    _asmOverlappingIsIndex[level - 1][vb_index][PBsize] = kkdof;
                    indexb[kkdof - DofOffset] = PBsize++;
                  }
                }
                else {
                  mymap[kkdof] = true;
                }
              }
            }
          }
        }
        //-----------------------------------------------------------------------------------------
      }

      // *** re-initialize indeces(a,c,d)
      for(PetscInt i = 0; i < PAsize; i++) {
        indexa[_asmLocalIsIndex[level - 1][vb_index][i] - DofOffset] = DofOffsetSize;
      }

      for(PetscInt i = 0; i < PBsize; i++) {
        indexb[_asmOverlappingIsIndex[level - 1][vb_index][i] - DofOffset] = DofOffsetSize;
      }

      for(PetscInt i = 0; i < Csize; i++) {
        indexc[indexci[i]] = ElemOffsetSize;
      }

      _asmLocalIsIndex[level - 1][vb_index].resize(PAsize);
      std::vector < PetscInt >(_asmLocalIsIndex[level - 1][vb_index]).swap(_asmLocalIsIndex[level - 1][vb_index]);

      _asmOverlappingIsIndex[level - 1][vb_index].resize(PBsize + mymap.size());
      int i = 0;

      for(std::map<int, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it, ++i) {
        _asmOverlappingIsIndex[level - 1][vb_index][PBsize + i] = it->first;
      }

      std::vector < PetscInt >(_asmOverlappingIsIndex[level - 1][vb_index]).swap(_asmOverlappingIsIndex[level - 1][vb_index]);


      mymap.clear();

      std::sort(_asmLocalIsIndex[level - 1][vb_index].begin(), _asmLocalIsIndex[level - 1][vb_index].end());
      std::sort(_asmOverlappingIsIndex[level - 1][vb_index].begin(), _asmOverlappingIsIndex[level - 1][vb_index].end());


    }

    //BEGIN Generate std::vector<IS> for ASM PC ***********
    _asmLocalIs.resize(level);
    _asmOverlappingIs.resize(level);
    _asmLocalIs[level - 1].resize(_asmLocalIsIndex[level - 1].size());
    _asmOverlappingIs[level - 1].resize(_asmOverlappingIsIndex[level - 1].size());

    for(unsigned vb_index = 0; vb_index < _asmLocalIsIndex[level - 1].size(); vb_index++) {
      ISCreateGeneral(MPI_COMM_SELF, _asmLocalIsIndex[level - 1][vb_index].size(), &_asmLocalIsIndex[level - 1][vb_index][0], PETSC_USE_POINTER, &_asmLocalIs[level - 1][vb_index]);
      ISCreateGeneral(MPI_COMM_SELF, _asmOverlappingIsIndex[level - 1][vb_index].size(), &_asmOverlappingIsIndex[level - 1][vb_index][0], PETSC_USE_POINTER, &_asmOverlappingIs[level - 1][vb_index]);
    }

    //END Generate std::vector<IS> for ASM PC ***********
    return;

  }

}

