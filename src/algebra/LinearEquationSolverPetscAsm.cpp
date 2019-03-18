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
#include "LinearEquationSolverPetscAsm.hpp"
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

  // ==============================================

  void LinearEquationSolverPetscAsm::SetElementBlockNumber(const char all[], const unsigned& overlap) {
    _elementBlockNumber[0] = _msh->GetNumberOfElements();
    _elementBlockNumber[1] = _msh->GetNumberOfElements();
    _elementBlockNumber[2] = _msh->GetNumberOfElements();
    _standardASM = 1;
    _overlap = overlap;
  }

  // =================================================

  void LinearEquationSolverPetscAsm::SetElementBlockNumber(const unsigned& block_elemet_number) {
    _elementBlockNumber[0] = block_elemet_number;
    _elementBlockNumber[1] = block_elemet_number;
    _elementBlockNumber[2] = block_elemet_number;
    _bdcIndexIsInitialized = 0;
    _standardASM = 0;
  }

  // =================================================

  void LinearEquationSolverPetscAsm::SetElementBlockNumberSolid(const unsigned& block_elemet_number, const unsigned& overlap) {
    _elementBlockNumber[0] = block_elemet_number;
    _bdcIndexIsInitialized = 0;
    _standardASM = 0;
    _overlap = overlap;
  }

  // =================================================

  void LinearEquationSolverPetscAsm::SetElementBlockNumberFluid(const unsigned& block_elemet_number, const unsigned& overlap) {
    _elementBlockNumber[2] = block_elemet_number;
    _bdcIndexIsInitialized = 0;
    _standardASM = 0;
    _overlap = overlap;
  }
  
  void LinearEquationSolverPetscAsm::SetElementBlockNumberPorous(const unsigned& block_elemet_number, const unsigned& overlap) {
    _elementBlockNumber[1] = block_elemet_number;
    _bdcIndexIsInitialized = 0;
    _standardASM = 0;
    _overlap = overlap;
  }

  // ==============================================

  void LinearEquationSolverPetscAsm::BuildAMSIndex(const vector <unsigned>& variable_to_be_solved) {

    unsigned nel = _msh->GetNumberOfElements();

    bool FastVankaBlock = true;

    if(_NSchurVar != 0) {
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

    meshasmpartitioning.DoPartition(_elementBlockNumber, block_elements, _blockTypeRange);

    vector <bool> ThisVaribaleIsNonSchur(_SolPdeIndex.size(), true);

    for(unsigned iind = variable_to_be_solved.size() - _NSchurVar; iind < variable_to_be_solved.size(); iind++) {
      unsigned PdeIndexSol = variable_to_be_solved[iind];
      ThisVaribaleIsNonSchur[PdeIndexSol] = false;
    }

    // *** Start Vanka Block ***

    _localIsIndex.resize(block_elements.size());
    _overlappingIsIndex.resize(block_elements.size());

    for(int vb_index = 0; vb_index < block_elements.size(); vb_index++) { //loop on the vanka-blocks
      _localIsIndex[vb_index].resize(DofOffsetSize);
      _overlappingIsIndex[vb_index].resize(DofOffsetSize);

      PetscInt PAsize = 0;
      PetscInt PBsize = 0;

      PetscInt Csize = 0;

      // ***************** NODE/ELEMENT SERCH *******************
      for(int kel = 0; kel < block_elements[vb_index].size(); kel++) { //loop on the vanka-block elements
        unsigned iel = block_elements[vb_index][kel];
	for(unsigned j = 0; j < _msh->el->GetElementNearElementSize(iel,!FastVankaBlock);j++){
	  unsigned jel = _msh->el->GetElementNearElement(iel,j);
	  if( jel >= ElemOffset && jel<ElemOffsetp1 ){
	    
            //add elements for velocity to be solved
            if(indexc[jel - ElemOffset] == ElemOffsetSize) {
              indexci[Csize] = jel - ElemOffset;
              indexc[jel - ElemOffset] = Csize++;

              //add non-schur variables to be solved
              for(int indexSol = 0; indexSol < _SolPdeIndex.size(); indexSol++) {
                if(ThisVaribaleIsNonSchur[indexSol]) {
                  unsigned SolPdeIndex = _SolPdeIndex[indexSol];
                  unsigned SolType = _SolType[SolPdeIndex];
                  unsigned nvej = _msh->GetElementDofNumber(jel, SolType);

                  for(unsigned jj = 0; jj < nvej; jj++) {
                    unsigned jdof = _msh->GetSolutionDof(jj, jel, SolType);
                    unsigned kkdof = GetSystemDof(SolPdeIndex, indexSol, jj, jel);

                    if(jdof >= _msh->_dofOffset[SolType][iproc] &&
                        jdof <  _msh->_dofOffset[SolType][iproc + 1]) {
                      if(indexa[kkdof - DofOffset] == DofOffsetSize && owned[kkdof - DofOffset] == false) {
                        owned[kkdof - DofOffset] = true;
                        _localIsIndex[vb_index][PAsize] = kkdof;
                        indexa[kkdof - DofOffset] = PAsize++;
                      }

                      if(indexb[kkdof - DofOffset] == DofOffsetSize) {
                        _overlappingIsIndex[vb_index][PBsize] = kkdof;
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
          for(int indexSol = 0; indexSol < _SolPdeIndex.size(); indexSol++) {
            if(!ThisVaribaleIsNonSchur[indexSol]) {
              unsigned SolPdeIndex = _SolPdeIndex[indexSol];
              unsigned SolType = _SolType[SolPdeIndex];
              unsigned nvei = _msh->GetElementDofNumber(iel, SolType);

              for(unsigned ii = 0; ii < nvei; ii++) {
                unsigned inode_Metis = _msh->GetSolutionDof(ii, iel, SolType);
                unsigned kkdof = GetSystemDof(SolPdeIndex, indexSol, ii, iel);

                if(inode_Metis >= _msh->_dofOffset[SolType][iproc] &&
                    inode_Metis <  _msh->_dofOffset[SolType][iproc + 1]) {
                  if(indexa[kkdof - DofOffset] == DofOffsetSize && owned[kkdof - DofOffset] == false) {
                    owned[kkdof - DofOffset] = true;
                    _localIsIndex[vb_index][PAsize] = kkdof;
                    indexa[kkdof - DofOffset] = PAsize++;
                  }

                  if(indexb[kkdof - DofOffset] == DofOffsetSize) {
                    _overlappingIsIndex[vb_index][PBsize] = kkdof;
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
        indexa[_localIsIndex[vb_index][i] - DofOffset] = DofOffsetSize;
      }

      for(PetscInt i = 0; i < PBsize; i++) {
        indexb[_overlappingIsIndex[vb_index][i] - DofOffset] = DofOffsetSize;
      }

      for(PetscInt i = 0; i < Csize; i++) {
        indexc[indexci[i]] = ElemOffsetSize;
      }

      _localIsIndex[vb_index].resize(PAsize);
      std::vector < PetscInt >(_localIsIndex[vb_index]).swap(_localIsIndex[vb_index]);

      _overlappingIsIndex[vb_index].resize(PBsize + mymap.size());
      int i = 0;
      for(std::map<int, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it, ++i) {
        _overlappingIsIndex[vb_index][PBsize + i] = it->first;
      }
      std::vector < PetscInt >(_overlappingIsIndex[vb_index]).swap(_overlappingIsIndex[vb_index]);


      mymap.clear();

      std::sort(_localIsIndex[vb_index].begin(), _localIsIndex[vb_index].end());
      std::sort(_overlappingIsIndex[vb_index].begin(), _overlappingIsIndex[vb_index].end());


    }

    //BEGIN Generate std::vector<IS> for ASM PC ***********
    _localIs.resize(_localIsIndex.size());
    _overlappingIs.resize(_overlappingIsIndex.size());

    for(unsigned vb_index = 0; vb_index < _localIsIndex.size(); vb_index++) {
      ISCreateGeneral(MPI_COMM_SELF, _localIsIndex[vb_index].size(), &_localIsIndex[vb_index][0],  PETSC_USE_POINTER , &_localIs[vb_index]);
      ISCreateGeneral(MPI_COMM_SELF, _overlappingIsIndex[vb_index].size(), &_overlappingIsIndex[vb_index][0],  PETSC_USE_POINTER , &_overlappingIs[vb_index]);
    }

    //END Generate std::vector<IS> for ASM PC ***********

    return;
  }

  // =================================================

  void LinearEquationSolverPetscAsm::SetPreconditioner(KSP& subksp, PC& subpc) {
    
    PetscPreconditioner::set_petsc_preconditioner_type(ASM_MULTIPLICATIVE_PRECOND, subpc);

    if(!_standardASM) {
      PCASMSetLocalSubdomains(subpc, _localIsIndex.size(), &_overlappingIs[0], &_localIs[0]);
    }
    PCASMSetOverlap(subpc, _overlap);
    
    KSPSetUp(subksp);

    KSP* subksps;
    PCASMGetSubKSP(subpc, &_nlocal, PETSC_NULL, &subksps);
    PetscReal epsilon = 1.e-16;

    if(!_standardASM) {
      for(int i = 0; i < _blockTypeRange[1]; i++) {
        PC subpcs;
        KSPGetPC(subksps[i], &subpcs);
        KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
        KSPSetFromOptions(subksps[i]);
        PetscPreconditioner::set_petsc_preconditioner_type(MLU_PRECOND, subpcs);
        PCFactorSetZeroPivot(subpcs, epsilon);
        PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
      }

      for(int i = _blockTypeRange[1]; i < _blockTypeRange[2]; i++) {
        PC subpcs;
        KSPGetPC(subksps[i], &subpcs);
        KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
        KSPSetFromOptions(subksps[i]);

        if(this->_preconditioner_type == ILU_PRECOND)
          PCSetType(subpcs, (char*) PCILU);
        else
          PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);

        PCFactorSetZeroPivot(subpcs, epsilon);
        PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
      }
    }
    else {
      for(int i = 0; i < _nlocal; i++) {
        PC subpcs;
        KSPGetPC(subksps[i], &subpcs);
        KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
        KSPSetFromOptions(subksps[i]);

        if(this->_preconditioner_type == ILU_PRECOND)
          PCSetType(subpcs, (char*) PCILU);
        else
          PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);

        PCFactorSetZeroPivot(subpcs, epsilon);
        PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
      }
    }
  }

} //end namespace femus

#endif

