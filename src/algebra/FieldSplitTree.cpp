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

#include "FieldSplitTree.hpp"
#include "FieldSplitPetscLinearEquationSolver.hpp"

namespace femus {

  FieldSplitTree::FieldSplitTree( const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, const std::vector < unsigned > & solutionType, std::string name ){
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
       
    if( preconditioner == ASM_PRECOND && _solutionType.size() != fields.size() ){
      std::cout << "Error! The Solution type with PCFieldSplit - ASM preconditioner has to be specified"<<std::endl;
      abort();
    }
    
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
  }


  FieldSplitTree::~FieldSplitTree() {
    if(_numberOfSplits > 1) {
      for(unsigned i = 0; i < _isSplit.size(); i++) {
        for(unsigned j = 0; j < _isSplit[i].size(); j++) {
          ISDestroy(&_isSplit[i][j]);
        }
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
      _MatrixOffset[level-1] = KKoffset;
      
    if ( GetNumberOfSplits() == 1) {
      if ( _preconditioner == ASM_PRECOND ){
	BuildASMIndexSet( level, solver );
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

    std::cout<<level<<std::endl;
    for(int i = 0; i < _MatrixOffset[level-1].size(); i++){
      for(int j = 0; j < _MatrixOffset[level-1][i].size(); j++){
	std::cout << _MatrixOffset[level-1][i][j] << " ";
      }
      std::cout <<std::endl;
    }
    std::cout <<std::endl;
    
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
      
      for(int i = 0; i < _solutionType.size(); i++){
	std::cout<<"solution Type " << i << " = " << _solutionType[i]<<std::endl;
      }
      
//       _standardASM = false;
//       
//       
//       
//       PetscPreconditioner::set_petsc_preconditioner_type(ASM_PRECOND, pc);
// 
//       if(!_standardASM) {
// 	PCASMSetLocalSubdomains(pc, _localIsIndex[level-1].size(), &_overlappingIs[level-1][0], &_localIs[level-1][0]);
//       }
// 
//       PCASMSetOverlap(pc, _overlap);
//       
// 
//       KSPSetUp(ksp);
// 
//       KSP* subksps;
//       PCASMGetSubKSP(pc, &_nlocal[level-1], PETSC_NULL, &subksps);
//       PetscReal epsilon = 1.e-16;
// 
//       if(!_standardASM) {
// 	for(int i = 0; i < _blockTypeRange[level-1][0]; i++) {
// 	  PC subpcs;
// 	  KSPGetPC(subksps[i], &subpcs);
// 	  KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
// 	  KSPSetFromOptions(subksps[i]);
// 	  PetscPreconditioner::set_petsc_preconditioner_type(MLU_PRECOND, subpcs);
// 	  PCFactorSetZeroPivot(subpcs, epsilon);
// 	  PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
// 	}
// 
// 	for(int i = _blockTypeRange[level-1][0]; i < _blockTypeRange[level-1][1]; i++) {
// 	  PC subpcs;
// 	  KSPGetPC(subksps[i], &subpcs);
// 	  KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
// 	  KSPSetFromOptions(subksps[i]);
// 
// 	  //if(this->_preconditioner_type == ILU_PRECOND)
// 	    PCSetType(subpcs, (char*) PCILU);
// 	  //else
// 	  //  PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);
// 
// 	  PCFactorSetZeroPivot(subpcs, epsilon);
// 	  PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
// 	}
//       }
//       else {
// 	for(int i = 0; i < _nlocal[level-1]; i++) {
// 	  PC subpcs;
// 	  KSPGetPC(subksps[i], &subpcs);
// 	  KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
// 	  KSPSetFromOptions(subksps[i]);
// 
// 	  //if(this->_preconditioner_type == ILU_PRECOND)
// 	    PCSetType(subpcs, (char*) PCILU);
// 	  //else
// 	  //  PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);
// 
// 	  PCFactorSetZeroPivot(subpcs, epsilon);
// 	  PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
// 	}
//       }
//     
      

      
      PetscPreconditioner::set_petsc_preconditioner_type(_preconditioner, pc);

      bool _standardASM = 1;
      PetscInt _nlocal;

      if(!_standardASM) {
        //PCASMSetLocalSubdomains(subpc, _localIsIndex.size(), &_overlappingIs[0], &_localIs[0]);
      }

      PCASMSetOverlap(pc, 0); //PCASMSetOverlap(subpc, _overlap);

      KSPSetUp(ksp);

      KSP* subksps;
      PCASMGetSubKSP(pc, &_nlocal, PETSC_NULL, &subksps);
      PetscReal epsilon = 1.e-16;

      if(!_standardASM) {
// 	for(int i = 0; i < _blockTypeRange[0]; i++) {
// 	  PC subpcs;
// 	  KSPGetPC(subksps[i], &subpcs);
// 	  KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
// 	  KSPSetFromOptions(subksps[i]);
// 	  PetscPreconditioner::set_petsc_preconditioner_type(MLU_PRECOND, subpcs);
// 	  PCFactorSetZeroPivot(subpcs, epsilon);
// 	  PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
// 	}
//
// 	for(int i = _blockTypeRange[0]; i < _blockTypeRange[1]; i++) {
// 	  PC subpcs;
// 	  KSPGetPC(subksps[i], &subpcs);
// 	  KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
// 	  KSPSetFromOptions(subksps[i]);
//
// 	  if(this->_preconditioner_type == ILU_PRECOND)
// 	    PCSetType(subpcs, (char*) PCILU);
// 	  else
// 	    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);
//
// 	  PCFactorSetZeroPivot(subpcs, epsilon);
// 	  PCFactorSetShiftType(subpcs, MAT_SHIFT_NONZERO);
// 	}
      }
      else {
        for(int i = 0; i < _nlocal; i++) {
          PC subpcs;
          KSPGetPC(subksps[i], &subpcs);
          KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
          KSPSetFromOptions(subksps[i]);


          PCSetType(subpcs, (char*) PCILU);

// 	  if(this->_preconditioner_type == ILU_PRECOND)
// 	    PCSetType(subpcs, (char*) PCILU);
// 	  else
// 	    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type, subpcs);

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
      SetPetscSolverType(ksp);
      PC pc;
      KSPGetPC(ksp, &pc);
      KSPSetTolerances(ksp, _rtol, _abstol, _dtol, _maxits);
      KSPSetFromOptions(ksp);
      PetscReal epsilon = 1.e-16;

      /* adjusted by Guoyi Ke*/
      if(_preconditioner == LSC_PRECOND) {
        if(_solver == PREONLY) {
          std::cout << "LSC does not allow to use PREONLY, and ksp switches to GMRES instead!" << std::endl;
          _solver = GMRES;
        }
      }
      /* adjusted by Guoyi Ke*/

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
  

  void FieldSplitTree::BuildASMIndexSet( const unsigned& level, const FieldSplitPetscLinearEquationSolver *solver){
    
    Mesh* msh = solver->_msh;
    
    unsigned nel = msh->GetNumberOfElements();

    bool FastVankaBlock = true;

    
    unsigned _NSchurVar = 1;
    
    if(_NSchurVar != 0) {
      FastVankaBlock = (_solutionType[ _solutionType.size() - _NSchurVar ] < 3 ) ? false : true;
    }

    unsigned iproc = solver->processor_id();

    unsigned DofOffset = _MatrixOffset[level-1][0][iproc];
    unsigned DofOffsetSize = _MatrixOffset[level-1][_solutionType.size()][iproc] - DofOffset;
        
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

//     MeshASMPartitioning meshasmpartitioning(*_msh);
// 
//     meshasmpartitioning.DoPartition(_elementBlockNumber, block_elements, _blockTypeRange);
// 
//     vector <bool> ThisVaribaleIsNonSchur(_SolPdeIndex.size(), true);
// 
//     for(unsigned iind = variable_to_be_solved.size() - _NSchurVar; iind < variable_to_be_solved.size(); iind++) {
//       unsigned PdeIndexSol = variable_to_be_solved[iind];
//       ThisVaribaleIsNonSchur[PdeIndexSol] = false;
//     }
// 
//     // *** Start Vanka Block ***
// 
//     _localIsIndex.resize(block_elements.size());
//     _overlappingIsIndex.resize(block_elements.size());
// 
//     for(int vb_index = 0; vb_index < block_elements.size(); vb_index++) { //loop on the vanka-blocks
//       _localIsIndex[vb_index].resize(DofOffsetSize);
//       _overlappingIsIndex[vb_index].resize(DofOffsetSize);
// 
//       PetscInt PAsize = 0;
//       PetscInt PBsize = 0;
// 
//       PetscInt Csize = 0;
// 
//       // ***************** NODE/ELEMENT SERCH *******************
//       for(int kel = 0; kel < block_elements[vb_index].size(); kel++) { //loop on the vanka-block elements
//         unsigned iel = block_elements[vb_index][kel];
// 
//         for(unsigned i = 0; i < _msh->GetElementDofNumber(iel, 0); i++) { //loop on the element vertices
//           unsigned inode = _msh->el->GetElementDofIndex(iel, i);
//           const std::vector < unsigned > & localElementNearVertexNumber = _msh->el->GetLocalElementNearVertex(inode);
//           //loop on the neighboring elemnets (!FastVankaBlock) or on iel only (FastVankaBlock) 
// 	  unsigned nve = (FastVankaBlock) ? 1 : localElementNearVertexNumber.size();
//           for(unsigned j = 0; j < nve; j++) { 
//             unsigned jel = (!FastVankaBlock) ? localElementNearVertexNumber[j] : iel;
// 
//             //add elements for velocity to be solved
//             if(indexc[jel - ElemOffset] == ElemOffsetSize) {
//               indexci[Csize] = jel - ElemOffset;
//               indexc[jel - ElemOffset] = Csize++;
// 
//               //add non-schur variables to be solved
//               for(int indexSol = 0; indexSol < _SolPdeIndex.size(); indexSol++) {
//                 if(ThisVaribaleIsNonSchur[indexSol]) {
//                   unsigned SolPdeIndex = _SolPdeIndex[indexSol];
//                   unsigned SolType = _SolType[SolPdeIndex];
//                   unsigned nvej = _msh->GetElementDofNumber(jel, SolType);
// 
//                   for(unsigned jj = 0; jj < nvej; jj++) {
//                     unsigned jdof = _msh->GetSolutionDof(jj, jel, SolType);
//                     unsigned kkdof = GetSystemDof(SolPdeIndex, indexSol, jj, jel);
// 
//                     if(jdof >= _msh->_dofOffset[SolType][iproc] &&
//                         jdof <  _msh->_dofOffset[SolType][iproc + 1]) {
//                       if(indexa[kkdof - DofOffset] == DofOffsetSize && owned[kkdof - DofOffset] == false) {
//                         owned[kkdof - DofOffset] = true;
//                         _localIsIndex[vb_index][PAsize] = kkdof;
//                         indexa[kkdof - DofOffset] = PAsize++;
//                       }
// 
//                       if(indexb[kkdof - DofOffset] == DofOffsetSize) {
//                         _overlappingIsIndex[vb_index][PBsize] = kkdof;
//                         indexb[kkdof - DofOffset] = PBsize++;
//                       }
//                     }
//                     else {
//                       mymap[kkdof] = true;
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
// 
//         //-----------------------------------------------------------------------------------------
//         //Add Schur nodes (generally pressure type variables) to be solved
//         {
//           for(int indexSol = 0; indexSol < _SolPdeIndex.size(); indexSol++) {
//             if(!ThisVaribaleIsNonSchur[indexSol]) {
//               unsigned SolPdeIndex = _SolPdeIndex[indexSol];
//               unsigned SolType = _SolType[SolPdeIndex];
//               unsigned nvei = _msh->GetElementDofNumber(iel, SolType);
// 
//               for(unsigned ii = 0; ii < nvei; ii++) {
//                 unsigned inode_Metis = _msh->GetSolutionDof(ii, iel, SolType);
//                 unsigned kkdof = GetSystemDof(SolPdeIndex, indexSol, ii, iel);
// 
//                 if(inode_Metis >= _msh->_dofOffset[SolType][iproc] &&
//                     inode_Metis <  _msh->_dofOffset[SolType][iproc + 1]) {
//                   if(indexa[kkdof - DofOffset] == DofOffsetSize && owned[kkdof - DofOffset] == false) {
//                     owned[kkdof - DofOffset] = true;
//                     _localIsIndex[vb_index][PAsize] = kkdof;
//                     indexa[kkdof - DofOffset] = PAsize++;
//                   }
// 
//                   if(indexb[kkdof - DofOffset] == DofOffsetSize) {
//                     _overlappingIsIndex[vb_index][PBsize] = kkdof;
//                     indexb[kkdof - DofOffset] = PBsize++;
//                   }
//                 }
//                 else {
//                   mymap[kkdof] = true;
//                 }
//               }
//             }
//           }
//         }
//         //-----------------------------------------------------------------------------------------
//       }
// 
//       // *** re-initialize indeces(a,c,d)
//       for(PetscInt i = 0; i < PAsize; i++) {
//         indexa[_localIsIndex[vb_index][i] - DofOffset] = DofOffsetSize;
//       }
// 
//       for(PetscInt i = 0; i < PBsize; i++) {
//         indexb[_overlappingIsIndex[vb_index][i] - DofOffset] = DofOffsetSize;
//       }
// 
//       for(PetscInt i = 0; i < Csize; i++) {
//         indexc[indexci[i]] = ElemOffsetSize;
//       }
// 
//       _localIsIndex[vb_index].resize(PAsize);
//       std::vector < PetscInt >(_localIsIndex[vb_index]).swap(_localIsIndex[vb_index]);
// 
//       _overlappingIsIndex[vb_index].resize(PBsize + mymap.size());
//       int i = 0;
//       for(std::map<int, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it, ++i) {
//         _overlappingIsIndex[vb_index][PBsize + i] = it->first;
//       }
//       std::vector < PetscInt >(_overlappingIsIndex[vb_index]).swap(_overlappingIsIndex[vb_index]);
//       
//      
//       mymap.clear();
// 
//       std::sort(_localIsIndex[vb_index].begin(), _localIsIndex[vb_index].end());
//       std::sort(_overlappingIsIndex[vb_index].begin(), _overlappingIsIndex[vb_index].end());
// 
// 
//     }
// 
//     //BEGIN Generate std::vector<IS> for ASM PC ***********
//     _localIs.resize(_localIsIndex.size());
//     _overlappingIs.resize(_overlappingIsIndex.size());
// 
//     for(unsigned vb_index = 0; vb_index < _localIsIndex.size(); vb_index++) {
//       ISCreateGeneral(MPI_COMM_SELF, _localIsIndex[vb_index].size(), &_localIsIndex[vb_index][0], PETSC_USE_POINTER, &_localIs[vb_index]);
//       ISCreateGeneral(MPI_COMM_SELF, _overlappingIsIndex[vb_index].size(), &_overlappingIsIndex[vb_index][0], PETSC_USE_POINTER, &_overlappingIs[vb_index]);
//     }
// 
//     //END Generate std::vector<IS> for ASM PC ***********

    return;

//    
    
    
    
  }
  
}

