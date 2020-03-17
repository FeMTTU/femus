/*=========================================================================

  Program: FEMUS
  Module: LinearEquationSolverPetscFieldSplit
  Authors: Eugenio Aulisa, Guoyi Ke

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/


#include "FemusConfig.hpp"

#ifdef HAVE_PETSC


#include "LinearEquationSolverPetscFieldSplit.hpp"



namespace femus {

  void LinearEquationSolverPetscFieldSplit::SetFieldSplitTree(FieldSplitTree* fieldSplitTree) {
    _fieldSplitTree = fieldSplitTree;
  }

  void LinearEquationSolverPetscFieldSplit::BuildBdcIndex(const vector <unsigned>& variable_to_be_solved) {
    if(_fieldSplitTree != NULL) _fieldSplitTree->BuildIndexSet(KKoffset, _iproc, _nprocs, _msh->GetLevel(), this);
    else FielSlipTreeIsNotDefined();
    LinearEquationSolverPetsc::BuildBdcIndex(variable_to_be_solved);
  }

  void LinearEquationSolverPetscFieldSplit::SetPreconditioner(KSP& subksp, PC& subpc) {
    if(_fieldSplitTree != NULL) _fieldSplitTree->SetPC(subksp, _msh->GetLevel());
    else FielSlipTreeIsNotDefined();
  }
  
  void LinearEquationSolverPetscFieldSplit::FielSlipTreeIsNotDefined(){
    std::cout << "Error! No FieldSplitTree object has been passed to the FEMuS_FIELDSPLIT system"<<std::endl;
    std::cout << "Define a FieldSplitTree object FS and pass it to the FEMuS_FIELDSPLIT system with"<<std::endl;
    std::cout << "system.SetFieldSplitTree(&FS);"<<std::endl;
    abort();
  }
  

} //end namespace femus

#endif

