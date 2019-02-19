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


#include "FemusConfig.hpp"

#ifdef HAVE_PETSC


#include "FieldSplitPetscLinearEquationSolver.hpp"



namespace femus {

  void FieldSplitPetscLinearEquationSolver::SetFieldSplitTree(FieldSplitTree* fieldSplitTree) {
    _fieldSplitTree = fieldSplitTree;
  }

  void FieldSplitPetscLinearEquationSolver::BuildBdcIndex(const vector <unsigned>& variable_to_be_solved) {
    std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAA\n"<<std::flush;
    _fieldSplitTree->BuildIndexSet(KKoffset, _iproc, _nprocs, _msh->GetLevel(), this);
    std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAA\n"<<std::flush;
    GmresPetscLinearEquationSolver::BuildBdcIndex(variable_to_be_solved);
  }

  void FieldSplitPetscLinearEquationSolver::SetPreconditioner(KSP& subksp, PC& subpc) {
    
    std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAA\n"<<std::flush;
    _fieldSplitTree->SetPC(subksp, _msh->GetLevel());
  }

} //end namespace femus

#endif

