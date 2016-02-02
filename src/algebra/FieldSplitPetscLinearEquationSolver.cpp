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
    _fieldSplitTree->BuildIndexSet(KKoffset, _iproc, _nprocs, _msh->GetLevel());
    GmresPetscLinearEquationSolver::BuildBdcIndex(variable_to_be_solved);
  }

  void FieldSplitPetscLinearEquationSolver::SetPreconditioner(KSP& subksp, PC& subpc) {
    unsigned n_split;
    _fieldSplitTree->SetPC(subksp, _msh->GetLevel());
  }

} //end namespace femus

#endif

