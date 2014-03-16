/*=========================================================================

 Program: FEMUS
 Module: NonLinearImplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "NonLinearImplicitSystem.hpp"
#include "LinearSolver.hpp"

// ------------------------------------------------------------
// NonLinearImplicitSystem implementation
NonLinearImplicitSystem::NonLinearImplicitSystem (NonLinearMultiLevelProblem& es,
				const std::string& name_in,
				const unsigned int number_in) :
  ImplicitSystem (es, name_in, number_in),
  _n_nonlinear_iterations   (0),
  _final_nonlinear_residual (1.e20)
{
}

NonLinearImplicitSystem::~NonLinearImplicitSystem() {
   this->clear(); 
}

void NonLinearImplicitSystem::clear() {
  for (unsigned ig=0; ig<_equation_systems.GetNumberOfGrid(); ig++) {
    _LinSolver[ig]->DeletePde();
    delete _LinSolver[ig];
  }
}

void NonLinearImplicitSystem::init() {
  CreateSystemPDEStructure();
}

void NonLinearImplicitSystem::CreateSystemPDEStructure() {
    _LinSolver.resize(_equation_systems.GetNumberOfGrid());
    for(unsigned i=0;i<_equation_systems.GetNumberOfGrid();i++){
      _LinSolver[i]=LinearSolver::build(i,_equation_systems._msh[i]).release();
    }
    
    for (unsigned i=0; i<_equation_systems.GetNumberOfGrid(); i++) {
      //_LinSolver[i]->InitPde(_SolPdeIndex[ipde],SolType,SolName,&_solution[i]->_Bdc,_gridr,_gridn);
      _LinSolver[i]->InitPde(_SolSystemPdeIndex,_equation_systems.SolType,_equation_systems.SolName,&_equation_systems._solution[i]->_Bdc,_equation_systems.GetNumberOfGridNotRefined(),_equation_systems.GetNumberOfGrid());
    }  
}





