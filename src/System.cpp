/*=========================================================================

 Program: FEMUS
 Module: System
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "System.hpp"

  /** Constructor.  Optionally initializes required data structures. */
  System::System (NonLinearMultiLevelProblem& es, const std::string& name_in, const unsigned int number_in) :
  _equation_systems                 (es),
//  _mesh                           (es.get_mesh()),
  _sys_name                         (name_in),
  _sys_number                       (number_in)
  {   
  }
  
  
System::~System() {
  this->clear();
}

void System::clear() {
  
}

void System::init() {
  
}

void System::AttachAssembleFunction(void fptr(NonLinearMultiLevelProblem &nl_ml_prob, unsigned level, 
				      const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix))
{
  assert(fptr);

  _assemble_system_function = fptr;
}
  
// *******************************************************
void System::AddSolutionToSytemPDE(const char solname[]){
  unsigned jsol=0;
  for(unsigned j=0;j<_SolSystemPdeIndex.size();j++){
    if(strcmp(_equation_systems.SolName[_SolSystemPdeIndex[j]],solname)) jsol++;
  }
  if(jsol==_SolSystemPdeIndex.size()){
    _SolSystemPdeIndex.resize(jsol+1u);
    _SolSystemPdeIndex[jsol]=_equation_systems.GetIndex(solname);
  }
}