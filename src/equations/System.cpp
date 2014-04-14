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


namespace femus {



  /** Constructor.  Optionally initializes required data structures. */
  System::System (MultiLevelProblem& ml_probl, const std::string& name_in, const unsigned int number_in) :
  _equation_systems                 (ml_probl),
  _sys_name                         (name_in),
  _sys_number                       (number_in),
  _gridn(ml_probl.GetNumberOfGrid()), 
  _gridr(ml_probl.GetNumberOfGridTotallyRefined()),
  _ml_sol(ml_probl._ml_sol)
{ 
  _msh.resize(_gridn);
  _solution.resize(_gridn);
  for(unsigned i=0;i<_gridn;i++){
    _msh[i]=_equation_systems._ml_msh->GetLevel(i);
    _solution[i]=_ml_sol->GetSolutionLevel(i);
  }
}
  
  
System::~System() {
  this->clear();
}

void System::clear() {
  
}

void System::init() {
  
}

void System::AttachAssembleFunction(void fptr(MultiLevelProblem &ml_prob, unsigned level, 
				      const unsigned &gridn, const bool &assembe_matrix))
{
  assert(fptr);

  _assemble_system_function = fptr;
}
  

void System::AddSolutionToSytemPDE(const char solname[]){
  unsigned jsol=0;
  for(unsigned j=0;j<_SolSystemPdeIndex.size();j++){
    if(strcmp(_ml_sol->GetSolutionName(_SolSystemPdeIndex[j]),solname)) jsol++;
  }
  if(jsol==_SolSystemPdeIndex.size()){
    _SolSystemPdeIndex.resize(jsol+1u);
    _SolSystemPdeIndex[jsol]=_ml_sol->GetIndex(solname);
  }
}


unsigned System::GetSolPdeIndex(const char solname[]) {
  //unsigned ipde=GetPdeIndex(pdename);  
  unsigned index=0;
  while (strcmp(_ml_sol->GetSolutionName(_SolSystemPdeIndex[index]),solname)) {
    index++;
    if (index==_SolSystemPdeIndex.size()) {
      std::cout<<"error! invalid name entry MultiLevelProblem::GetSolPdeIndex(const char pdename[], const char solname[])"<<std::endl;
      exit(0);
    }
  }
  return index;
}


} //end namespace femus


