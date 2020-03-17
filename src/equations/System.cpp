/*=========================================================================

 Program: FEMUS
 Module: System
 Authors: Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "System.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "Assemble_unknown.hpp"

#include <sstream>

namespace femus {



  /** Constructor.  Optionally initializes required data structures. */
  System::System (MultiLevelProblem& ml_probl, const std::string& name_in, const unsigned int number_in, const LinearEquationSolverType & smoother_type) :
  _equation_systems                 (ml_probl),
  _sys_name                         (name_in),
  _sys_number                       (number_in),
  _gridn(ml_probl.GetNumberOfLevels()),
  _ml_sol(ml_probl._ml_sol),
  _ml_msh(ml_probl._ml_msh),
  _buildSolver( true )
{
  _msh.resize(_gridn);
  _solution.resize(_gridn);
  for(unsigned i=0;i<_gridn;i++){
    _msh[i]=_equation_systems._ml_msh->GetLevel(i);
    _solution[i]=_ml_sol->GetSolutionLevel(i);
  }
}


System::~System() {
}

void System::init() {

}

   System::AssembleFunctionType  System::GetAssembleFunction() {
  return _assemble_system_function;
}

void System::SetAssembleFunction(void fptr(MultiLevelProblem &ml_prob))
{
  assert(fptr);

  _assemble_system_function = fptr;
}

void System::AddSolutionToSystemPDEVector(const unsigned n_components, const std::string name) {

    for (unsigned i=0; i<n_components; i++) {
      std::ostringstream name_cmp; name_cmp << name << i;
       AddSolutionToSystemPDE(name_cmp.str().c_str());
     }

}


void System::AddSolutionToSystemPDE(const char solname[]){
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
      std::cout<<"error! invalid name entry" << std::endl;
      abort();
    }
  }
  return index;
}


const unsigned System::GetSolPdeIndex(const char solname[]) const {
  //unsigned ipde=GetPdeIndex(pdename);
  unsigned index=0;
  while (strcmp(_ml_sol->GetSolutionName(_SolSystemPdeIndex[index]),solname)) {
    index++;
    if (index==_SolSystemPdeIndex.size()) {
      std::cout<<"error! invalid name entry" << std::endl;
      abort();
    }
  }
  return index;
}


 void System::set_unknown_list_for_assembly(const std::vector< Unknown > unknown_in ) {
    _unknown_list_for_assembly = unknown_in; 
 }

 const std::vector< Unknown > System::get_unknown_list_for_assembly() const {
        return  _unknown_list_for_assembly; 
    }
    
 void System::assemble_call(const unsigned int n_times) const {
     
   for (unsigned it = 0; it < n_times; it++) {

            _assemble_system_function (_equation_systems);
   }
   
 }
    


} //end namespace femus


