/*=========================================================================

 Program: FEMUS


} //end namespace femus


 Module: ImplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "ImplicitSystem.hpp"


namespace femus {



// ------------------------------------------------------------
// ExplicitSystem implementation
ImplicitSystem::ImplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in,const LinearEquationSolverType & smoother_type) :
  ExplicitSystem (ml_probl, name_in, number_in,smoother_type)
{
}

ImplicitSystem::~ImplicitSystem() {
}

void ImplicitSystem::init() {
  
}


} //end namespace femus






