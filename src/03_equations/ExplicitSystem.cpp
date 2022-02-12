/*=========================================================================

 Program: FEMUS


} //end namespace femus


 Module: ExplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "ExplicitSystem.hpp"


namespace femus {



// ------------------------------------------------------------
// ExplicitSystem implementation
ExplicitSystem::ExplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in, const LinearEquationSolverType & smoother_type) :
  System (ml_probl, name_in, number_in, smoother_type)
{
}

ExplicitSystem::~ExplicitSystem() {
}
  

void ExplicitSystem::init() {
  
}



} //end namespace femus


