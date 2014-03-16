/*=========================================================================

 Program: FEMUS
 Module: ExplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "ExplicitSystem.hpp"

// ------------------------------------------------------------
// ExplicitSystem implementation
ExplicitSystem::ExplicitSystem (NonLinearMultiLevelProblem& es,
				const std::string& name_in,
				const unsigned int number_in) :
  System (es, name_in, number_in)
{
}

ExplicitSystem::~ExplicitSystem() {
  this->clear(); 
}
  
void ExplicitSystem::clear() {
  
}

void ExplicitSystem::init() {
  
}