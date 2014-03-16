/*=========================================================================

 Program: FEMUS
 Module: ImplicitSystem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "ImplicitSystem.hpp"

// ------------------------------------------------------------
// ExplicitSystem implementation
ImplicitSystem::ImplicitSystem (NonLinearMultiLevelProblem& es,
				const std::string& name_in,
				const unsigned int number_in) :
  ExplicitSystem (es, name_in, number_in)
{

}

ImplicitSystem::~ImplicitSystem() {
  this->clear(); 
}

void ImplicitSystem::clear() {
  
}

void ImplicitSystem::init() {
  
}