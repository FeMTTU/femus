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
  _npre = 1;
  _npost = 1;
  _VankaIsSet = false;
  _NSchurVar = 1;
  _Schur = false;
}

ImplicitSystem::~ImplicitSystem() {
  this->clear(); 
}

void ImplicitSystem::clear() {
  
}

void ImplicitSystem::init() {
  
}




