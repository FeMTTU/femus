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
#include "LinearEquationSolver.hpp"

// ------------------------------------------------------------
// NonLinearImplicitSystem implementation
NonLinearImplicitSystem::NonLinearImplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in) :
  LinearImplicitSystem (ml_probl, name_in, number_in),
  _n_nonlinear_iterations   (0),
  _n_max_nonlinear_iterations (15),
  _final_nonlinear_residual (1.e20)
{
}

NonLinearImplicitSystem::~NonLinearImplicitSystem() {
   this->clear(); 
}

void NonLinearImplicitSystem::clear() {
}

void NonLinearImplicitSystem::init() {
  Parent::init();
}

void NonLinearImplicitSystem::solve() {
  Parent::solve();
}





