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
#include "NumericVector.hpp"

// ------------------------------------------------------------
// NonLinearImplicitSystem implementation
NonLinearImplicitSystem::NonLinearImplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in) :
  LinearImplicitSystem (ml_probl, name_in, number_in),
  _n_nonlinear_iterations   (0),
  _n_max_nonlinear_iterations (15),
  _final_nonlinear_residual (1.e20),
  _max_nonlinear_convergence_tolerance(1.e-6)
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
  
  clock_t start_mg_time = clock();
  
  bool full_cycle;
  unsigned igrid0; 
  
  if(_mg_type == F_CYCLE) {
    full_cycle=1;
    igrid0=1;
  }
  else if(_mg_type == V_CYCLE){
    full_cycle=0;
    igrid0=_gridn;
  }
  else {
    full_cycle=0;
    igrid0=_gridr;
  }
  
  std::pair<int, double> solver_info;
     
  for ( unsigned igridn=igrid0; igridn <= _gridn; igridn++) {   //_igridn
    
    std::cout << std::endl << "    ************* Level Max: " << igridn << " *************\n" << std::endl;

     
    for ( _n_nonlinear_iterations = 0; _n_nonlinear_iterations < _n_max_nonlinear_iterations; _n_nonlinear_iterations++ ) { //non linear cycle
      clock_t start_time = clock();
      std::cout << std::endl << "    ************** NonLinear Cycle: " << _n_nonlinear_iterations + 1 << " ****************\n" << std::endl;
      // ============== Fine level Assembly ==============
      _LinSolver[igridn-1u]->SetResZero();
      _LinSolver[igridn-1u]->SetEpsZero();
      bool assemble_matrix = true; //Be carefull!!!! this is needed in the _assemble_function
      
      /// Be careful !!!! adesso stiamo usando _sys_number invece che ipde, da togliere al + presto
      _assemble_system_function(_equation_systems, igridn-1u, igridn-1u, assemble_matrix);    
      
      std::cout << "Grid: " << igridn-1 << "\t        ASSEMBLY TIME:\t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
 
      for(_n_linear_iterations = 0; _n_linear_iterations < _n_max_linear_iterations; _n_linear_iterations++) { //linear cycle
	
	bool ksp_clean=!_n_linear_iterations;
	
	for (unsigned ig = igridn-1u; ig > 0; ig--) {
	  
	  // ============== Presmoothing ============== 
	  for (unsigned k = 0; k < _npre; k++) {
	    solver_info = (_VankaIsSet) ? _LinSolver[ig]->solve(_VankaIndex, _NSchurVar, _Schur) : _LinSolver[ig]->solve(ksp_clean*(!k));
	  }
	  // ============== Non-Standard Multigrid Restriction ==============
	  start_time = clock();
	  Restrictor(ig, igridn, _n_nonlinear_iterations, _n_linear_iterations, full_cycle);
	  std::cout << "Grid: " << ig << "-->" << ig-1 << "  RESTRICTION TIME:\t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
	}
       
 	// ============== Coarse Direct Solver ==============
 	solver_info = ( _VankaIsSet ) ? _LinSolver[0]->solve(_VankaIndex, _NSchurVar, _Schur) : _LinSolver[0]->solve(ksp_clean);
 	
             
 	for (unsigned ig = 1; ig < igridn; ig++) {
 	  
 	  // ============== Standard Prolongation ==============
 	  start_time=clock();
 	  Prolongator(ig);
 	  std::cout << "Grid: " << ig-1 << "-->" << ig << " PROLUNGATION TIME:\t" << static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
 
 	  // ============== PostSmoothing ==============    
 	  for (unsigned k = 0; k < _npost; k++) {
 	    solver_info = ( _VankaIsSet ) ? _LinSolver[ig]->solve(_VankaIndex, _NSchurVar, _Schur) : _LinSolver[ig]->solve(ksp_clean*(!_npre)*(!k));
 	  }
 	}
 	// ============== Update Solution ( _gridr-1 <= ig <= igridn-2 ) ==============
 	for (unsigned ig = _gridr-1; ig < igridn-1; ig++) {  // _gridr
 	  _solution[ig]->SumEpsToSol(_SolSystemPdeIndex, _LinSolver[ig]->_EPS, _LinSolver[ig]->_RES, _LinSolver[ig]->KKoffset );	
 	}
 	
 	_final_linear_residual = solver_info.second;
	// ============== Test for linear Convergence (now we are using only the absolute convergence tolerance)==============
 	if(_final_linear_residual < _absolute_convergence_tolerance) 
	  break;
      }
      std::cout <<"GRID: "<<igridn-1<< "\t    FINAL LINEAR RESIDUAL:\t"<< _final_linear_residual << std::endl;
      // ============== Update Solution ( ig = igridn )==============
      _solution[igridn-1]->SumEpsToSol(_SolSystemPdeIndex, _LinSolver[igridn-1]->_EPS, 
							 _LinSolver[igridn-1]->_RES, _LinSolver[igridn-1]->KKoffset );
      // ============== Test for non-linear Convergence ==============
      bool conv = CheckConvergence(_sys_name.c_str(), igridn-1);
      if (conv == true) _n_nonlinear_iterations = _n_max_nonlinear_iterations+1;
     
      std::cout << std::endl;
      std::cout << "COMPUTATION RESIDUAL: \t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
    }
    // ==============  Solution Prolongation ==============
    if (igridn < _gridn) {
      ProlongatorSol(igridn);
    }
  }

  std::cout << "SOLVER TIME:   \t\t\t"<<static_cast<double>((clock()-start_mg_time))/CLOCKS_PER_SEC << std::endl;
  
}



bool NonLinearImplicitSystem::CheckConvergence(const char pdename[], const unsigned igridn){

  bool conv=true;
  double ResMax;
  double L2normEps;
 
  
  //for debugging purpose
  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    unsigned indexSol=_SolSystemPdeIndex[k];
    
    L2normEps    = _solution[igridn]->_Eps[indexSol]->l2_norm();
    ResMax       = _solution[igridn]->_Res[indexSol]->linfty_norm();

    std::cout << "level=" << igridn<< "\tLinftynormRes" << _ml_sol->GetSolutionName(indexSol) << "=" << ResMax    <<std::endl;
    std::cout << "level=" << igridn<< "\tL2normEps"     << _ml_sol->GetSolutionName(indexSol) << "=" << L2normEps <<std::endl;
    
    if (L2normEps <_max_nonlinear_convergence_tolerance && conv==true) {
      conv=true;
    } 
    else {
      conv=false;
    }
  }
  return conv;
}



