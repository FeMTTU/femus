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
  _max_nonlinear_convergence_tolerance(1.e-4)
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
    igrid0=_equation_systems.GetNumberOfGrid();
  }
  else {
    full_cycle=0;
    igrid0=_equation_systems.GetNumberOfGridNotRefined();
  }
  
  unsigned ipde = _equation_systems.GetPdeIndex(_sys_name.c_str());  // non esiste + la chiamata alla pde, siamo già dentro ad una pde specifica
    
  std::pair<int, double> solver_info;
     
  for ( unsigned igridn=igrid0; igridn <= _equation_systems.GetNumberOfGrid(); igridn++) {   //_igridn
    
    std::cout << std::endl << "    ************* Level Max: " << igridn << " *************\n" << std::endl;

     
    for ( _n_nonlinear_iterations = 0; _n_nonlinear_iterations < _n_max_nonlinear_iterations; _n_nonlinear_iterations++ ) { //non linear cycle
      clock_t start_time = clock();
      std::cout << std::endl << "    ************** NonLinear Cycle: " << _n_nonlinear_iterations + 1 << " ****************\n" << std::endl;
      // ============== Fine level Assembly ==============
      _LinSolver[igridn-1u]->SetResZero();
      _LinSolver[igridn-1u]->SetEpsZero();
      bool assemble_matrix = true; //Be carefull!!!! this is needed in the _assemble_function
      
      /// Be careful !!!! adesso stiamo usando _sys_number invece che ipde, da togliere al + presto
      _assemble_system_function(_equation_systems,igridn-1u, igridn-1u,ipde, assemble_matrix);    
      
      std::cout << "Grid: " << igridn-1 << "\t        ASSEMBLY TIME:\t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
 
      for(_n_linear_iterations = 0; _n_linear_iterations < _n_max_linear_iterations; _n_linear_iterations++) { //linear cycle
 
	for (unsigned ig = igridn-1u; ig > 0; ig--) {
	  
	  // ============== Presmoothing ============== 
	  for (unsigned k = 0; k < _npre; k++) {
	    solver_info = (_VankaIsSet) ? _LinSolver[ig]->solve(_VankaIndex, _NSchurVar, _Schur) : _LinSolver[ig]->solve();
	  }
	  // ============== Non-Standard Multigrid Restriction ==============
	  start_time = clock();
	  Restrictor(ig, igridn, _n_nonlinear_iterations, _n_linear_iterations, full_cycle);
	  std::cout << "Grid: " << ig << "-->" << ig-1 << "  RESTRICTION TIME:\t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
	}
       
 	// ============== Coarse Direct Solver ==============
 	solver_info = ( _VankaIsSet ) ? _LinSolver[0]->solve(_VankaIndex, _NSchurVar, _Schur) : _LinSolver[0]->solve();
 	
             
 	for (unsigned ig = 1; ig < igridn; ig++) {
 	  
 	  // ============== Standard Prolongation ==============
 	  start_time=clock();
 	  Prolongator(ig);
 	  std::cout << "Grid: " << ig-1 << "-->" << ig << " PROLUNGATION TIME:\t" << static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
 
 	  // ============== PostSmoothing ==============    
 	  for (unsigned k = 0; k < _npost; k++) {
 	    solver_info = ( _VankaIsSet ) ? _LinSolver[ig]->solve(_VankaIndex, _NSchurVar, _Schur) : _LinSolver[ig]->solve();
 	  }
 	}
 	// ============== Update Solution ( _gridr-1 <= ig <= igridn-2 ) ==============
 	for (unsigned ig = _equation_systems.GetNumberOfGridNotRefined()-1; ig < igridn-1; ig++) {  // _gridr
 	  _equation_systems._solution[ig]->SumEpsToSol(_equation_systems._SolPdeIndex[ipde], _LinSolver[ig]->_EPS, _LinSolver[ig]->_RES, _LinSolver[ig]->KKoffset );	
 	}
 	
 	_final_linear_residual = solver_info.second;
	// ============== Test for linear Convergence (now we are using only the absolute convergence tolerance)==============
 	if(_final_linear_residual < _absolute_convergence_tolerance) 
	   break;
      }
      
      // ============== Update Solution ( ig = igridn )==============
      _equation_systems._solution[igridn-1]->SumEpsToSol(_equation_systems._SolPdeIndex[ipde], _LinSolver[igridn-1]->_EPS, 
							 _LinSolver[igridn-1]->_RES, _LinSolver[igridn-1]->KKoffset );
      // ============== Test for non-linear Convergence ==============
      bool conv = CheckConvergence(_sys_name.c_str(), igridn-1);
      if (conv == true) _n_nonlinear_iterations = _n_max_nonlinear_iterations+1;
     
      std::cout << std::endl;
      std::cout << "COMPUTATION RESIDUAL: \t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
    }
    // ==============  Solution Prolongation ==============
    if (igridn < _equation_systems.GetNumberOfGrid()) {
      _equation_systems.ProlongatorSol(_sys_name.c_str(), igridn);
    }
  }

  std::cout << "SOLVER TIME:   \t\t\t"<<static_cast<double>((clock()-start_mg_time))/CLOCKS_PER_SEC << std::endl;
  
}



bool NonLinearImplicitSystem::CheckConvergence(const char pdename[], const unsigned gridn){

  unsigned ipde=_equation_systems.GetPdeIndex(pdename);
  bool conv=true;
  double ResMax;
  double L2normEps;
  Solution *solution=_equation_systems._solution[gridn];
  
  //for debugging purpose
  for (unsigned k=0; k<_equation_systems._SolPdeIndex[ipde].size(); k++) {
    unsigned indexSol=_equation_systems._SolPdeIndex[ipde][k];
    
    L2normEps    = solution->_Eps[indexSol]->l2_norm();
    ResMax       = solution->_Res[indexSol]->linfty_norm();

    std::cout << "level=" << gridn<< "\tLinftynormRes" << _equation_systems.SolName[indexSol] << "=" << ResMax    <<std::endl;
    std::cout << "level=" << gridn<< "\tL2normEps"     << _equation_systems.SolName[indexSol] << "=" << L2normEps <<std::endl;
    
    if (L2normEps <_max_nonlinear_convergence_tolerance && conv==true) {
      conv=true;
    } 
    else {
      conv=false;
    }
  }
  return conv;
}



