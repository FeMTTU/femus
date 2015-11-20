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
#include "iomanip"

namespace femus {



// ------------------------------------------------------------
// NonLinearImplicitSystem implementation
NonLinearImplicitSystem::NonLinearImplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in, const MgSmoother & smoother_type) :
  LinearImplicitSystem (ml_probl, name_in, number_in, smoother_type),
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
    std::cout<< std::endl<<" *** Start MultiLevel Full-Cycle ***" << std::endl;
    full_cycle=1;
    igrid0=1;
  }
  else if(_mg_type == V_CYCLE){
    std::cout<< std::endl<<" *** Start MultiLevel V-Cycle ***" << std::endl;
    full_cycle=0;
    igrid0=_gridn;
  }
  else {
    std::cout<< std::endl<<" *** Start MultiLevel AMR-Cycle ***" << std::endl;
    full_cycle=0;
    igrid0=_gridr;
  }

  unsigned AMR_counter=0;

  for ( unsigned igridn=igrid0; igridn <= _gridn; igridn++) {   //_igridn

    std::cout << std::endl << " ****** Start Level Max " << igridn << " ******" << std::endl;
    clock_t start_nl_time = clock();

    bool ThisIsAMR = (_mg_type == F_CYCLE && _AMRtest &&  AMR_counter<_maxAMRlevels && igridn==_gridn)?1:0;
    if(ThisIsAMR) _solution[igridn-1]->InitAMREps();


    for ( unsigned nonLinearIterator = 0; nonLinearIterator < _n_max_nonlinear_iterations; nonLinearIterator++ ) { //non linear cycle
      std::cout << std::endl << " ********* Nonlinear iteration " << nonLinearIterator + 1 << " *********" << std::endl;

      Vcycle(igridn, full_cycle, nonLinearIterator );

      // ============== Test for non-linear Convergence ==============
      bool isnonlinearconverged = IsNonLinearConverged(igridn-1);
      if (isnonlinearconverged)
	nonLinearIterator = _n_max_nonlinear_iterations+1;
    }

    if(ThisIsAMR){
      bool conv_test=0;
      if(_AMRnorm==0){
	conv_test=_solution[_gridn-1]->FlagAMRRegionBasedOnl2(_SolSystemPdeIndex,_AMRthreshold);
      }
      else if (_AMRnorm==1){
	conv_test=_solution[_gridn-1]->FlagAMRRegionBasedOnSemiNorm(_SolSystemPdeIndex,_AMRthreshold);
      }
      if(conv_test==0){
	_ml_msh->AddAMRMeshLevel();
	_ml_sol->AddSolutionLevel();
	AddSystemLevel();
	AMR_counter++;
      }
      else{
	_maxAMRlevels=AMR_counter;
	std::cout<<"The AMR solver has converged after "<<AMR_counter<<" refinements.\n";
      }
    }

    if (igridn < _gridn) {
      ProlongatorSol(igridn);
    }

    std::cout << std::endl << " ****** Nonlinear-Cycle TIME: " << std::setw(11) << std::setprecision(6) << std::fixed
    <<static_cast<double>((clock()-start_nl_time))/CLOCKS_PER_SEC << std::endl;

    std::cout << std::endl << " ****** End Level Max "<< igridn << " ******" << std::endl;


  }

  std::cout << std::endl << " *** MultiGrid TIME: " << std::setw(11) << std::setprecision(6) << std::fixed
  <<static_cast<double>((clock()-start_mg_time))/CLOCKS_PER_SEC << std::endl;


}



bool NonLinearImplicitSystem::IsNonLinearConverged(const unsigned igridn){
  bool conv=true;
  double ResMax;
  double L2normEps;
  std::cout << std::endl;
  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    unsigned indexSol=_SolSystemPdeIndex[k];
    L2normEps    = _solution[igridn]->_Eps[indexSol]->l2_norm();
    std::cout << " ********* Level Max " << igridn+1<< " Nonlinear Eps L2norm" << std::scientific << _ml_sol->GetSolutionName(indexSol) << " = " << L2normEps <<std::endl;
    if (L2normEps <_max_nonlinear_convergence_tolerance && conv==true) {
      conv=true;
    }
    else {
      conv=false;
    }
  }
  return conv;
}

void NonLinearImplicitSystem::MGsolve (const MgSmootherType& mgSmootherType){

  clock_t start_mg_time = clock();

  unsigned igrid0;

  if(_mg_type == F_CYCLE) {
    std::cout<< std::endl<<" *** Start Multigrid Full-Cycle ***" << std::endl;
    igrid0=1;
  }
  else if(_mg_type == V_CYCLE){
    std::cout<< std::endl<<" *** Start Multigrid V-Cycle ***" << std::endl;
    igrid0=_gridn;
  }
  else {
    std::cout<<"AMR-cycle not yet implemented in MGsolve"<<std::endl;
    abort();
  }

  for ( unsigned igridn = igrid0; igridn <= _gridn; igridn++) {   //_igridn
    std::cout << std::endl << " ****** Start Level Max " << igridn << " ******" << std::endl;
    clock_t start_nl_time = clock();
    for ( unsigned nonLinearIterator = 0; nonLinearIterator < _n_max_nonlinear_iterations; nonLinearIterator++ ) {
      std::cout << std::endl << " ********* Nonlinear iteration " << nonLinearIterator + 1 << " *********" << std::endl;

      _MGmatrixFineReuse = (0 == nonLinearIterator ) ? false : true;
      _MGmatrixCoarseReuse = ( igridn - igrid0 > 0 )? _MGmatrixFineReuse + true : _MGmatrixFineReuse + false;

      MGVcycle (igridn, mgSmootherType);

      bool isnonlinearconverged = IsNonLinearConverged(igridn-1);
      if (isnonlinearconverged)
        nonLinearIterator = _n_max_nonlinear_iterations+1;
    }

    if (igridn < _gridn) {
      ProlongatorSol(igridn);
    }
    std::cout << std::endl << " ****** Nonlinear-Cycle TIME: " << std::setw(11) << std::setprecision(6) << std::fixed
    <<static_cast<double>((clock()-start_nl_time))/CLOCKS_PER_SEC << std::endl;

    std::cout << std::endl << " ****** End Level Max "<< igridn << " ******" << std::endl;
  }
  std::cout << std::endl << " *** MultiGrid TIME: " << std::setw(11) << std::setprecision(6) << std::fixed
  <<static_cast<double>((clock()-start_mg_time))/CLOCKS_PER_SEC << std::endl;

}



} //end namespace femus

