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
  NonLinearImplicitSystem::NonLinearImplicitSystem(MultiLevelProblem& ml_probl,
      const std::string& name_in,
      const unsigned int number_in, const MgSmoother& smoother_type) :
    LinearImplicitSystem(ml_probl, name_in, number_in, smoother_type),
    _n_max_nonlinear_iterations(15),
    _final_nonlinear_residual(1.e20),
    _max_nonlinear_convergence_tolerance(1.e-6)
  {
    
  }

  NonLinearImplicitSystem::~NonLinearImplicitSystem() {
    this->clear();
  }

  void NonLinearImplicitSystem::clear() {
  }

  // ********************************************

  void NonLinearImplicitSystem::init() {
    Parent::init();
  }

  // ************************MG********************

  bool NonLinearImplicitSystem::IsNonLinearConverged(const unsigned igridn) {
    bool conv = true;
    double ResMax;
    double L2normEps;
//     std::cout << std::endl;

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned indexSol = _SolSystemPdeIndex[k];
      L2normEps    = _solution[igridn]->_Eps[indexSol]->l2_norm();
      std::cout << " Level " << igridn + 1 << " Eps L2norm" << std::scientific << _ml_sol->GetSolutionName(indexSol) << " = " << L2normEps << std::endl;

      if (L2normEps < _max_nonlinear_convergence_tolerance && conv == true) {
        conv = true;
      }
      else {
        conv = false;
      }
    }

    return conv;
  }

  // ********************************************

  void NonLinearImplicitSystem::solve(const MgSmootherType& mgSmootherType) {

    clock_t start_mg_time = clock();

    unsigned grid0;

    if (_mg_type == F_CYCLE) {
      std::cout << std::endl << " *** Start Nonlinear " << _solverType << " Full-Cycle ***" << std::endl;
      grid0 = 1;
    }
    else if (_mg_type == V_CYCLE) {
      std::cout << std::endl << " *** Start Nonlinear  " << _solverType << " V-Cycle ***" << std::endl;
      grid0 = _gridn;
    }
    else if (_mg_type == M_CYCLE) {
      std::cout << std::endl << " *** Start Nonlinear  " << _solverType << " Mixed-Cycle ***" << std::endl;
      grid0 = _gridr;
    }
    else {
      std::cout << "wrong mg_type for this solver " << std::endl;
      abort();
    }

    unsigned AMRCounter = 0;

    for (unsigned igridn = grid0; igridn <= _gridn; igridn++) {    //_igridn
      std::cout << std::endl << " ****** Start Level Max " << igridn << " ******" << std::endl;
      clock_t start_nl_time = clock();

      bool ThisIsAMR = (_mg_type == F_CYCLE && _AMRtest &&  AMRCounter < _maxAMRlevels && igridn == _gridn) ? 1 : 0;

      if (ThisIsAMR) _solution[igridn - 1]->InitAMREps();

      for (unsigned nonLinearIterator = 0; nonLinearIterator < _n_max_nonlinear_iterations; nonLinearIterator++) {
        std::cout << std::endl << " ********* Nonlinear iteration " << nonLinearIterator + 1 << " *********" << std::endl;

        _MGmatrixFineReuse = (0 == nonLinearIterator) ? false : true;
        _MGmatrixCoarseReuse = (igridn - grid0 > 0) ?  true : _MGmatrixFineReuse;

        if (_MGsolver) MGVcycle(igridn, mgSmootherType);
        else MLVcycle(igridn);

        bool nonLinearIsConverged = IsNonLinearConverged(igridn - 1);

        if (nonLinearIsConverged) break;
      }

      if (ThisIsAMR) AddAMRLevel(AMRCounter);

      if (igridn < _gridn) ProlongatorSol(igridn);

      std::cout << std::endl << " ****** Nonlinear-Cycle TIME: " << std::setw(11) << std::setprecision(6) << std::fixed
                << static_cast<double>((clock() - start_nl_time)) / CLOCKS_PER_SEC << std::endl;

      std::cout << std::endl << " ****** End Level Max " << igridn << " ******" << std::endl;
    }

    std::cout << std::endl << " *** Nonlinear "<<_solverType<<" TIME: " << std::setw(11) << std::setprecision(6) << std::fixed
              <<static_cast<double>((clock()-start_mg_time))/CLOCKS_PER_SEC << std::endl;

  }



} //end namespace femus

