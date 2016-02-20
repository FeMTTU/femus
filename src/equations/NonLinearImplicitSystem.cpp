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
  NonLinearImplicitSystem::NonLinearImplicitSystem( MultiLevelProblem& ml_probl,
      const std::string& name_in,
      const unsigned int number_in, const MgSmoother& smoother_type ) :
    LinearImplicitSystem( ml_probl, name_in, number_in, smoother_type ),
    _n_max_nonlinear_iterations( 15 ),
    _final_nonlinear_residual( 1.e20 ),
    _max_nonlinear_convergence_tolerance( 1.e-6 ),
    _maxNumberOfResidualUpdateIterations( 1 )
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

  bool NonLinearImplicitSystem::IsNonLinearConverged( const unsigned igridn, double &nonLinearEps ) {
    bool conv = true;
    double L2normEps, L2normSol, L2normEpsDividedSol;

    nonLinearEps = 0.;
    const double absMinNonlinearEps = 1.e-50;
    const double absMinNormSol = 1.e-15;
    const double mindeltaNormSol = 1.e-50;

    for( unsigned k = 0; k < _SolSystemPdeIndex.size(); k++ ) {
      unsigned indexSol = _SolSystemPdeIndex[k];
      L2normEps    = _solution[igridn]->_Eps[indexSol]->l2_norm();
      L2normSol    = _solution[igridn]->_Sol[indexSol]->l2_norm();
      L2normEpsDividedSol = L2normEps/(L2normSol+mindeltaNormSol);

      std::cout << "     ********* Level Max " << igridn + 1 << " Nonlinear Eps_l2norm/Sol_l2norm " << \
        std::scientific << _ml_sol->GetSolutionName( indexSol ) << "= " << L2normEpsDividedSol << \
        "  ** Eps_l2norm= " << L2normEps << "  ** Sol_l2norm= " << L2normSol << std::endl;
      nonLinearEps = ( nonLinearEps > L2normEpsDividedSol ) ? nonLinearEps : L2normEpsDividedSol;

      if( (L2normEpsDividedSol < _max_nonlinear_convergence_tolerance || L2normEps < absMinNonlinearEps || L2normSol < absMinNormSol) && conv == true ) {
        conv = true;
      }
      else {
        conv = false;
      }
    }

    return conv;
  }

  // ********************************************

  void NonLinearImplicitSystem::solve( const MgSmootherType& mgSmootherType ) {

    clock_t start_mg_time = clock();

    unsigned grid0;

    if( _mg_type == F_CYCLE ) {
      std::cout << std::endl << " *** Start Nonlinear " << _solverType << " Full-Cycle ***" << std::endl;
      grid0 = 1;
    }
    else if( _mg_type == V_CYCLE ) {
      std::cout << std::endl << " *** Start Nonlinear  " << _solverType << " V-Cycle ***" << std::endl;
      grid0 = _gridn;
    }
    else if( _mg_type == M_CYCLE ) {
      std::cout << std::endl << " *** Start Nonlinear  " << _solverType << " Mixed-Cycle ***" << std::endl;
      grid0 = _gridr;
    }
    else {
      std::cout << "wrong " << _solverType << " type for this solver " << std::endl;
      abort();
    }

    unsigned AMRCounter = 0;

    for( unsigned igridn = grid0; igridn <= _gridn; igridn++ ) {   //_igridn
      std::cout << std::endl << "   ****** Start Level Max " << igridn << " ******" << std::endl;
      clock_t start_nl_time = clock();

      bool ThisIsAMR = ( _mg_type == F_CYCLE && _AMRtest &&  AMRCounter < _maxAMRlevels && igridn == _gridn ) ? 1 : 0;

      if( ThisIsAMR ) _solution[igridn - 1]->InitAMREps();

      for( unsigned nonLinearIterator = 0; nonLinearIterator < _n_max_nonlinear_iterations; nonLinearIterator++ ) {

        std::cout << std::endl << "   ********* Nonlinear iteration " << nonLinearIterator + 1 << " *********" << std::endl;

        clock_t start_assembly_time = clock();
        _levelToAssemble = igridn - 1u; //Be carefull!!!! this is needed in the _assemble_function
        _LinSolver[igridn - 1u]->SetResZero();
        _assembleMatrix = _buildSolver;
        _assemble_system_function( _equation_systems );

        if( _buildSolver ){
          _MGmatrixFineReuse = ( 0 == nonLinearIterator ) ? false : true;
          _MGmatrixCoarseReuse = ( igridn - grid0 > 0 ) ?  true : _MGmatrixFineReuse;

	  clock_t mg_proj_mat_time = clock();
          for( unsigned i = igridn - 1u; i > 0; i-- ) {
            if( _RR[i] ) {
              if( i == igridn - 1u )
                _LinSolver[i - 1u]->_KK->matrix_ABC( *_RR[i], *_LinSolver[i]->_KK, *_PP[i], _MGmatrixFineReuse );
              else
                _LinSolver[i - 1u]->_KK->matrix_ABC( *_RR[i], *_LinSolver[i]->_KK, *_PP[i], _MGmatrixCoarseReuse );
            }
            else {
              if( i == igridn - 1u )
                _LinSolver[i - 1u]->_KK->matrix_PtAP( *_PP[i], *_LinSolver[i]->_KK, _MGmatrixFineReuse );
              else
                _LinSolver[i - 1u]->_KK->matrix_PtAP( *_PP[i], *_LinSolver[i]->_KK, _MGmatrixCoarseReuse );
            }
          }
            std::cout << "   ********* Level Max " << igridn << " MG PROJECTION MATRICES TIME:\t" \
	    << static_cast<double>( ( clock() - mg_proj_mat_time ) ) / CLOCKS_PER_SEC << std::endl;
          
          clock_t mg_init_time = clock();  
	  if( _MGsolver ) {
            _LinSolver[igridn - 1u]->MGInit( mgSmootherType, igridn, _outer_ksp_solver.c_str() );

            for( unsigned i = 0; i < igridn; i++ ) {
              if( _RR[i] )
                _LinSolver[i]->MGSetLevel( _LinSolver[igridn - 1u], igridn - 1u, _VariablesToBeSolvedIndex, _PP[i], _RR[i], _npre, _npost );
              else
                _LinSolver[i]->MGSetLevel( _LinSolver[igridn - 1u], igridn - 1u, _VariablesToBeSolvedIndex, _PP[i], _PP[i], _npre, _npost );
            }
          }
          std::cout << "   ********* Level Max " << igridn << " MGINIT TIME:\t" \
	    << static_cast<double>( ( clock() - mg_init_time ) ) / CLOCKS_PER_SEC << std::endl;
        }
        std::cout << "   ********* Level Max " << igridn << " ASSEMBLY TIME:\t" << \
          static_cast<double>( ( clock() - start_assembly_time ) ) / CLOCKS_PER_SEC << std::endl;
        clock_t startUpdateResidualTime = clock();

        for( unsigned updateResidualIterator = 0; updateResidualIterator < _maxNumberOfResidualUpdateIterations; updateResidualIterator++ ) {

          std::cout << "     ********* Linear Cycle + Residual Update iteration " << updateResidualIterator + 1 << std::endl;

          bool thisIsConverged;

          if( _MGsolver ) thisIsConverged = MGVcycle( igridn, mgSmootherType );
          else thisIsConverged = MLVcycle( igridn );

          if( thisIsConverged || updateResidualIterator == _maxNumberOfResidualUpdateIterations - 1 ) break;

          _LinSolver[igridn - 1u]->SetResZero();
          _assembleMatrix = false;
          _assemble_system_function( _equation_systems );

        }

        if( _MGsolver * _buildSolver ) _LinSolver[igridn - 1u]->MGClear();

        double nonLinearEps;
        bool nonLinearIsConverged = IsNonLinearConverged( igridn - 1, nonLinearEps );

        std::cout << "     ********* Linear Cycle + Residual Update-Cycle TIME:\t" << std::setw( 11 ) << std::setprecision( 6 ) << std::fixed
        << static_cast<double>( ( clock() - startUpdateResidualTime ) ) / CLOCKS_PER_SEC << std::endl;

        if( nonLinearIsConverged ) break;

      }

      if( ThisIsAMR ) AddAMRLevel( AMRCounter );

      if( igridn < _gridn ) ProlongatorSol( igridn );

      std::cout << std::endl << "   ****** Nonlinear-Cycle TIME: " << std::setw( 11 ) << std::setprecision( 6 ) << std::fixed
                << static_cast<double>( ( clock() - start_nl_time ) ) / CLOCKS_PER_SEC << std::endl;

      std::cout << std::endl << "   ****** End Level Max " << igridn << " ******" << std::endl;
    }

    std::cout << std::endl << "   *** Nonlinear " << _solverType << " TIME: " << std::setw( 11 ) << std::setprecision( 6 ) << std::fixed
              << static_cast<double>( ( clock() - start_mg_time ) ) / CLOCKS_PER_SEC << std::endl;

  }



} //end namespace femus

