/*=========================================================================

 Program: FEMUS
 Module: NonLinearImplicitSystemWithPrimalDualActiveSetMethod
 Authors: Giorgio Bornia, Saikanth Ratnavale

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <iomanip>
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "LinearEquationSolver.hpp"
#include "NumericVector.hpp"



namespace femus {



// ------------------------------------------------------------
// NonLinearImplicitSystemWithPrimalDualActiveSetMethod implementation
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod::NonLinearImplicitSystemWithPrimalDualActiveSetMethod (MultiLevelProblem& ml_probl,
      const std::string& name_in,
      const unsigned int number_in,
      const LinearEquationSolverType & smoother_type) :
    NonLinearImplicitSystem (ml_probl, name_in, number_in, smoother_type)   {

  }


  // ********************************************

  void NonLinearImplicitSystemWithPrimalDualActiveSetMethod::MGsolve (const MgSmootherType& mgSmootherType) {

    _bitFlipCounter = 0;

    clock_t start_mg_time = clock();

    double totalAssembyTime = 0.;

    unsigned grid0;

    if (_mg_type == F_CYCLE) {
      std::cout << std::endl << " *** Start Nonlinear Full-Cycle ***" << std::endl;
      grid0 = 0;
    }
    else if (_mg_type == V_CYCLE) {
      std::cout << std::endl << " *** Start Nonlinear V-Cycle ***" << std::endl;
      grid0 = _gridn - 1;
    }
    else {
      std::cout << "wrong CYCLE type for this solver " << std::endl;
      abort();
    }

    unsigned AMRCounter = 0;

    for (unsigned igridn = grid0; igridn < _gridn; igridn++) {    //_igridn
      std::cout << std::endl << "   ****** Start Level Max " << igridn + 1 << " ******" << std::endl;
      clock_t start_nl_time = clock();

      bool ThisIsAMR = (_mg_type == F_CYCLE && _AMRtest &&  AMRCounter < _maxAMRlevels && igridn == _gridn - 1u) ? 1 : 0;
    restart:
      if (ThisIsAMR) _solution[igridn]->InitAMREps();


      for (unsigned nonLinearIterator = 0; nonLinearIterator < _n_max_nonlinear_iterations; nonLinearIterator++) {

        _nonliniteration = nonLinearIterator;

        if (_debug_nonlinear)  {
          _eps_fine.push_back (NumericVector::build().release());
          _eps_fine[_eps_fine.size() - 1]->init (*_LinSolver[_gridn - 1]->_EPS);
        }

        std::cout << std::endl << "   ********* Nonlinear iteration " << nonLinearIterator + 1 << " *********" << std::endl;

        clock_t start_preparation_time = clock();
        clock_t start_assembly_time = clock();
        _levelToAssemble = igridn; //Be carefull!!!! this is needed in the _assemble_function
        _LinSolver[igridn]->SetResZero();
        _assembleMatrix = _buildSolver;
        _assemble_system_function (_equation_systems);
        std::cout << "   ********* Level Max " << igridn + 1 << " ASSEMBLY TIME:\t" << \
                  static_cast<double> ( (clock() - start_assembly_time)) / CLOCKS_PER_SEC << std::endl;

        if (!_ml_msh->GetLevel (igridn)->GetIfHomogeneous()) {
          if (!_RRamr[igridn]) {
            (_LinSolver[igridn]->_RESC)->matrix_mult_transpose (*_LinSolver[igridn]->_RES, *_PPamr[igridn]);
          }
          else {
            (_LinSolver[igridn]->_RESC)->matrix_mult (*_LinSolver[igridn]->_RES, *_RRamr[igridn]);
          }
          * (_LinSolver[igridn]->_RES) = * (_LinSolver[igridn]->_RESC);
        }

        if (_buildSolver) {

          _MGmatrixFineReuse = (0 == nonLinearIterator) ? false : true;
          _MGmatrixCoarseReuse = (igridn - grid0 > 0) ?  true : _MGmatrixFineReuse;

          if (!_ml_msh->GetLevel (igridn)->GetIfHomogeneous()) {
            _LinSolver[igridn]->SwapMatrices();
            if (!_RRamr[igridn]) {
              _LinSolver[igridn]->_KK->matrix_PtAP (*_PPamr[igridn], *_LinSolver[igridn]->_KKamr, _MGmatrixFineReuse);
            }
            else {
              _LinSolver[igridn]->_KK->matrix_ABC (*_RRamr[igridn], *_LinSolver[igridn]->_KKamr, *_PPamr[igridn], _MGmatrixFineReuse);
            }
          }

          clock_t mg_proj_mat_time = clock();
          for (unsigned i = igridn; i > 0; i--) {
            if (_RR[i]) {
              if (i == igridn)
                _LinSolver[i - 1u]->_KK->matrix_ABC (*_RR[i], *_LinSolver[i]->_KK, *_PP[i], _MGmatrixFineReuse);
              else {
                _LinSolver[i - 1u]->_KK->matrix_ABC (*_RR[i], *_LinSolver[i]->_KK, *_PP[i], _MGmatrixCoarseReuse);
                if (_LinSolver[i - 1u]->_KKamr) {
                  delete _LinSolver[i - 1u]->_KKamr;
                  _LinSolver[i - 1u]->_KKamr = NULL;
                }
              }
            }
            else {
              if (i == igridn)
                _LinSolver[i - 1u]->_KK->matrix_PtAP (*_PP[i], *_LinSolver[i]->_KK, _MGmatrixFineReuse);
              else {
                _LinSolver[i - 1u]->_KK->matrix_PtAP (*_PP[i], *_LinSolver[i]->_KK, _MGmatrixCoarseReuse);
                if (_LinSolver[i - 1u]->_KKamr) {
                  delete _LinSolver[i - 1u]->_KKamr;
                  _LinSolver[i - 1u]->_KKamr = NULL;
                }
              }
            }
          }
          std::cout << "   ********* Level Max " << igridn + 1 << " MG PROJECTION MATRICES TIME:\t" \
                    << static_cast<double> ( (clock() - mg_proj_mat_time)) / CLOCKS_PER_SEC << std::endl;

          clock_t mg_init_time = clock();

          _LinSolver[igridn]->MGInit (mgSmootherType, igridn + 1, _mgOuterSolver);

          for (unsigned i = 0; i <= igridn; i++) {
            if (_RR[i])
              _LinSolver[i]->MGSetLevel (_LinSolver[igridn], igridn, _VariablesToBeSolvedIndex, _PP[i], _RR[i], _npre, _npost);
            else
              _LinSolver[i]->MGSetLevel (_LinSolver[igridn], igridn, _VariablesToBeSolvedIndex, _PP[i], _PP[i], _npre, _npost);
          }

          std::cout << "   ********* Level Max " << igridn + 1 << " MGINIT TIME:\t" \
                    << static_cast<double> ( (clock() - mg_init_time)) / CLOCKS_PER_SEC << std::endl;
        }
        totalAssembyTime += static_cast<double> ( (clock() - start_assembly_time)) / CLOCKS_PER_SEC;
        std::cout << "   ********* Level Max " << igridn + 1 << " PREPARATION TIME:\t" << \
                  static_cast<double> ( (clock() - start_preparation_time)) / CLOCKS_PER_SEC << std::endl;
        clock_t startUpdateResidualTime = clock();

        for (unsigned updateResidualIterator = 0; updateResidualIterator < _maxNumberOfResidualUpdateIterations; updateResidualIterator++) {

          std::cout << "     ********* Linear Cycle + Residual Update iteration " << updateResidualIterator + 1 << std::endl;

          bool thisHasConverged;

          thisHasConverged = Vcycle (igridn, mgSmootherType);

          if (thisHasConverged || updateResidualIterator == _maxNumberOfResidualUpdateIterations - 1) break;

          _LinSolver[igridn]->SetResZero();
          _assembleMatrix = false;
          _assemble_system_function (_equation_systems);
          if (!_ml_msh->GetLevel (igridn)->GetIfHomogeneous()) {
            if (!_RRamr[igridn]) {
              (_LinSolver[igridn]->_RESC)->matrix_mult_transpose (*_LinSolver[igridn]->_RES, *_PPamr[igridn]);
            }
            else {
              (_LinSolver[igridn]->_RESC)->matrix_mult (*_LinSolver[igridn]->_RES, *_RRamr[igridn]);
            }
            * (_LinSolver[igridn]->_RES) = * (_LinSolver[igridn]->_RESC);
          }
          if (_bitFlipOccurred) break;
        }

        if (_buildSolver) {
          if (!_ml_msh->GetLevel (igridn)->GetIfHomogeneous()) {
            _LinSolver[igridn]->SwapMatrices();
          }
          _LinSolver[igridn]->MGClear();

        }

        double nonLinearEps;
        bool nonLinearIsConverged = HasNonLinearConverged (igridn, nonLinearEps);

        std::cout << "     ********* Linear Cycle + Residual Update-Cycle TIME:\t" << std::setw (11) << std::setprecision (6) << std::fixed
                  << static_cast<double> ( (clock() - startUpdateResidualTime)) / CLOCKS_PER_SEC << std::endl;


        if (_debug_nonlinear)  {
          std::vector < std::string > variablesToBePrinted;
          variablesToBePrinted.push_back ("All");
          std::ostringstream output_file_name_stream;
          output_file_name_stream << "biquadratic" << "." << std::setfill ('0') << std::setw (2)   << nonLinearIterator; // the "." after biquadratic is needed to see the sequence of files in Paraview as "time steps"

          std::string out_path;
          if (this->GetMLProb().GetFilesHandler() != NULL)  out_path = this->GetMLProb().GetFilesHandler()->GetOutputPath();
          else                                              out_path = DEFAULT_OUTPUTDIR;

          //print all variables to file
          this->GetMLProb()._ml_sol->GetWriter()->Write (out_path, output_file_name_stream.str().c_str(), variablesToBePrinted);

          //do desired additional computations at the end of each nonlinear iteration
          if (_debug_function_is_initialized) _debug_function (this->GetMLProb());

        }


        // ***************** check active flag sets *******************
        Solution*                sol = this->GetMLProb()._ml_sol->GetSolutionLevel (_levelToAssemble);   // pointer to the solution (level) object
        unsigned int solIndex_act_flag = this->GetMLProb()._ml_sol->GetIndex (_active_flag_name.c_str());

        int compare_return = ( (sol->_SolOld[solIndex_act_flag])->compare (* (sol->_Sol[solIndex_act_flag])));
        bool compare_bool = false;
        if (compare_return == -1) compare_bool = true;
        if (compare_bool && (_nonliniteration  > 0)) {
          std::cout << "Active set did not change at iteration " << _nonliniteration << std::endl;
          break;
        }


        if (nonLinearIsConverged || _bitFlipOccurred) break;

      }


      if (_bitFlipOccurred && _bitFlipCounter == 1) {
        goto restart;
      }

      if (igridn + 1 < _gridn) ProlongatorSol (igridn + 1);

      if (ThisIsAMR) AddAMRLevel (AMRCounter);


      std::cout << std::endl << "   ****** Nonlinear-Cycle TIME: " << std::setw (11) << std::setprecision (6) << std::fixed
                << static_cast<double> ( (clock() - start_nl_time)) / CLOCKS_PER_SEC << std::endl;

      std::cout << std::endl << "   ****** End Level Max " << igridn + 1 << " ******" << std::endl;
    }

    double totalSolverTime = static_cast<double> ( (clock() - start_mg_time)) / CLOCKS_PER_SEC;
    std::cout << std::endl << "   *** Nonlinear Solver TIME: " << std::setw (11) << std::setprecision (6) << std::fixed
              << totalSolverTime <<  " = assembly TIME( " << totalAssembyTime << " ) + "
              << " solver TIME( " << totalSolverTime - totalAssembyTime << " ) " << std::endl;

    _totalAssemblyTime += totalAssembyTime;
    _totalSolverTime += totalSolverTime - totalAssembyTime;
  }



} //end namespace femus



