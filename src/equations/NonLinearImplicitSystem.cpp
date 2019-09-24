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

#include <iomanip>
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "NumericVector.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"

namespace femus {



// ------------------------------------------------------------
// NonLinearImplicitSystem implementation
  NonLinearImplicitSystem::NonLinearImplicitSystem(MultiLevelProblem& ml_probl,
      const std::string& name_in,
      const unsigned int number_in, const LinearEquationSolverType& smoother_type) :
    LinearImplicitSystem(ml_probl, name_in, number_in, smoother_type),
    _final_nonlinear_residual(1.e20),
    _n_max_nonlinear_iterations(15),
    _max_nonlinear_convergence_tolerance(1.e-6),
    _maxNumberOfResidualUpdateIterations(1),
    _debug_nonlinear(false),
    _debug_function(NULL),
    _debug_function_is_initialized(false)
  {

  }

  NonLinearImplicitSystem::~NonLinearImplicitSystem() {
  }

  // ********************************************

  void NonLinearImplicitSystem::init() {
    Parent::init();
  }
  
  
  // ********************************************
  void NonLinearImplicitSystem::SetDebugNonlinear(const bool my_value) {
      
        if ( this->GetMLProb()._ml_sol->GetWriter() != NULL)        _debug_nonlinear = my_value;
        else {std::cout << "SetWriter first" << std::endl; abort(); }
        
 }
  

  // ************************MG********************

  bool NonLinearImplicitSystem::HasNonLinearConverged(const unsigned igridn, double &nonLinearEps) {
      
    bool conv = true;
    double L2normEps;
    double L2normSol;
    double L2normEpsDividedSol;
    double L2normRes;

    nonLinearEps = 0.;
    const double absMinNonlinearEps = 1.e-12;
    const double absMinNormSol = 1.e-12;
    const double mindeltaNormSol = 1.e-50;

// we need to store the global vector here
     if (_debug_nonlinear)  {           *(_eps_fine[_nonliniteration]) = *(_LinSolver[_gridn-1]->_EPS);    }
    
    for(unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
        
      unsigned indexSol = _SolSystemPdeIndex[k];
      L2normRes    = _solution[igridn]->_Res[indexSol]->l2_norm();
      L2normEps    = _solution[igridn]->_Eps[indexSol]->l2_norm();
      L2normSol    = _solution[igridn]->_Sol[indexSol]->l2_norm();
      L2normEpsDividedSol = L2normEps / (L2normSol + mindeltaNormSol);

      std::cout << "     ********* Level Max " << igridn + 1 << " Nonlinear Eps_l2norm/Sol_l2norm " << \
                std::scientific << _ml_sol->GetSolutionName(indexSol) << "= " << L2normEpsDividedSol << \
                "  ** Eps_l2norm= " << L2normEps << "  ** Sol_l2norm= " << L2normSol << std::endl;
      nonLinearEps = (nonLinearEps > L2normEpsDividedSol) ? nonLinearEps : L2normEpsDividedSol;

      if((L2normEpsDividedSol < _max_nonlinear_convergence_tolerance || L2normEps < absMinNonlinearEps || L2normSol < absMinNormSol ) && conv == true) {
        conv = true;
      }
      else {
        conv = false;
      }
      
    }
    

    return conv;
  }

  // ********************************************

  void NonLinearImplicitSystem::MGsolve(const MgSmootherType& mgSmootherType) {

    _bitFlipCounter = 0;
    
    clock_t start_mg_time = clock();

    double totalAssembyTime = 0.;

    unsigned grid0;

    if(_mg_type == F_CYCLE) {
      std::cout << std::endl << " *** Start Nonlinear Full-Cycle ***" << std::endl;
      grid0 = 0;
    }
    else if(_mg_type == V_CYCLE) {
      std::cout << std::endl << " *** Start Nonlinear V-Cycle ***" << std::endl;
      grid0 = _gridn - 1;
    }
    else {
      std::cout << "wrong CYCLE type for this solver " << std::endl;
      abort();
    }

    unsigned AMRCounter = 0;

    for(unsigned igridn = grid0; igridn < _gridn; igridn++) {     //_igridn
      std::cout << std::endl << "   ****** Start Level Max " << igridn + 1 << " ******" << std::endl;
      clock_t start_nl_time = clock();

      bool ThisIsAMR = (_mg_type == F_CYCLE && _AMRtest &&  AMRCounter < _maxAMRlevels && igridn == _gridn - 1u) ? 1 : 0;
      
restart:
      if(ThisIsAMR) _solution[igridn]->InitAMREps();

      
      for(unsigned nonLinearIterator = 0; nonLinearIterator < _n_max_nonlinear_iterations; nonLinearIterator++) {

        _nonliniteration = nonLinearIterator;
        
       if (_debug_nonlinear)  {
                   _eps_fine.push_back(NumericVector::build().release());
                   _eps_fine.back()->init(*_LinSolver[_gridn-1]->_EPS);  //I'd say init also fills the vector
                   *(_eps_fine.back()) = *(_LinSolver[_gridn-1]->_EPS);
            }
      
        std::cout << std::endl << "   ********* Nonlinear iteration " << nonLinearIterator + 1 << " *********" << std::endl;

        clock_t start_preparation_time = clock();
        clock_t start_assembly_time = clock();
        _levelToAssemble = igridn; //Be carefull!!!! this is needed in the _assemble_function
        _LinSolver[igridn]->SetResZero();
        _assembleMatrix = _buildSolver;
        _assemble_system_function(_equation_systems);
        std::cout << "   ********* Level Max " << igridn + 1 << " ASSEMBLY TIME:\t" << \
                  static_cast<double>((clock() - start_assembly_time)) / CLOCKS_PER_SEC << std::endl;
	
        if(!_ml_msh->GetLevel(igridn)->GetIfHomogeneous()) {
          if(!_RRamr[igridn]) {
            (_LinSolver[igridn]->_RESC)->matrix_mult_transpose(*_LinSolver[igridn]->_RES, *_PPamr[igridn]);
          }
          else {
            (_LinSolver[igridn]->_RESC)->matrix_mult(*_LinSolver[igridn]->_RES, *_RRamr[igridn]);
          }
          *(_LinSolver[igridn]->_RES) = *(_LinSolver[igridn]->_RESC);
        }

        if(_buildSolver) {

          _MGmatrixFineReuse = (0 == nonLinearIterator) ? false : true;
          _MGmatrixCoarseReuse = (igridn - grid0 > 0) ?  true : _MGmatrixFineReuse;

          if(!_ml_msh->GetLevel(igridn)->GetIfHomogeneous()) {
            _LinSolver[igridn]->SwapMatrices();
            if(!_RRamr[igridn]) {
              _LinSolver[igridn]->_KK->matrix_PtAP(*_PPamr[igridn], *_LinSolver[igridn]->_KKamr, _MGmatrixFineReuse);
            }
            else {
              _LinSolver[igridn]->_KK->matrix_ABC(*_RRamr[igridn], *_LinSolver[igridn]->_KKamr, *_PPamr[igridn], _MGmatrixFineReuse);
            }
          }

          clock_t mg_proj_mat_time = clock();
          for(unsigned i = igridn; i > 0; i--) {
            if(_RR[i]) {
              if(i == igridn)
                _LinSolver[i - 1u]->_KK->matrix_ABC(*_RR[i], *_LinSolver[i]->_KK, *_PP[i], _MGmatrixFineReuse);
              else {
                _LinSolver[i - 1u]->_KK->matrix_ABC(*_RR[i], *_LinSolver[i]->_KK, *_PP[i], _MGmatrixCoarseReuse);
                if(_LinSolver[i - 1u]->_KKamr) {
                  delete _LinSolver[i - 1u]->_KKamr;
                  _LinSolver[i - 1u]->_KKamr = NULL;
                }
              }
            }
            else {
              if(i == igridn)
                _LinSolver[i - 1u]->_KK->matrix_PtAP(*_PP[i], *_LinSolver[i]->_KK, _MGmatrixFineReuse);
              else {
                _LinSolver[i - 1u]->_KK->matrix_PtAP(*_PP[i], *_LinSolver[i]->_KK, _MGmatrixCoarseReuse);
                if(_LinSolver[i - 1u]->_KKamr) {
                  delete _LinSolver[i - 1u]->_KKamr;
                  _LinSolver[i - 1u]->_KKamr = NULL;
                }
              }
            }
          }
          std::cout << "   ********* Level Max " << igridn + 1 << " MG PROJECTION MATRICES TIME:\t" \
                    << static_cast<double>((clock() - mg_proj_mat_time)) / CLOCKS_PER_SEC << std::endl;

          clock_t mg_init_time = clock();
          
            _LinSolver[igridn]->MGInit(mgSmootherType, igridn + 1, _mgOuterSolver);

            for(unsigned i = 0; i <= igridn; i++) {
              unsigned npre = (i == 0)? _npre0 : _npre;  
              unsigned npost = (i == 0)? 0 : _npost;  
              if(_RR[i])
                _LinSolver[i]->MGSetLevel(_LinSolver[igridn], igridn, _VariablesToBeSolvedIndex, _PP[i], _RR[i], npre, npost);
              else
                _LinSolver[i]->MGSetLevel(_LinSolver[igridn], igridn, _VariablesToBeSolvedIndex, _PP[i], _PP[i], npre, npost);
            }
         
          std::cout << "   ********* Level Max " << igridn + 1 << " MGINIT TIME:\t" \
                    << static_cast<double>((clock() - mg_init_time)) / CLOCKS_PER_SEC << std::endl;
        }
        
        totalAssembyTime += static_cast<double>((clock() - start_assembly_time)) / CLOCKS_PER_SEC;
        std::cout << "   ********* Level Max " << igridn + 1 << " PREPARATION TIME:\t" << \
                  static_cast<double>((clock() - start_preparation_time)) / CLOCKS_PER_SEC << std::endl;
        clock_t startUpdateResidualTime = clock();

        for(unsigned updateResidualIterator = 0; updateResidualIterator < _maxNumberOfResidualUpdateIterations; updateResidualIterator++) {

          std::cout << "     ********* Linear Cycle + Residual Update iteration " << updateResidualIterator + 1 << std::endl;

          bool thisHasConverged;
          
          thisHasConverged = Vcycle(igridn, mgSmootherType);
          
          if(thisHasConverged || updateResidualIterator == _maxNumberOfResidualUpdateIterations - 1) break;

          _LinSolver[igridn]->SetResZero();
          _assembleMatrix = false;
          _assemble_system_function(_equation_systems);
          if(!_ml_msh->GetLevel(igridn)->GetIfHomogeneous()) {
            if(!_RRamr[igridn]) {
              (_LinSolver[igridn]->_RESC)->matrix_mult_transpose(*_LinSolver[igridn]->_RES, *_PPamr[igridn]);
            }
            else {
              (_LinSolver[igridn]->_RESC)->matrix_mult(*_LinSolver[igridn]->_RES, *_RRamr[igridn]);
            }
            *(_LinSolver[igridn]->_RES) = *(_LinSolver[igridn]->_RESC);
          }
          if(_bitFlipOccurred) break;
        }

        if(_buildSolver) {
          if(!_ml_msh->GetLevel(igridn)->GetIfHomogeneous()) {
            _LinSolver[igridn]->SwapMatrices();
          }
          _LinSolver[igridn]->MGClear();
        }

        double nonLinearEps;
        bool nonLinearIsConverged = HasNonLinearConverged(igridn, nonLinearEps);

        std::cout << "     ********* Linear Cycle + Residual Update-Cycle TIME:\t" << std::setw(11) << std::setprecision(6) << std::fixed
                  << static_cast<double>((clock() - startUpdateResidualTime)) / CLOCKS_PER_SEC << std::endl;

                  
       if (_debug_nonlinear)  {
          std::vector < std::string > variablesToBePrinted;
          variablesToBePrinted.push_back("All");
          std::ostringstream output_file_name_stream; output_file_name_stream << "biquadratic" << "." << std::setfill('0') << std::setw(2)   << nonLinearIterator; // the "." after biquadratic is needed to see the sequence of files in Paraview as "time steps"

          std::string out_path;
           if (this->GetMLProb().GetFilesHandler() != NULL)  out_path = this->GetMLProb().GetFilesHandler()->GetOutputPath();
	       else                                              out_path = DEFAULT_OUTPUTDIR;
	       
           //print all variables to file
           this->GetMLProb()._ml_sol->GetWriter()->Write(out_path,output_file_name_stream.str().c_str(),variablesToBePrinted);
	       
           //do desired additional computations at the end of each nonlinear iteration
	      if (_debug_function_is_initialized) _debug_function(this->GetMLProb());
          
        }
        
    
        if(nonLinearIsConverged || _bitFlipOccurred) break;

      }  //end nonlinear iterations
      
      _last_nonliniteration = _nonliniteration;
      
      
      if(_bitFlipOccurred && _bitFlipCounter == 1){
	goto restart;
      }
      
      if(igridn + 1 < _gridn) ProlongatorSol(igridn + 1);

      if(ThisIsAMR) AddAMRLevel(AMRCounter);


      std::cout << std::endl << "   ****** Nonlinear-Cycle TIME: " << std::setw(11) << std::setprecision(6) << std::fixed
                << static_cast<double>((clock() - start_nl_time)) / CLOCKS_PER_SEC << std::endl;

      std::cout << std::endl << "   ****** End Level Max " << igridn + 1 << " ******" << std::endl;
    }

    double totalSolverTime = static_cast<double>((clock() - start_mg_time)) / CLOCKS_PER_SEC;
    std::cout << std::endl << "   *** Nonlinear Solver TIME: " << std::setw(11) << std::setprecision(6) << std::fixed
              << totalSolverTime <<  " = assembly TIME( " << totalAssembyTime << " ) + "
              << " solver TIME( " << totalSolverTime - totalAssembyTime << " ) " << std::endl;

    _totalAssemblyTime += totalAssembyTime;
    _totalSolverTime += totalSolverTime - totalAssembyTime;
  }

  
  
  void NonLinearImplicitSystem::compute_convergence_rate() const {
      
      
           const unsigned index_upper = _last_nonliniteration;

    for(unsigned nonLinearIterator = 0; nonLinearIterator < index_upper; nonLinearIterator++) {
         
          
            NumericVector*    eps_fine_temp = NumericVector::build().release();
                   eps_fine_temp->init(*_LinSolver[_gridn-1]->_EPS);
                   
                   eps_fine_temp->close();
                   eps_fine_temp->zero();
         
           const unsigned index_lower = nonLinearIterator + 1;
         
               for(unsigned n = index_upper; n >= index_lower; n--)  *(eps_fine_temp) += *(_eps_fine[n]);
               
          const double  numerator = eps_fine_temp->/*linfty_norm*/l2_norm();
          *(eps_fine_temp) += *(_eps_fine[index_lower - 1]);
          const double denominator = eps_fine_temp->/*linfty_norm*/l2_norm();

         std::cout <<  std::setw(16) << std::setprecision(16) << std::scientific << nonLinearIterator << " " <<  numerator / (denominator * denominator)  << std::endl;
         
     }
      
      
  }
  


} //end namespace femus



