/*=========================================================================

 Program: FEMUS
 Module: LinearImplicitSystem
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <iomanip>
#include "LinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "ElemType.hpp"
#include "MeshRefinement.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"



namespace femus {

  // ********************************************

  LinearImplicitSystem::LinearImplicitSystem (MultiLevelProblem& ml_probl,
                                              const std::string& name_in,
                                              const unsigned int number_in, const LinearEquationSolverType& smoother_type) :
    ImplicitSystem (ml_probl, name_in, number_in, smoother_type),
    _debug_linear (false),
    _n_max_linear_iterations (3),
    _final_linear_residual (1.e20),
    _linearAbsoluteConvergenceTolerance (1.e-08),
    _mg_type (F_CYCLE),
    _npre (1u),
    _npre0 (1u),
    _npost (1u),
    _AMRtest (0),
    _maxAMRlevels (0),
    _AMRnorm (0),
    _AMReighborThresholdValue (0.),
    _smootherType (smoother_type),
    _includeCoarseLevelSmoother (INCLUDE_COARSE_LEVEL_FALSE),
    _MGmatrixFineReuse (false),
    _MGmatrixCoarseReuse (false),
    _printSolverInfo (false),
    _assembleMatrix (true),
    _sparsityPatternMinimumSize (1u),
    _numberOfGlobalVariables (0u) {
    _SparsityPattern.resize (0);
    _mgOuterSolver = GMRES;
    _totalAssemblyTime = 0.;
    _totalSolverTime = 0.;
    

  }

  // ********************************************

  LinearImplicitSystem::~LinearImplicitSystem() {

    for (unsigned ig = 0; ig < _LinSolver.size(); ig++) {
      _LinSolver[ig]->DeletePde();
      delete _LinSolver[ig];
      delete _PP[ig];
      delete _RR[ig];
      delete _PPamr[ig];
      delete _RRamr[ig];
    }

    _NSchurVar_test = 0;
    _numblock_test = 0;
    _numblock_all_test = 0;
    _numberOfGlobalVariables = 0;


  }

  
  
  void LinearImplicitSystem::SetDebugLinear(const bool my_value) {
      
        if ( this->GetMLProb()._ml_sol->GetWriter() != NULL)        _debug_linear = my_value;
        else {std::cout << "SetWriter first" << std::endl; abort(); }
        
  }
  
  
  // ********************************************

  void LinearImplicitSystem::SetSparsityPattern (vector< bool > other_sparcity_pattern) {
    unsigned SolPdeSize2 = _SolSystemPdeIndex.size() * _SolSystemPdeIndex.size();

    if (other_sparcity_pattern.size() != SolPdeSize2) {
      std::cout << "Error! Sparsity Pattern size ( " << other_sparcity_pattern.size() << " ) does not match system PDE size" << std::endl;
      exit (0);
    }

    _SparsityPattern.resize (SolPdeSize2);

    for (int i = 0; i < SolPdeSize2; i++) _SparsityPattern[i] = other_sparcity_pattern[i];
  }

  // ********************************************

  void LinearImplicitSystem::SetSparsityPatternMinimumSize (const unsigned &minimumSize) {
    _sparsityPatternMinimumSize = (minimumSize < 2u) ? 1u : minimumSize;
  }

  // ******************************************

  void LinearImplicitSystem::init() {

    _LinSolver.resize (_gridn);


    if (_includeCoarseLevelSmoother == INCLUDE_COARSE_LEVEL_TRUE) {
      _LinSolver[0] = LinearEquationSolver::build (0, _solution[0], _smootherType).release();
    }
    else {
      _LinSolver[0] = LinearEquationSolver::build (0, _solution[0], FEMuS_DEFAULT).release();
    }
    for (unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i] = LinearEquationSolver::build (i, _solution[i], _smootherType).release();
    }

    if (_sparsityPatternMinimumSize != 1u) {
      for (unsigned i = 0; i < _gridn; i++) {
        _LinSolver[i]->SetSparsityPatternMinimumSize (_sparsityPatternMinimumSize);
      }
    }

    for (unsigned i = 0; i < _gridn; i++) {
      _LinSolver[i]->SetNumberOfGlobalVariables (_numberOfGlobalVariables);
      _LinSolver[i]->InitPde (_SolSystemPdeIndex, _ml_sol->GetSolType(),
                              _ml_sol->GetSolName(), &_solution[i]->_Bdc, _gridn, _SparsityPattern);
    }

    _PP.resize (_gridn);
    _RR.resize (_gridn);
    for (unsigned i = 0; i < _gridn; i++) {
      _PP[i] = NULL;
      _RR[i] = NULL;
    }

    for (unsigned ig = 1; ig < _gridn; ig++) {
      BuildProlongatorMatrix (ig);
    }

    _PPamr.resize (_gridn);
    _RRamr.resize (_gridn);
    for (unsigned i = 0; i < _gridn; i++) {
      _PPamr[i] = NULL;
      _RRamr[i] = NULL;
    }
    for (unsigned ig = 0; ig < _gridn; ig++) {
      if (!_ml_msh->GetLevel (ig)->GetIfHomogeneous()) {
        BuildAmrProlongatorMatrix (ig);
      }
    }
    for (unsigned ig = 1; ig < _gridn; ig++) {
      if (!_ml_msh->GetLevel (ig - 1)->GetIfHomogeneous()) {
        _PP[ig]->matrix_RightMatMult (*_PPamr[ig - 1]);
        if (_RR[ig]) _RR[ig]->matrix_LeftMatMult (*_RRamr[ig - 1]);
      }
    }
    for (unsigned ig = 1; ig < _gridn; ig++) {
      ZeroInterpolatorDirichletNodes (ig);
    }


    _NSchurVar_test = 0;
    _numblock_test = 0;
    _numblock_all_test = 0;
    _richardsonScaleFactorIsSet = false;
    // By default we solve for all the PDE variables
    ClearVariablesToBeSolved();
    AddVariableToBeSolved ("All");
  }

  // ********************************************

  void LinearImplicitSystem::MGsolve (const MgSmootherType& mgSmootherType) {

    _bitFlipCounter = 0;

    clock_t start_mg_time = clock();

    unsigned grid0;

    if (_mg_type == F_CYCLE) {
      std::cout << std::endl << " *** Start Linear Full-Cycle ***" << std::endl;
      grid0 = 0;
    }
    else if (_mg_type == V_CYCLE) {
      std::cout << std::endl << " *** Start Linear V-Cycle ***" << std::endl;
      grid0 = _gridn - 1;
    }
    else {
      std::cout << "wrong CYCLE type for this solver " << std::endl;
      abort();
    }

    unsigned AMRCounter = 0;

    for (unsigned igridn = grid0; igridn < _gridn; igridn++) {    //_igridn
      std::cout << std::endl << " ****** Start Level Max " << igridn + 1 << " ******" << std::endl;


      bool ThisIsAMR = (_mg_type == F_CYCLE && _AMRtest &&  AMRCounter < _maxAMRlevels && igridn == _gridn - 1) ? 1 : 0;
    restart:
      if (ThisIsAMR) _solution[igridn]->InitAMREps();

      clock_t start_preparation_time = clock();

      _levelToAssemble = igridn; //Be carefull!!!! this is needed in the _assemble_function
      _LinSolver[igridn]->SetResZero();
      _assembleMatrix = true;
      clock_t start_assembly_time = clock();
      _assemble_system_function (_equation_systems);
      std::cout << std::endl << " ****** Level Max " << igridn + 1 << " ASSEMBLY TIME:\t" << static_cast<double> ( (clock() - start_assembly_time)) / CLOCKS_PER_SEC << std::endl;


      if (!_ml_msh->GetLevel (igridn)->GetIfHomogeneous()) {
        if (!_RRamr[igridn]) {
          (_LinSolver[igridn]->_RESC)->matrix_mult_transpose (*_LinSolver[igridn]->_RES, *_PPamr[igridn]);
          * (_LinSolver[igridn]->_RES) = * (_LinSolver[igridn]->_RESC);
          _LinSolver[igridn]->SwapMatrices();
          _LinSolver[igridn]->_KK->matrix_PtAP (*_PPamr[igridn], *_LinSolver[igridn]->_KKamr, false);
        }
        else {
          (_LinSolver[igridn]->_RESC)->matrix_mult (*_LinSolver[igridn]->_RES, *_RRamr[igridn]);
          * (_LinSolver[igridn]->_RES) = * (_LinSolver[igridn]->_RESC);
          _LinSolver[igridn]->SwapMatrices();
          _LinSolver[igridn]->_KK->matrix_ABC (*_RRamr[igridn], *_LinSolver[igridn]->_KKamr, *_PPamr[igridn], false);
        }
      }


      _MGmatrixFineReuse = false;
      _MGmatrixCoarseReuse = (igridn - grid0 > 0) ?  true : _MGmatrixFineReuse;
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

      std::cout << std::endl << " ****** Level Max " << igridn + 1 << " PREPARATION TIME:\t" << static_cast<double> ( (clock() - start_preparation_time)) / CLOCKS_PER_SEC << std::endl;

      _LinSolver[igridn]->MGInit (mgSmootherType, igridn + 1, _mgOuterSolver);

      for (unsigned i = 0; i < igridn + 1; i++) {
        unsigned npre = (i == 0) ? _npre0 : _npre;
        unsigned npost = (i == 0) ? 0 : _npost;
        if (_RR[i])
          _LinSolver[i]->MGSetLevel (_LinSolver[igridn], igridn, _VariablesToBeSolvedIndex, _PP[i], _RR[i], npre, npost);
        else
          _LinSolver[i]->MGSetLevel (_LinSolver[igridn], igridn, _VariablesToBeSolvedIndex, _PP[i], _PP[i], npre, npost);
      }

      Vcycle (igridn, mgSmootherType);

      _LinSolver[igridn]->MGClear();


      if (!_ml_msh->GetLevel (igridn)->GetIfHomogeneous()) {
        _LinSolver[igridn]->SwapMatrices();
      }

      if (_bitFlipOccurred && _bitFlipCounter == 1) {
        goto restart;
      }

      if (igridn + 1 < _gridn) ProlongatorSol (igridn + 1);

      if (ThisIsAMR) AddAMRLevel (AMRCounter);

      std::cout << std::endl << " ****** End Level Max " << igridn << " ******" << std::endl;

    }

    std::cout << std::endl << " *** Linear Solver TIME: " << std::setw (11) << std::setprecision (6) << std::fixed
              << static_cast<double> ( (clock() - start_mg_time)) / CLOCKS_PER_SEC << std::endl;

    _totalAssemblyTime += 0.;
    _totalSolverTime += static_cast<double> ( (clock() - start_mg_time)) / CLOCKS_PER_SEC;
  }

  // ********************************************

  bool LinearImplicitSystem::IsLinearConverged (const unsigned igridn) {

    _bitFlipOccurred = false;

    bool conv = true;
    double L2normRes;
//     std::cout << std::endl;

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned indexSol = _SolSystemPdeIndex[k];
      L2normRes       = _solution[igridn]->_Res[indexSol]->l2_norm();
      std::cout << "       *************** Level Max " << igridn + 1 << "  Linear Res  L2norm " << std::scientific << _ml_sol->GetSolutionName (indexSol) << " = " << L2normRes << std::endl;
      if (isnan (L2normRes)) {
        std::cout << "Warning the linear solver did not converge.\n";
        std::cout << "A bit flip may have occurred, let's try to restart the solver!" << std::endl;
        _bitFlipOccurred = true;
      }
      if (L2normRes < _linearAbsoluteConvergenceTolerance && conv == true) {
        conv = true;
      }
      else {
        conv = false;
      }
    }

    if (_bitFlipOccurred) {
      _bitFlipCounter += 1;
    }
    else {
      _bitFlipCounter = 0;
    }

    return conv;

  }

  // *******************************************************

  void LinearImplicitSystem::ProlongatorSol (unsigned gridf) {

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {

      unsigned SolIndex = _SolSystemPdeIndex[k];
      unsigned solType = _ml_sol->GetSolutionType (SolIndex);

      _solution[gridf]->_Sol[SolIndex]->matrix_mult (*_solution[gridf - 1]->_Sol[SolIndex],
                                                     *_msh[gridf]->GetCoarseToFineProjection (solType));
      _solution[gridf]->_Sol[SolIndex]->close();
    }
  }

  // ********************************************

  bool LinearImplicitSystem::Vcycle (const unsigned& level, const MgSmootherType& mgSmootherType) {

    clock_t start_mg_time = clock();

    _LinSolver[level]->SetEpsZero();

    bool linearIsConverged;

    for (unsigned linearIterator = 0; linearIterator < _n_max_linear_iterations; linearIterator++) {  //linear cycle

      std::cout << "       *************** Linear iteration " << linearIterator + 1 << " ***********" << std::endl;
      bool ksp_clean = !linearIterator * _assembleMatrix;
      _LinSolver[level]->MGSolve (ksp_clean);
      _solution[level]->UpdateRes (_SolSystemPdeIndex, _LinSolver[level]->_RES, _LinSolver[level]->KKoffset);
      linearIsConverged = IsLinearConverged (level);

      if (linearIsConverged || _bitFlipOccurred)  break;
    }

    if (!_bitFlipOccurred) {
      if (!_ml_msh->GetLevel (level)->GetIfHomogeneous()) {
        (_LinSolver[level]->_EPSC)->matrix_mult (*_LinSolver[level]->_EPS, *_PPamr[level]);
        * (_LinSolver[level]->_EPS) = * (_LinSolver[level]->_EPSC);
      }
      _solution[level]->UpdateSol (_SolSystemPdeIndex, _LinSolver[level]->_EPS, _LinSolver[level]->KKoffset);
    }
    std::cout << "       *************** Linear-Cycle TIME:\t" << std::setw (11) << std::setprecision (6) << std::fixed
              << static_cast<double> ( (clock() - start_mg_time)) / CLOCKS_PER_SEC << std::endl;
    return linearIsConverged;
  }

  // ********************************************

  void LinearImplicitSystem::Restrictor (const unsigned& gridf) {

    _LinSolver[gridf - 1u]->SetEpsZero();
    _LinSolver[gridf - 1u]->SetResZero();

    if (_RR[gridf]) {
      _LinSolver[gridf - 1u]->_RESC->matrix_mult (*_LinSolver[gridf]->_RES, *_RR[gridf]);
    }
    else {
      _LinSolver[gridf - 1u]->_RESC->matrix_mult_transpose (*_LinSolver[gridf]->_RES, *_PP[gridf]);
    }

    *_LinSolver[gridf - 1u]->_RES += *_LinSolver[gridf - 1u]->_RESC;

  }

  // *******************************************************

  void LinearImplicitSystem::Prolongator (const unsigned& gridf) {

    _LinSolver[gridf]->_EPSC->matrix_mult (*_LinSolver[gridf - 1]->_EPS, *_PP[gridf]);
    _LinSolver[gridf]->UpdateResidual();
    _LinSolver[gridf]->SumEpsCToEps();

  }

  // ********************************************

  void LinearImplicitSystem::AddAMRLevel (unsigned& AMRCounter) {
    bool conv_test = true;

    if (_gridn == 1) {
      MeshRefinement meshcoarser (*_msh[0]);
      meshcoarser.FlagAllElementsToBeRefined();
      conv_test = false;
    }
    else {
      conv_test = _solution[_gridn - 1]->FlagAMRRegionBasedOnErroNormAdaptive (_SolSystemPdeIndex, _AMRthreshold, _AMRnorm, _AMReighborThresholdValue);
    }

    //     if( _AMRnorm == 0 ) {
//       conv_test = _solution[_gridn - 1]->FlagAMRRegionBasedOnl2( _SolSystemPdeIndex, _AMRthreshold );
//     }
//     else if( _AMRnorm == 1 ) {
//       conv_test = _solution[_gridn - 1]->FlagAMRRegionBasedOnSemiNorm( _SolSystemPdeIndex, _AMRthreshold );
//     }

    if (conv_test == false) {
      _ml_msh->AddAMRMeshLevel();
      _ml_sol->AddSolutionLevel();
      AddSystemLevel();
      AMRCounter++;
    }
    else {
      _maxAMRlevels = AMRCounter;
      std::cout << "The AMR solver has converged after " << AMRCounter << " refinements\n";
    }
  }

  // ********************************************

  void LinearImplicitSystem::AddSystemLevel() {

    _equation_systems.AddLevel();

    _msh.resize (_gridn + 1);
    _solution.resize (_gridn + 1);
    _msh[_gridn] = _equation_systems._ml_msh->GetLevel (_gridn);
    _solution[_gridn] = _ml_sol->GetSolutionLevel (_gridn);

    for (int i = 0; i < _gridn; i++) {
      _LinSolver[i]->AddLevel();
    }

    _LinSolver.resize (_gridn + 1);

    _LinSolver[_gridn] = LinearEquationSolver::build (_gridn, _solution[_gridn], _smootherType).release();

    _LinSolver[_gridn]->InitPde (_SolSystemPdeIndex, _ml_sol->GetSolType(),
                                 _ml_sol->GetSolName(), &_solution[_gridn]->_Bdc,  _gridn + 1, _SparsityPattern);

    _PP.resize (_gridn + 1);
    _RR.resize (_gridn + 1);
    _PP[_gridn] = NULL;
    _RR[_gridn] = NULL;
    BuildProlongatorMatrix (_gridn);
    if (!_ml_msh->GetLevel (_gridn - 1)->GetIfHomogeneous()) {
      _PP[_gridn]->matrix_RightMatMult (*_PPamr[_gridn - 1]);
      if (_RR[_gridn]) _RR[_gridn]->matrix_LeftMatMult (*_RRamr[_gridn - 1]);
    }

    _PPamr.resize (_gridn + 1);
    _RRamr.resize (_gridn + 1);
    _PPamr[_gridn] = NULL;
    _RRamr[_gridn] = NULL;
    if (!_ml_msh->GetLevel (_gridn)->GetIfHomogeneous()) {
      BuildAmrProlongatorMatrix (_gridn);
    }

    ZeroInterpolatorDirichletNodes (_gridn);

    _LinSolver[_gridn]->set_solver_type (_finegridsolvertype);
    _LinSolver[_gridn]->SetTolerances (_rtol, _atol, _divtol, _maxits, _restart);
    _LinSolver[_gridn]->set_preconditioner_type (_finegridpreconditioner);
    _LinSolver[_gridn]->PrintSolverInfo (_printSolverInfo);

    if (_numblock_test) {
      unsigned num_block2 = std::min (_num_block, _msh[_gridn]->GetNumberOfElements());
      _LinSolver[_gridn]->SetElementBlockNumber (num_block2);
    }
    else if (_numblock_all_test) {
      _LinSolver[_gridn]->SetElementBlockNumber ("All", _overlap);
    }

    if (_NSchurVar_test) {
      _LinSolver[_gridn]->SetNumberOfSchurVariables (_NSchurVar);
    }

    if (_richardsonScaleFactorIsSet) {
      _LinSolver[_gridn]->SetRichardsonScaleFactor (_richardsonScaleFactor);
      //_LinSolver[_gridn]->SetRichardsonScaleFactor(_richardsonScaleFactor + _richardsonScaleFactorDecrease * (_gridn - 1));
    }

    _gridn++;
  }

  // ********************************************


//---------------------------------------------------------------------------------------------
// This is function sets the AMR options
//---------------------------------------------------------------------------------------------

  void LinearImplicitSystem::SetAMRSetOptions (const std::string& AMR, const unsigned& AMRlevels,
                                               const std::string& AMRnorm, const double& AMRthreshold,
                                               bool (* SetRefinementFlag) (const std::vector < double >& x,
                                                                           const int& ElemGroupNumber, const int& level)) {
    if (!strcmp ("yes", AMR.c_str()) || !strcmp ("YES", AMR.c_str()) || !strcmp ("Yes", AMR.c_str())) {
      _AMRtest = 1;
    }

    _maxAMRlevels = AMRlevels;
    _AMRthreshold.resize (1);
    _AMRthreshold[0] = AMRthreshold;

    if (!strcmp ("L2", AMRnorm.c_str()) || !strcmp ("l2", AMRnorm.c_str())  || !strcmp ("h0", AMRnorm.c_str()) || !strcmp ("H0", AMRnorm.c_str())) {
      _AMRnorm = 0;
    }
    else if (!strcmp ("H1", AMRnorm.c_str()) || !strcmp ("h1", AMRnorm.c_str())) {
      _AMRnorm = 1;
    }
    else {
      std::cout << AMRnorm << " invalid AMRnorm type \n set to default H1 norm" << std::endl;
      _AMRnorm = 1;
    }

    if (SetRefinementFlag == NULL) {
    }
    else {
      _msh[0]->Mesh::_SetRefinementFlag = SetRefinementFlag;
      _msh[0]->Mesh::_IsUserRefinementFunctionDefined = true;
    }
  }


//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the FE matrix to finer grids.
// This is a virtual function overloaded in the class MonolithicFSINonLinearImplicitSystem.
//---------------------------------------------------------------------------------------------------

  void LinearImplicitSystem::BuildProlongatorMatrix (unsigned gridf) {

    if (gridf < 1) {
      std::cout << "Error! In function \"BuildProlongatorMatrix\" argument less then 1" << std::endl;
      exit (0);
    }

    int iproc;
    MPI_Comm_rank (MPI_COMM_WORLD, &iproc);

    LinearEquationSolver* LinSolf = _LinSolver[gridf];
    LinearEquationSolver* LinSolc = _LinSolver[gridf - 1];
    Mesh* mshc = _msh[gridf - 1];
    int nf = LinSolf->KKIndex[LinSolf->KKIndex.size() - 1u];
    int nc = LinSolc->KKIndex[LinSolc->KKIndex.size() - 1u];
    int nf_loc = LinSolf->KKoffset[LinSolf->KKIndex.size() - 1][iproc] - LinSolf->KKoffset[0][iproc];
    int nc_loc = LinSolc->KKoffset[LinSolc->KKIndex.size() - 1][iproc] - LinSolc->KKoffset[0][iproc];

    NumericVector* NNZ_d = NumericVector::build().release();
    NNZ_d->init (*LinSolf->_EPS);
    NNZ_d->zero();

    NumericVector* NNZ_o = NumericVector::build().release();
    NNZ_o->init (*LinSolf->_EPS);
    NNZ_o->zero();

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned SolIndex = _SolSystemPdeIndex[k];
      unsigned  SolType = _ml_sol->GetSolutionType (SolIndex);

      // loop on the coarse grid
      for (int isdom = iproc; isdom < iproc + 1; isdom++) {
        for (int iel = mshc->_elementOffset[isdom]; iel < mshc->_elementOffset[isdom + 1]; iel++) {
          short unsigned ielt = mshc->GetElementType (iel);
          mshc->_finiteElement[ielt][SolType]->GetSparsityPatternSize (*LinSolf, *LinSolc, iel, NNZ_d, NNZ_o, SolIndex, k);
        }
      }
    }

    NNZ_d->close();
    NNZ_o->close();

    unsigned offset = LinSolf->KKoffset[0][iproc];
    vector <int> nnz_d (nf_loc);
    vector <int> nnz_o (nf_loc);

    for (int i = 0; i < nf_loc; i++) {
      nnz_d[i] = static_cast <int> ( (*NNZ_d) (offset + i));
      nnz_o[i] = static_cast <int> ( (*NNZ_o) (offset + i));
    }

    delete NNZ_d;
    delete NNZ_o;

    _PP[gridf] = SparseMatrix::build().release();
    _PP[gridf]->init (nf, nc, nf_loc, nc_loc, nnz_d, nnz_o);

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned SolIndex = _SolSystemPdeIndex[k];
      unsigned  SolType = _ml_sol->GetSolutionType (SolIndex);

      // loop on the coarse grid
      for (int isdom = iproc; isdom < iproc + 1; isdom++) {
        for (int iel = mshc->_elementOffset[isdom]; iel < mshc->_elementOffset[isdom + 1]; iel++) {
          short unsigned ielt = mshc->GetElementType (iel);
          mshc->_finiteElement[ielt][SolType]->BuildProlongation (*LinSolf, *LinSolc, iel, _PP[gridf], SolIndex, k);
        }
      }
    }

    _PP[gridf]->close();
  }


  void LinearImplicitSystem::BuildAmrProlongatorMatrix (unsigned level) {

    int iproc;
    MPI_Comm_rank (MPI_COMM_WORLD, &iproc);

    LinearEquationSolver* LinSol = _LinSolver[level];

    Mesh* mesh = _msh[level];
    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > &amrRestriction = mesh->GetAmrRestrictionMap();

    int n = LinSol->KKIndex[LinSol->KKIndex.size() - 1u];
    int n_loc = LinSol->KKoffset[LinSol->KKIndex.size() - 1][iproc] - LinSol->KKoffset[0][iproc];

    NumericVector* NNZ_d = NumericVector::build().release();
    NNZ_d->init (*LinSol->_EPS);
    NNZ_d->zero();

    NumericVector* NNZ_o = NumericVector::build().release();
    NNZ_o->init (*LinSol->_EPS);
    NNZ_o->zero();

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {

      unsigned solIndex = _SolSystemPdeIndex[k];
      unsigned solType = _ml_sol->GetSolutionType (solIndex);
      unsigned kOffset = LinSol->KKoffset[k][iproc];

      unsigned solOffset = mesh->_dofOffset[solType][iproc];
      unsigned solOffsetp1 = mesh->_dofOffset[solType][iproc + 1];
      for (unsigned i = solOffset; i < solOffsetp1; i++) {
        if (solType > 2 || amrRestriction[solType].find (i) == amrRestriction[solType].end()) {
          NNZ_d->set (kOffset + (i - solOffset), 1);
        }
        else {
          double cnt_d = 0;
          double cnt_o = 0;
          for (std::map<unsigned, double> ::iterator it = amrRestriction[solType][i].begin(); it != amrRestriction[solType][i].end(); it++) {
            if (it->first >= solOffset && it->first < solOffsetp1) {
              cnt_d++;
            }
            else {
              cnt_o++;
            }
          }
          NNZ_d->set (kOffset + (i - solOffset), cnt_d);
          NNZ_o->set (kOffset + (i - solOffset), cnt_o);
        }
      }
    }

    NNZ_d->close();
    NNZ_o->close();

    unsigned offset = LinSol->KKoffset[0][iproc];
    vector <int> nnz_d (n_loc);
    vector <int> nnz_o (n_loc);

    for (int i = 0; i < n_loc; i++) {
      nnz_d[i] = static_cast <int> ( (*NNZ_d) (offset + i));
      nnz_o[i] = static_cast <int> ( (*NNZ_o) (offset + i));
    }

    delete NNZ_d;
    delete NNZ_o;

    _PPamr[level] = SparseMatrix::build().release();
    _PPamr[level]->init (n, n, n_loc, n_loc, nnz_d, nnz_o);

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned solIndex = _SolSystemPdeIndex[k];
      unsigned  solType = _ml_sol->GetSolutionType (solIndex);

      unsigned kOffset = LinSol->KKoffset[k][iproc];

      unsigned solOffset = mesh->_dofOffset[solType][iproc];
      unsigned solOffsetp1 = mesh->_dofOffset[solType][iproc + 1];

      for (unsigned i = solOffset; i < solOffsetp1; i++) {
        unsigned irow = kOffset + (i - solOffset);
        if (solType > 2 || amrRestriction[solType].find (i) == amrRestriction[solType].end()) {
          std::vector <int> col (1, irow);
          double value = 1.;
          _PPamr[level]->insert_row (irow, 1, col, &value);
        }
        else {
          int ncols = amrRestriction[solType][i].size();
          std::vector <int> col (ncols);
          std::vector <double> value (ncols);
          unsigned j = 0;
          for (std::map<unsigned, double> ::iterator it = amrRestriction[solType][i].begin(); it != amrRestriction[solType][i].end(); it++) {
            if (it->first >= solOffset && it->first < solOffsetp1) {
              col[j] = kOffset + (it->first - solOffset);
            }
            else {
              unsigned jproc = _msh[level]->IsdomBisectionSearch (it->first, solType);
              col[j] = LinSol->KKoffset[k][jproc] + (it->first - mesh->_dofOffset[solType][jproc]);
            }
            value[j] = it->second;
            j++;
          }
          _PPamr[level]->insert_row (irow, ncols, col, &value[0]);
        }
      }
    }
    _PPamr[level]->close();
    _PPamr[level]->get_transpose (*_PPamr[level]);
  }

  void LinearImplicitSystem::ZeroInterpolatorDirichletNodes (const unsigned &level) {

    int iproc;
    MPI_Comm_rank (MPI_COMM_WORLD, &iproc);

    // Delete the Dirichlet nodes of the fine level (level):
    // set to zero all the corresponding rows for _PP[level] and columns for _RR[level]

    LinearEquationSolver* LinSol = _LinSolver[level];
    Mesh* mesh = _msh[level];
    Solution* solution = _solution[level];

    unsigned BDCIndexSize = LinSol->KKoffset[LinSol->KKIndex.size() - 1][iproc] - LinSol->KKoffset[0][iproc];
    std::vector < int > dirichletNodeIndex (BDCIndexSize);

    unsigned count = 0;

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned solIndex = _SolSystemPdeIndex[k];
      unsigned  solType = _ml_sol->GetSolutionType (solIndex);

      for (unsigned inode_mts = mesh->_dofOffset[solType][iproc]; inode_mts < mesh->_dofOffset[solType][iproc + 1]; inode_mts++) {
        int local_mts = inode_mts - mesh->_dofOffset[solType][iproc];
        int idof_kk = LinSol->KKoffset[k][iproc] + local_mts;
        double bcvalue = (*solution->_Bdc[solIndex]) (inode_mts);
        if (bcvalue < 1.5) {
          dirichletNodeIndex[count] = idof_kk;
          count++;
        }
      }
    }

    dirichletNodeIndex.resize (count);
    std::vector < PetscInt > (dirichletNodeIndex).swap (dirichletNodeIndex);
    std::sort (dirichletNodeIndex.begin(), dirichletNodeIndex.end());
    _PP[level]->mat_zero_rows (dirichletNodeIndex, 0);

    if (_RR[level]) {
      SparseMatrix *RRt;
      RRt = SparseMatrix::build().release();
      _RR[level]->get_transpose (*RRt);
      RRt->mat_zero_rows (dirichletNodeIndex, 0);
      RRt->get_transpose (*_RR[level]);
      delete RRt;
    }

    // Delete the Dirichlet nodes of the coarse level (level-1):
    // set to zero all the corresponding columns for _PP[level] and rows for _RR[level]

    LinSol = _LinSolver[level - 1];
    mesh = _msh[level - 1];
    solution = _solution[level - 1];

    BDCIndexSize = LinSol->KKoffset[LinSol->KKIndex.size() - 1][iproc] - LinSol->KKoffset[0][iproc];
    dirichletNodeIndex.resize (BDCIndexSize);

    count = 0;

    for (unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned solIndex = _SolSystemPdeIndex[k];
      unsigned  solType = _ml_sol->GetSolutionType (solIndex);

      for (unsigned inode_mts = mesh->_dofOffset[solType][iproc]; inode_mts < mesh->_dofOffset[solType][iproc + 1]; inode_mts++) {
        int local_mts = inode_mts - mesh->_dofOffset[solType][iproc];
        int idof_kk = LinSol->KKoffset[k][iproc] + local_mts;
        double bcvalue = (*solution->_Bdc[solIndex]) (inode_mts);
        if (bcvalue < 1.5) {
          dirichletNodeIndex[count] = idof_kk;
          count++;
        }
      }
    }

    dirichletNodeIndex.resize (count);
    std::vector < PetscInt > (dirichletNodeIndex).swap (dirichletNodeIndex);
    std::sort (dirichletNodeIndex.begin(), dirichletNodeIndex.end());

    SparseMatrix *PPt;
    PPt = SparseMatrix::build().release();
    _PP[level]->get_transpose (*PPt);
    PPt->mat_zero_rows (dirichletNodeIndex, 0);
    PPt->get_transpose (*_PP[level]);
    delete PPt;

    if (_RR[level]) {
      _RR[level]->mat_zero_rows (dirichletNodeIndex, 0);
    }

  }

  // ********************************************

  void LinearImplicitSystem::SetDirichletBCsHandling (const DirichletBCType DirichletMode) {

    if (DirichletMode == PENALTY) {
      _DirichletBCsHandlingMode = 0;
    }
    else { // elimination
      _DirichletBCsHandlingMode = 1;
    }
  }

  // ********************************************

  void LinearImplicitSystem::AddVariableToBeSolved (const char solname[]) {

    if (!strcmp (solname, "All") || !strcmp (solname, "ALL") || !strcmp (solname, "all")) {
      _VariablesToBeSolvedIndex.resize (_SolSystemPdeIndex.size());

      for (unsigned i = 0; i < _SolSystemPdeIndex.size(); i++) {
        _VariablesToBeSolvedIndex[i] = i;
      }
    }
    else {
      unsigned n = _VariablesToBeSolvedIndex.size();
      _VariablesToBeSolvedIndex.resize (n + 1u);
      unsigned varind = _ml_sol->GetIndex (solname);

      for (unsigned i = 0; i < _SolSystemPdeIndex.size(); i++) {
        if (_SolSystemPdeIndex[i] == varind) {
          _VariablesToBeSolvedIndex[n] = i;
          break;
        }

        if (_SolSystemPdeIndex.size() - 1u == i) {
          std::cout << "Error! The variable " << solname << " cannot be added to AddVariableToBeSolved"
                    << " Index because it is not included in the solution variable set." << std::endl;
          std::exit (0);
        }
      }
    }
  }

  // ********************************************

  void LinearImplicitSystem::ClearVariablesToBeSolved() {
    _VariablesToBeSolvedIndex.clear();
  }

  // ********************************************

  void LinearImplicitSystem::SetLinearEquationSolverType (const LinearEquationSolverType LinearEquationSolverType, const CoarseLevelInclude &includeCoarseLevelSmoother) {
    _includeCoarseLevelSmoother = includeCoarseLevelSmoother;
    _smootherType = LinearEquationSolverType;
  }


  void LinearImplicitSystem::PrintSolverInfo (const bool & printInfo) {

    _printSolverInfo = printInfo;

    for (unsigned i = 0; i < _gridn; i++) {
      _LinSolver[i]->PrintSolverInfo (_printSolverInfo);
    }
  }


  // ********************************************

  void LinearImplicitSystem::SetElementBlockNumber (unsigned const& dim_block) {
    _numblock_test = 1;
    const unsigned dim = _msh[0]->GetDimension();
    const unsigned base = pow (2, dim);
    _num_block = pow (base, dim_block);

    for (unsigned i = 1; i < _gridn; i++) {
      unsigned num_block2 = std::min (_num_block, _msh[i]->GetNumberOfElements());
      _LinSolver[i]->SetElementBlockNumber (num_block2);
    }
  }

  // ********************************************

  void LinearImplicitSystem::SetElementBlockNumber (const char all[], const unsigned& overlap) {
    _numblock_all_test = 1;
    _overlap = overlap;

    for (unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i]->SetElementBlockNumber (all, overlap);
    }
  }

  // ********************************************

  void LinearImplicitSystem::SetSolverCoarseGrid (const SolverType &coarseGridSolver) {
    _LinSolver[0]->set_solver_type (coarseGridSolver);
  }


  void LinearImplicitSystem::SetSolverFineGrids (const SolverType &fineGridSolver) {
    _finegridsolvertype = fineGridSolver;

    for (unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i]->set_solver_type (_finegridsolvertype);
    }
  }

  // ********************************************




  void LinearImplicitSystem::SetPreconditionerCoarseGrid (const PreconditionerType &coarseGridPreconditioner) {
    _LinSolver[0]->set_preconditioner_type (coarseGridPreconditioner);
  }


  void LinearImplicitSystem::SetPreconditionerFineGrids (const PreconditionerType &fineGridPreconditioner) {
    _finegridpreconditioner = fineGridPreconditioner;

    for (unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i]->set_preconditioner_type (_finegridpreconditioner);
    }
  }

  // ********************************************

  void LinearImplicitSystem::SetTolerances (const double &rtol, const double &atol,
                                            const double &divtol, const unsigned &maxits, const unsigned &restart) {
    _rtol = rtol;
    _atol = atol;
    _divtol = divtol;
    _maxits = maxits;
    _restart = restart;

    for (unsigned i = 0; i < _gridn; i++) {
      _LinSolver[i]->SetTolerances (_rtol, _atol, _divtol, _maxits, _restart);
    }
  }

  // ********************************************

  void LinearImplicitSystem::SetNumberOfSchurVariables (const unsigned short& NSchurVar) {
    _NSchurVar_test = 1;
    _NSchurVar = NSchurVar;

    for (unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i]->SetNumberOfSchurVariables (_NSchurVar);
    }
  }

  // ********************************************

  void LinearImplicitSystem::SetFieldSplitTree (FieldSplitTree *fieldSplitTree) {
    for (unsigned i = 0; i < _gridn; i++) {
      _LinSolver[i]->SetFieldSplitTree (fieldSplitTree);
    }
  };


  /// @deprecated
// this function is like init but it doesn't call InitPDE
  void LinearImplicitSystem::init_two() {

    _LinSolver.resize (_gridn);

    _LinSolver[0] = LinearEquationSolver::build (0, _solution[0], FEMuS_DEFAULT).release();

    for (unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i] = LinearEquationSolver::build (i, _solution[i], _smootherType).release();
    }

//     for (unsigned i=0; i<_gridn; i++) {
//       _LinSolver[i]->InitPde(_SolSystemPdeIndex,_ml_sol->GetSolType(),
//                           _ml_sol->GetSolName(),&_solution[i]->_Bdc,_gridr,_gridn,_SparsityPattern);
//     }


    _PP.resize (_gridn);
    _RR.resize (_gridn);

    for (unsigned i = 0; i < _gridn; i++) {
      _PP[i] = NULL;
      _RR[i] = NULL;
    }

//
//     for (unsigned ig=1; ig<_gridn; ig++) {
//       BuildProlongatorMatrix(ig);
//     }

    _NSchurVar_test = 0;
    _numblock_test = 0;
    _numblock_all_test = 0;
    _richardsonScaleFactorIsSet = false;
    // By default we solved for all the PDE variables
    ClearVariablesToBeSolved();
    AddVariableToBeSolved ("All");
  }







// =================================================================
//   std::vector<NumericVector *> _b;   //// LinearEquation (each level)   _RESC
//   std::vector<NumericVector *> _res; //// LinearEquation (each level)   _RES
//   std::vector<NumericVector *> _x;   //// LinearEquation (each level)   _EPS
//   std::vector<NumericVector *> _x_old; //// LinearEquation (each level)  _EPSC

//if the residual norm is small enough,exit the cycle, and so also the MGSolve
//this is also the check for the single level solver
//what is the meaning of having multiple cycles for single-grid solver?
//they are all like just a single linear solver loop, where convergence has already been reached,
// and you do another check after the previous one in the linear solver loop
//the big question is:
///@todo why dont you do "res_fine < Eps1"
// instead of  "res_fine < Eps1*(1.+ bNorm_fine)" ???
//because one is for the absolute error and another one is for the relative error
/// This function solves the discrete problem with multigrid solver
  void LinearImplicitSystem::MGSolve (double Eps1,         // tolerance for the linear solver
                                      int MaxIter,           // n iterations
                                      const uint Gamma,     // Control V W cycle
                                      const uint Nc_pre,    // n pre-smoothing cycles
                                      const uint Nc_coarse, // n coarse cycles
                                      const uint Nc_post    // n post-smoothing cycles
                                     ) {

#ifdef DEFAULT_PRINT_INFO
    std::cout << "######### BEGIN MG SOLVE ########" << std::endl;
#endif
    double res_fine;

    _LinSolver[GetGridn() - 1]->_RESC->close();
    double bNorm_fine =     _LinSolver[GetGridn() - 1]->_RESC->l2_norm();
    _LinSolver[GetGridn() - 1]->_EPSC->close();
    double x_old_fine = _LinSolver[GetGridn() - 1]->_EPSC->l2_norm();

#ifdef DEFAULT_PRINT_INFO
    std::cout << " bNorm_fine l2 "     <<  bNorm_fine                     << std::endl;
    std::cout << " bNorm_fine linfty " << _LinSolver[GetGridn() - 1]->_RESC->linfty_norm()  << std::endl;
    std::cout << " xold_fine l2 "      <<  x_old_fine                     << std::endl;
#endif

    // FAS Multigrid (Nested) ---------
    bool NestedMG = false;

    if (NestedMG) {
      _LinSolver[0]->_EPS->zero();
      MGStep (0, 1.e-20, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post);

      //smooth on the coarse level WITH PHYSICAL b !
      //and compute the residual

      for (uint Level = 1; Level < GetGridn(); Level++) {

        _LinSolver[Level]->_EPS->matrix_mult (*_LinSolver[Level - 1]->_EPS, *_PP[Level]);  //**** project the solution

        res_fine = MGStep (Level, Eps1, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post);

      }
    } // NestedMG

    // V or W cycle
    int cycle = 0;
    bool exit_mg = false;

    while (!exit_mg && cycle < MaxIter) {

///std::cout << "@@@@@@@@@@ BEGIN MG CYCLE @@@@@@@@"<< std::endl;

///std::cout << "@@@@@@@@@@ start on the finest level @@@@@@@@"<< std::endl;

      res_fine = MGStep (GetGridn() - 1, Eps1, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post);

///std::cout << "@@@@@@@@@@ back to the finest level @@@@@@@@"<< std::endl;

///std::cout << "@@@@@@@@@@ END MG CYCLE @@@@@@@@"<< std::endl;

      std::cout << "@@@@@@@@@@ CHECK THE RESIDUAL NORM OF THE FINEST LEVEL @@@@@@@@" << std::endl;

      std::cout << "res_fine: " << res_fine << std::endl;
      std::cout << "bNorm_fine: " << bNorm_fine << std::endl;

      if (res_fine < Eps1 * (1. + bNorm_fine)) exit_mg = true;

      cycle++;

#ifdef DEFAULT_PRINT_INFO
      std::cout << " cycle= " << cycle   << " residual= " << res_fine << " \n";
#endif

    }

#ifdef DEFAULT_PRINT_INFO
    std::cout << "######### END MG SOLVE #######" << std::endl;
#endif
    return;
  }

// ====================================================================
/// This function does one multigrid step
// solve Ax=b
// compute the residual res=b- Ax
// restrict the residual R*res
//notice that the A and b for the POST-smoothing are the same
//as for the pre-smoothing


  double LinearImplicitSystem::MGStep (int Level,           // Level
                                       double Eps1,          // Tolerance
                                       int MaxIter,          // n iterations - number of mg cycles
                                       const uint Gamma,     // Control V W cycle
                                       const uint Nc_pre,    // n pre-smoothing smoother iterations
                                       const uint Nc_coarse, // n coarse smoother iterations
                                       const uint Nc_post    // n post-smoothing smoother iterations
                                      ) {


    std::pair<uint, double> rest;

    if (Level == 0) {
///  std::cout << "************ REACHED THE BOTTOM *****************"<< std::endl;

#ifdef DEFAULT_PRINT_CONV
      _LinSolver[Level]->_EPS->close();
      double xNorm0 = _LinSolver[Level]->_EPS->linfty_norm();
      _LinSolver[Level]->_RESC->close();
      double bNorm0 = _LinSolver[Level]->_RESC->linfty_norm();
      _LinSolver[Level]->_KK->close();
      double ANorm0 = _LinSolver[Level]->_KK->l1_norm();
      std::cout << "Level " << Level << " ANorm l1 " << ANorm0 << " bNorm linfty " << bNorm0  << " xNormINITIAL linfty " << xNorm0 << std::endl;
#endif

      rest = _LinSolver[Level]->solve (*_LinSolver[Level]->_KK, *_LinSolver[Level]->_KK, *_LinSolver[Level]->_EPS, *_LinSolver[Level]->_RESC, DEFAULT_EPS_LSOLV_C, Nc_coarse);  //****** smooth on the coarsest level

#ifdef DEFAULT_PRINT_CONV
      std::cout << " Coarse sol : res-norm: " << rest.second << " n-its: " << rest.first << std::endl;
      _LinSolver[Level]->_EPS->close();
      std::cout << " Norm of x after the coarse solution " << _LinSolver[Level]->_EPS->linfty_norm() << std::endl;
#endif

      _LinSolver[Level]->_RES->resid (*_LinSolver[Level]->_RESC, *_LinSolver[Level]->_EPS, *_LinSolver[Level]->_KK);   //************ compute the coarse residual

      _LinSolver[Level]->_RES->close();
      std::cout << "COARSE Level " << Level << " res linfty " << _LinSolver[Level]->_RES->linfty_norm() << " res l2 " << _LinSolver[Level]->_RES->l2_norm() << std::endl;

    }

    else {

///  std::cout << "************ BEGIN ONE PRE-SMOOTHING *****************"<< std::endl;
#ifdef DEFAULT_PRINT_TIME
      std::clock_t start_time = std::clock();
#endif
#ifdef DEFAULT_PRINT_CONV
      _LinSolver[Level]->_EPS->close();
      double xNormpre = _LinSolver[Level]->_EPS->linfty_norm();
      _LinSolver[Level]->_RESC->close();
      double bNormpre = _LinSolver[Level]->_RESC->linfty_norm();
      _LinSolver[Level]->_KK->close();
      double ANormpre = _LinSolver[Level]->_KK->l1_norm();
      std::cout << "Level " << Level << " ANorm l1 " << ANormpre << " bNorm linfty " << bNormpre  << " xNormINITIAL linfty " << xNormpre << std::endl;
#endif

      rest = _LinSolver[Level]->solve (*_LinSolver[Level]->_KK, *_LinSolver[Level]->_KK, *_LinSolver[Level]->_EPS, *_LinSolver[Level]->_RESC, DEFAULT_EPS_PREPOST, Nc_pre);  //****** smooth on the finer level

#ifdef DEFAULT_PRINT_CONV
      std::cout << " Pre Lev: " << Level << ", res-norm: " << rest.second << " n-its: " << rest.first << std::endl;
#endif
#ifdef DEFAULT_PRINT_TIME
      std::clock_t end_time = std::clock();
      std::cout << " time =" << double (end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif

      _LinSolver[Level]->_RES->resid (*_LinSolver[Level]->_RESC, *_LinSolver[Level]->_EPS, *_LinSolver[Level]->_KK);  //********** compute the residual

///    std::cout << "************ END ONE PRE-SMOOTHING *****************"<< std::endl;

///    std::cout << ">>>>>>>> BEGIN ONE DESCENT >>>>>>>>>>"<< std::endl;

      _LinSolver[Level - 1]->_RESC->matrix_mult (*_LinSolver[Level]->_RES, *_RR[Level - 1]);  //****** restrict the residual from the finer grid ( new rhs )
      _LinSolver[Level - 1]->_EPS->close();                                //initial value of x for the presmoothing iterations
      _LinSolver[Level - 1]->_EPS->zero();                                //initial value of x for the presmoothing iterations

      for (uint g = 1; g <= Gamma; g++)
        MGStep (Level - 1, Eps1, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post);  //***** call MGStep for another possible descent

//at this point you have certainly reached the COARSE level

///      std::cout << ">>>>>>>> BEGIN ONE ASCENT >>>>>>>>"<< std::endl;
#ifdef DEFAULT_PRINT_CONV
      _LinSolver[Level - 1]->_RES->close();
      std::cout << "BEFORE PROL Level " << Level << " res linfty " << _LinSolver[Level - 1]->_RES->linfty_norm() << " res l2 " << _LinSolver[Level - 1]->_RES->l2_norm() << std::endl;
#endif

      _LinSolver[Level]->_RES->matrix_mult (*_LinSolver[Level - 1]->_EPS, *_PP[Level]);  //******** project the dx from the coarser grid
#ifdef DEFAULT_PRINT_CONV
      _LinSolver[Level]->_RES->close();
      //here, the _res contains the prolongation of dx, so the new x is x_old + P dx
      std::cout << "AFTER PROL Level " << Level << " res linfty " << _LinSolver[Level]->_RES->linfty_norm() << " res l2 " << _LinSolver[Level]->_RES->l2_norm() << std::endl;
#endif
      _LinSolver[Level]->_EPS->add (*_LinSolver[Level]->_RES);  // adding the coarser residual to x
      //initial value of x for the post-smoothing iterations
      //_b is the same as before
///   std::cout << "************ BEGIN ONE POST-SMOOTHING *****************"<< std::endl;
      // postsmooting (Nc_post)
#ifdef DEFAULT_PRINT_TIME
      start_time = std::clock();
#endif
#ifdef DEFAULT_PRINT_CONV
      _LinSolver[Level]->_EPS->close();
      double xNormpost = _LinSolver[Level]->_EPS->linfty_norm();
      _LinSolver[Level]->_RESC->close();
      double bNormpost = _LinSolver[Level]->_RESC->linfty_norm();
      _LinSolver[Level]->_KK->close();
      double ANormpost = _LinSolver[Level]->_KK->l1_norm();
      std::cout << "Level " << Level << " ANorm l1 " << ANormpost << " bNorm linfty " << bNormpost << " xNormINITIAL linfty " << xNormpost << std::endl;
#endif

      rest = _LinSolver[Level]->solve (*_LinSolver[Level]->_KK, *_LinSolver[Level]->_KK, *_LinSolver[Level]->_EPS, *_LinSolver[Level]->_RESC, DEFAULT_EPS_PREPOST, Nc_post);  //***** smooth on the coarser level

#ifdef DEFAULT_PRINT_CONV
      std::cout << " Post Lev: " << Level << ", res-norm: " << rest.second << " n-its: " << rest.first << std::endl;
#endif
#ifdef DEFAULT_PRINT_TIME
      end_time = std::clock();
      std::cout << " time =" << double (end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif

      _LinSolver[Level]->_RES->resid (*_LinSolver[Level]->_RESC, *_LinSolver[Level]->_EPS, *_LinSolver[Level]->_KK);  //*******  compute the residual

///    std::cout << "************ END ONE POST-SMOOTHING *****************"<< std::endl;

    }

    _LinSolver[Level]->_RES->close();

    return  rest.second;  //it returns the residual norm of whatever level you are in
    // if this is the l2_norm then also the nonlinear solver is computed in the l2 norm
    //WHAT NORM is THIS?!? l2, but PRECONDITIONED!!!

  }




} //end namespace femus






