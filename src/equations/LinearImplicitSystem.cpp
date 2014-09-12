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

#include "LinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "ElemType.hpp"
#include <iomanip>


namespace femus {



// ------------------------------------------------------------
// LinearImplicitSystem implementation
LinearImplicitSystem::LinearImplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in, const MgSmoother & smoother_type) :
  ImplicitSystem (ml_probl, name_in, number_in, smoother_type),
  _n_linear_iterations   (0),
  _n_max_linear_iterations (3),
  _final_linear_residual (1.e20),
  _absolute_convergence_tolerance (1.e-08),
  _mg_type(F_CYCLE),
  _npre(1),
  _npost(1),
  _AMRtest(0),
  _maxAMRlevels(0),
  _AMRnorm(0),
  _AMRthreshold(0.01),
  _SmootherType(smoother_type)
 {
  
 }

LinearImplicitSystem::~LinearImplicitSystem() {
   this->clear(); 
}

void LinearImplicitSystem::clear() {
    for (unsigned ig=0; ig<_gridn; ig++) {
      _LinSolver[ig]->DeletePde();
      delete _LinSolver[ig];
    }  
    _NSchurVar_test=0;
    _numblock_test=0;
    _numblock_all_test=0;
}

void LinearImplicitSystem::init() {
  
    _LinSolver.resize(_gridn);
    
    _LinSolver[0]=LinearEquationSolver::build(0,_msh[0],GMRES_SMOOTHER).release();
    for(unsigned i=1;i<_gridn;i++){
      _LinSolver[i]=LinearEquationSolver::build(i,_msh[i],_SmootherType).release();
    }
    
    for (unsigned i=0; i<_gridn; i++) {
      _LinSolver[i]->InitPde(_SolSystemPdeIndex,_ml_sol->GetSolType(),
			     _ml_sol->GetSolName(),&_solution[i]->_Bdc,_gridr,_gridn);
    }  
    
    for (unsigned ig=1; ig<_gridn; ig++) {
      BuildProlongatorMatrix(ig);
    }
    
    _NSchurVar_test=0;
    _numblock_test=0;   
    _numblock_all_test=0;
    // By default we solved for all the PDE variables
    ClearVariablesToBeSolved();
    AddVariableToBeSolved("All");
}


void LinearImplicitSystem::AddSystemLevel() {
  
    _equation_systems.AddLevel();
    
    
    _msh.resize(_gridn+1);
    _solution.resize(_gridn+1);
    _msh[_gridn]=_equation_systems._ml_msh->GetLevel(_gridn);
    _solution[_gridn]=_ml_sol->GetSolutionLevel(_gridn);
        
    for(int i=0;i<_gridn;i++){
      _LinSolver[i]->AddLevel();
    }
    _LinSolver.resize(_gridn+1);
    
    _LinSolver[_gridn]=LinearEquationSolver::build(_gridn,_msh[_gridn],_SmootherType).release();
    
    _LinSolver[_gridn]->InitPde(_SolSystemPdeIndex,_ml_sol->GetSolType(),
				_ml_sol->GetSolName(),&_solution[_gridn]->_Bdc,_gridr,_gridn+1);    
    BuildProlongatorMatrix(_gridn);
   
    
    _LinSolver[_gridn]->set_solver_type(_finegridsolvertype);
    _LinSolver[_gridn]->set_tolerances(_rtol,_atol,_divtol,_maxits);  
    _LinSolver[_gridn]->set_preconditioner_type(_finegridpreconditioner);
    _LinSolver[_gridn]->SetDirichletBCsHandling(_DirichletBCsHandlingMode);
    
    if(_numblock_test){
      unsigned num_block2 = std::min(_num_block,_msh[_gridn]->GetElementNumber());
      _LinSolver[_gridn]->SetElementBlockNumber(num_block2);
    }
    else if(_numblock_all_test){
      _LinSolver[_gridn]->SetElementBlockNumber("All", _overlap);
    }
    
    if(_NSchurVar_test) {
      _LinSolver[_gridn]->SetNumberOfSchurVariables(_NSchurVar);
    }

    _gridn++;
}



//---------------------------------------------------------------------------------------------------

void LinearImplicitSystem::solve() {
  
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
  
  unsigned AMR_counter=0;
  
  std::pair<int, double> solver_info;
     
  for ( unsigned igridn=igrid0; igridn <= _gridn; igridn++) {   //_igridn
    
    std::cout << std::endl << " ************* Level : " << igridn -1 << " *************\n" << std::endl;

     
    int nonlinear_cycle = 0; // da eliminare anche questo parametro!!!!
    
      // ============== Fine level Assembly ==============
      clock_t start_time = clock();
      _LinSolver[igridn-1u]->SetResZero();
      _LinSolver[igridn-1u]->SetEpsZero();
      bool assemble_matrix = true; //Be carefull!!!! this is needed in the _assemble_function
      
      /// Be careful !!!! adesso stiamo usando _sys_number invece che ipde, da togliere al + presto
      _assemble_system_function(_equation_systems, igridn-1u, igridn-1u, assemble_matrix);  
      
#ifndef NDEBUG       
      std::cout << "Grid: " << igridn-1 << "\t        ASSEMBLY TIME:\t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
#endif
 
      for(_n_linear_iterations = 0; _n_linear_iterations < _n_max_linear_iterations; _n_linear_iterations++) { //linear cycle
	
	//std::cout << std::endl;
	std::cout << " ************* V-Cycle : "<< _n_linear_iterations << " *************" << std::endl;
	
	bool ksp_clean=!_n_linear_iterations;
	
	for (unsigned ig = igridn-1u; ig > 0; ig--) {
	  
	  // ============== Presmoothing ============== 
	  for (unsigned k = 0; k < _npre; k++) {
// 	    solver_info = (_VankaIsSet) ? _LinSolver[ig]->solve(_VankaIndex, _NSchurVar, _Schur, ksp_clean*(!k)) : _LinSolver[ig]->solve(ksp_clean*(!k));
	    solver_info = _LinSolver[ig]->solve(_VariablesToBeSolvedIndex, ksp_clean*(!k));
	  }
	  // ============== Non-Standard Multigrid Restriction ==============
	  start_time = clock();
	  Restrictor(ig, igridn, nonlinear_cycle, _n_linear_iterations, full_cycle);
	  
#ifndef NDEBUG 
	  std::cout << "Grid: " << ig << "-->" << ig-1 << "   RESTRICTION TIME:\t       "<< std::setw(11) << std::setprecision(6) << std::fixed
	  << static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
#endif
	}
       
 	// ============== Coarse Direct Solver ==============
 	//solver_info = ( _VankaIsSet ) ? _LinSolver[0]->solve(_VankaIndex, _NSchurVar, _Schur, ksp_clean) : _LinSolver[0]->solve(ksp_clean);
 	solver_info = _LinSolver[0]->solve(_VariablesToBeSolvedIndex, ksp_clean);
             
 	for (unsigned ig = 1; ig < igridn; ig++) {
 	  
 	  // ============== Standard Prolongation ==============
 	  start_time=clock();
 	  Prolongator(ig);

#ifndef NDEBUG 
 	  std::cout << "Grid: " << ig-1 << "-->" << ig << "  PROLUNGATION TIME:\t       " << std::setw(11) << std::setprecision(6) << std::fixed
 	  << static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
#endif 
	  
 	  // ============== PostSmoothing ==============    
 	  for (unsigned k = 0; k < _npost; k++) {
 	    //solver_info = ( _VankaIsSet ) ? _LinSolver[ig]->solve(_VankaIndex, _NSchurVar, _Schur,ksp_clean * (!_npre) * (!k)) : _LinSolver[ig]->solve(ksp_clean * (!_npre) * (!k) );
	    solver_info = _LinSolver[ig]->solve(_VariablesToBeSolvedIndex, ksp_clean * (!_npre) * (!k));
 	  }
 	}
 	// ============== Update Solution ( _gridr-1 <= ig <= igridn-2 ) ==============
 	for (unsigned ig = _gridr-1; ig < igridn-1; ig++) {  // _gridr
 	  _solution[ig]->SumEpsToSol(_SolSystemPdeIndex, _LinSolver[ig]->_EPS, _LinSolver[ig]->_RES, _LinSolver[ig]->KKoffset );	
 	}
 	
 	_final_linear_residual = solver_info.second;
	
	//std::cout << std::endl;
	std::cout << "Grid: " << igridn-1 << "      RESIDUAL:\t\t      " << std::setw(11) << std::setprecision(6) << std::scientific << 
	_final_linear_residual << std::endl << std::endl;
	// ============== Test for linear Convergence (now we are using only the absolute convergence tolerance)==============
 	if(_SmootherType != VANKA_SMOOTHER){
	  if(_final_linear_residual < _absolute_convergence_tolerance) 
	  break;
	}
      }
      
      // ============== Update Solution ( ig = igridn )==============
      _solution[igridn-1]->SumEpsToSol(_SolSystemPdeIndex, _LinSolver[igridn-1]->_EPS, 
					    _LinSolver[igridn-1]->_RES, _LinSolver[igridn-1]->KKoffset );
   
//       std::cout << std::endl;
//       std::cout <<"GRID: "<<igridn-1<< "\t    EAR RESIDUAL:\t"<< _final_linear_residual << std::endl;

    // ==============  Solution Prolongation ==============
       
    if(_mg_type == F_CYCLE && _AMRtest &&  AMR_counter<_maxAMRlevels && igridn==_gridn){
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
  }

  //std::cout << std::endl;
  std::cout << "\t     SOLVER TIME:\t       " << std::setw(11) << std::setprecision(6) << std::fixed 
  <<static_cast<double>((clock()-start_mg_time))/CLOCKS_PER_SEC << std::endl;
  
}

//---------------------------------------------------------------------------------------------
// This is function sets the AMR options
//---------------------------------------------------------------------------------------------

void LinearImplicitSystem::SetAMRSetOptions(const std::string& AMR, const unsigned &AMRlevels,
					    const std::string& AMRnorm, const double &AMRthreshold,
					    bool (* SetRefinementFlag)(const double &x, const double &y, const double &z,
                                       const int &ElemGroupNumber,const int &level)){
  if ( !strcmp("yes",AMR.c_str()) || !strcmp("YES",AMR.c_str()) || !strcmp("Yes",AMR.c_str()) ) {
    _AMRtest=1;
  }
  _maxAMRlevels = AMRlevels;  
  _AMRthreshold = AMRthreshold;
  if ( !strcmp("l2",AMRnorm.c_str()) ){
    _AMRnorm=0;
  }
  else if ( !strcmp("seminorm",AMRnorm.c_str()) ){
    _AMRnorm=1;
  }
  else {
    std::cout<<AMRnorm<<" invalid AMRnorm type \n set to default l2 norm" <<std::endl;
    _AMRnorm=0;
  }
  
  if(SetRefinementFlag==NULL){    
    }
    else{
      _msh[0]->mesh::_SetRefinementFlag = SetRefinementFlag;
      _msh[0]->mesh::_TestSetRefinementFlag=1;
    }
  
}

//---------------------------------------------------------------------------------------------
// This is a virtual function overloaded in the class MonolithicFSINonLinearImplicitSystem.
//---------------------------------------------------------------------------------------------
void LinearImplicitSystem::Restrictor(const unsigned &gridf, const unsigned &gridn, 
					    const unsigned &non_linear_iteration, const unsigned &linear_iteration, const bool &full_cycle){
      
  _LinSolver[gridf-1u]->SetEpsZero();
  _LinSolver[gridf-1u]->SetResZero();
  
  bool assemble_matrix = (linear_iteration == 0) ? true : false;  //Be carefull!!!! this is needed in the _assemble_function      
  if (gridf>=_gridr) {   //_gridr
    _assemble_system_function(_equation_systems, gridf-1, gridn-1u, assemble_matrix);
  }
     
  bool matrix_reuse=true;
  if(assemble_matrix){
    if (gridf>=_gridr) {  //_gridr
      if (!_LinSolver[gridf-1]->_CC_flag) {
	_LinSolver[gridf-1]->_CC_flag=1;
	_LinSolver[gridf-1]->_CC->matrix_PtAP(*_LinSolver[gridf]->_PP,*_LinSolver[gridf]->_KK,!matrix_reuse);
      } 
      else{
	_LinSolver[gridf-1]->_CC->matrix_PtAP(*_LinSolver[gridf]->_PP,*_LinSolver[gridf]->_KK,matrix_reuse);
      }
      _LinSolver[gridf-1u]->_KK->matrix_add(1.,*_LinSolver[gridf-1u]->_CC,"different_nonzero_pattern");
    } 
    else { //Projection of the Matrix on the lower level
      if (non_linear_iteration==0 && ( full_cycle*(gridf==gridn-1u) || !full_cycle )) {
	_LinSolver[gridf-1]->_KK->matrix_PtAP(*_LinSolver[gridf]->_PP,*_LinSolver[gridf]->_KK,!matrix_reuse);
      }
      else{
	_LinSolver[gridf-1]->_KK->matrix_PtAP(*_LinSolver[gridf]->_PP,*_LinSolver[gridf]->_KK,matrix_reuse);
      }    
    }
  }
      
  _LinSolver[gridf-1u]->_RESC->matrix_mult_transpose(*_LinSolver[gridf]->_RES, *_LinSolver[gridf]->_PP);
  *_LinSolver[gridf-1u]->_RES += *_LinSolver[gridf-1u]->_RESC;
  
}

// *******************************************************
void LinearImplicitSystem::Prolongator(const unsigned &gridf) {
  
  _LinSolver[gridf]->_EPSC->matrix_mult(*_LinSolver[gridf-1]->_EPS,*_LinSolver[gridf]->_PP);
  _LinSolver[gridf]->UpdateResidual();
  _LinSolver[gridf]->SumEpsCToEps();
}

// *******************************************************
void LinearImplicitSystem::ProlongatorSol(unsigned gridf) {

  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    
    unsigned SolIndex=_SolSystemPdeIndex[k];
    unsigned Typeindex=_ml_sol->GetSolutionType(SolIndex);
    
    _solution[gridf]->_Sol[SolIndex]->matrix_mult(*_solution[gridf-1]->_Sol[SolIndex],*_solution[gridf]->_ProjMat[Typeindex]);
    _solution[gridf]->_Sol[SolIndex]->close();     
  }
}

//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the FE matrix to finer grids.
// This is a virtual function overloaded in the class MonolithicFSINonLinearImplicitSystem.
//---------------------------------------------------------------------------------------------------

void LinearImplicitSystem::BuildProlongatorMatrix(unsigned gridf) {

  if (gridf<1) {
    std::cout<<"Error! In function \"BuildProlongatorMatrix\" argument less then 1"<<std::endl;
    exit(0);
  }
  
  int iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  
  int nf= _LinSolver[gridf]->KKIndex[_LinSolver[gridf]->KKIndex.size()-1u];
  int nc= _LinSolver[gridf-1]->KKIndex[_LinSolver[gridf-1]->KKIndex.size()-1u];
  int nf_loc = _LinSolver[gridf]->KKoffset[_LinSolver[gridf]->KKIndex.size()-1][iproc]-_LinSolver[gridf]->KKoffset[0][iproc];
  int nc_loc = _LinSolver[gridf-1]->KKoffset[_LinSolver[gridf-1]->KKIndex.size()-1][iproc]-_LinSolver[gridf-1]->KKoffset[0][iproc];
  _LinSolver[gridf]->_PP = SparseMatrix::build().release();
  _LinSolver[gridf]->_PP->init(nf,nc,nf_loc,nc_loc,27,27);
  
  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    unsigned SolIndex=_SolSystemPdeIndex[k];
    
    // loop on the coarse grid 
    for(int isdom=iproc; isdom<iproc+1; isdom++) {
      for (int iel_mts=_msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom]; 
	   iel_mts < _msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = _msh[gridf-1]->IS_Mts2Gmt_elem[iel_mts];
	if(_msh[gridf-1]->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
    
	  short unsigned ielt=_msh[gridf-1]->el->GetElementType(iel);
	  
	  _equation_systems._ml_msh->_type_elem[ielt][_ml_sol->GetSolutionType(SolIndex)]->BuildProlongation(*_LinSolver[gridf],*_LinSolver[gridf-1],iel,
								 _LinSolver[gridf]->_PP,SolIndex,k);
	
	}
      }
    }
  }
  _LinSolver[gridf]->_PP->close();
}


void LinearImplicitSystem::SetDirichletBCsHandling(const DirichletBCType DirichletMode) {
    
  if (DirichletMode == PENALTY) {
    _DirichletBCsHandlingMode = 0;   
  }
  else { // elimination
    _DirichletBCsHandlingMode = 1;   
  } 
  
  for (unsigned i=0; i<_gridn; i++) {
    _LinSolver[i]->SetDirichletBCsHandling(_DirichletBCsHandlingMode);
  }
}


void LinearImplicitSystem::AddVariableToBeSolved(const char solname[]) {
  
  if(!strcmp(solname,"All") || !strcmp(solname,"ALL") || !strcmp(solname,"all")){
    _VariablesToBeSolvedIndex.resize(_SolSystemPdeIndex.size());
    for (unsigned i=0; i<_SolSystemPdeIndex.size(); i++) {
      _VariablesToBeSolvedIndex[i]=i;
    }
  }
  else{
    unsigned n=_VariablesToBeSolvedIndex.size();
    _VariablesToBeSolvedIndex.resize(n+1u);
    unsigned varind=_ml_sol->GetIndex(solname);

    for (unsigned i=0; i<_SolSystemPdeIndex.size(); i++) {
      if (_SolSystemPdeIndex[i]==varind) {
	_VariablesToBeSolvedIndex[n]=i;
	break;
      }
      if (_SolSystemPdeIndex.size()-1u==i) {
	std::cout<<"Error! The variable "<<solname<<" cannot be added to AddVariableToBeSolved" 
		<<" Index because it is not included in the solution variable set."<<std::endl;
	std::exit(0);
      }
    }
  }
}

void LinearImplicitSystem::ClearVariablesToBeSolved() {
  _VariablesToBeSolvedIndex.clear();
}


void LinearImplicitSystem::SetMgSmoother(const MgSmoother mgsmoother) {
  _SmootherType = mgsmoother;
}
  
  

void LinearImplicitSystem::SetElementBlockNumber(unsigned const dim_block) {
  _numblock_test=1;
  const unsigned dim = _msh[0]->GetDimension();
  const unsigned base = pow(2,dim);
  _num_block = pow(base,dim_block);

  for (unsigned i=1; i<_gridn; i++) {
    unsigned num_block2 = std::min(_num_block,_msh[i]->GetElementNumber());
    _LinSolver[i]->SetElementBlockNumber(num_block2);
  }
}

void LinearImplicitSystem::SetElementBlockNumber(const char all[], const unsigned & overlap) {
  _numblock_all_test=1;
  _overlap=overlap;
  for (unsigned i=1; i<_gridn; i++) {
    _LinSolver[i]->SetElementBlockNumber(all, overlap);
  }
}



void LinearImplicitSystem::SetSolverFineGrids(const SolverType finegridsolvertype) {
  _finegridsolvertype=finegridsolvertype;
  for (unsigned i=1; i<_gridn; i++) {
    _LinSolver[i]->set_solver_type(_finegridsolvertype);
  }
}


void LinearImplicitSystem::SetPreconditionerFineGrids(const PreconditionerType finegridpreconditioner) {
  _finegridpreconditioner=finegridpreconditioner;
  for (unsigned i=1; i<_gridn; i++) {
    _LinSolver[i]->set_preconditioner_type(_finegridpreconditioner);
  }
}


void LinearImplicitSystem::SetTolerances(const double rtol, const double atol,
					       const double divtol, const unsigned maxits) {     
  _rtol=rtol;
  _atol=atol;
  _divtol=divtol;
  _maxits=maxits;  
  for (unsigned i=1; i<_gridn; i++) {
    _LinSolver[i]->set_tolerances(_rtol,_atol,_divtol,_maxits);  
  }
}


void LinearImplicitSystem::SetNumberOfSchurVariables(const unsigned short &NSchurVar){
   _NSchurVar_test=1;
   _NSchurVar=NSchurVar;
    for (unsigned i=1; i<_gridn; i++) {
     _LinSolver[i]->SetNumberOfSchurVariables(_NSchurVar);
   }
   
   
}

} //end namespace femus



