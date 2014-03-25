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
#include "string.h"
#include "ElemType.hpp"

// ------------------------------------------------------------
// LinearImplicitSystem implementation
LinearImplicitSystem::LinearImplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in) :
  ImplicitSystem (ml_probl, name_in, number_in),
  _n_linear_iterations   (0),
  _n_max_linear_iterations (3),
  _final_linear_residual (1.e20),
  _absolute_convergence_tolerance (1.e-08),
  _mg_type(F_CYCLE)
{
}

LinearImplicitSystem::~LinearImplicitSystem() {
   this->clear(); 
}

void LinearImplicitSystem::clear() {
    for (unsigned ig=0; ig<_equation_systems.GetNumberOfGrid(); ig++) {
      _LinSolver[ig]->DeletePde();
      delete _LinSolver[ig];
    }  
}

void LinearImplicitSystem::init() {
  CreateSystemPDEStructure();
}

void LinearImplicitSystem::CreateSystemPDEStructure() {
    _LinSolver.resize(_equation_systems.GetNumberOfGrid());
    for(unsigned i=0;i<_equation_systems.GetNumberOfGrid();i++){
      _LinSolver[i]=LinearEquationSolver::build(i,_equation_systems._msh[i]).release();
    }
    
    for (unsigned i=0; i<_equation_systems.GetNumberOfGrid(); i++) {
      //_LinSolver[i]->InitPde(_SolPdeIndex[ipde],SolType,SolName,&_solution[i]->_Bdc,_gridr,_gridn);
      _LinSolver[i]->InitPde(_SolSystemPdeIndex,_equation_systems.SolType,_equation_systems.SolName,&_equation_systems._solution[i]->_Bdc,_equation_systems.GetNumberOfGridNotRefined(),_equation_systems.GetNumberOfGrid());
    }  
    
    for (unsigned ig=1; ig<_equation_systems.GetNumberOfGrid(); ig++) {
      BuildProlongatorMatrix(ig,_sys_name.c_str());
    }
}


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
    igrid0=_equation_systems.GetNumberOfGrid();
  }
  else {
    full_cycle=0;
    igrid0=_equation_systems.GetNumberOfGridNotRefined();
  }
    
  std::pair<int, double> solver_info;
     
  for ( unsigned igridn=igrid0; igridn <= _equation_systems.GetNumberOfGrid(); igridn++) {   //_igridn
    
    std::cout << std::endl << "    ************* Level Max: " << igridn << " *************\n" << std::endl;

     
    int nonlinear_cycle = 0; // da eliminare anche questo parametro!!!!
    
      // ============== Fine level Assembly ==============
      clock_t start_time = clock();
      _LinSolver[igridn-1u]->SetResZero();
      _LinSolver[igridn-1u]->SetEpsZero();
      bool assemble_matrix = true; //Be carefull!!!! this is needed in the _assemble_function
      
      /// Be careful !!!! adesso stiamo usando _sys_number invece che ipde, da togliere al + presto
      _assemble_system_function(_equation_systems, igridn-1u, igridn-1u, assemble_matrix);    
      
      std::cout << "Grid: " << igridn-1 << "\t        ASSEMBLY TIME:\t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;
 
      for(_n_linear_iterations = 0; _n_linear_iterations < _n_max_linear_iterations; _n_linear_iterations++) { //linear cycle
 
	for (unsigned ig = igridn-1u; ig > 0; ig--) {
	  
	  // ============== Presmoothing ============== 
	  for (unsigned k = 0; k < _npre; k++) {
	    solver_info = (_VankaIsSet) ? _LinSolver[ig]->solve(_VankaIndex, _NSchurVar, _Schur) : _LinSolver[ig]->solve();
	  }
	  // ============== Non-Standard Multigrid Restriction ==============
	  start_time = clock();
	  Restrictor(ig, igridn, nonlinear_cycle, _n_linear_iterations, full_cycle);
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
 	  _equation_systems._solution[ig]->SumEpsToSol(_SolSystemPdeIndex, _LinSolver[ig]->_EPS, _LinSolver[ig]->_RES, _LinSolver[ig]->KKoffset );	
 	}
 	
 	_final_linear_residual = solver_info.second;
	// ============== Test for linear Convergence (now we are using only the absolute convergence tolerance)==============
 	if(_final_linear_residual < _absolute_convergence_tolerance) 
	  
	  break;
      }
      
      // ============== Update Solution ( ig = igridn )==============
      _equation_systems._solution[igridn-1]->SumEpsToSol(_SolSystemPdeIndex, _LinSolver[igridn-1]->_EPS, 
							 _LinSolver[igridn-1]->_RES, _LinSolver[igridn-1]->KKoffset );
            
     
      std::cout << std::endl;
      std::cout << "COMPUTATION RESIDUAL: \t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << std::endl;

    // ==============  Solution Prolongation ==============
    if (igridn < _equation_systems.GetNumberOfGrid()) {
      ProlongatorSol(igridn);
    }
  }

  std::cout << "SOLVER TIME:   \t\t\t"<<static_cast<double>((clock()-start_mg_time))/CLOCKS_PER_SEC << std::endl;
  
}

void LinearImplicitSystem::Restrictor(const unsigned &gridf, const unsigned &gridn, 
					    const unsigned &non_linear_iteration, const unsigned &linear_iteration, const bool &full_cycle){
    
  _LinSolver[gridf-1u]->SetEpsZero();
  _LinSolver[gridf-1u]->SetResZero();
  
  bool assemble_matrix = (linear_iteration == 0) ? true : false;  //Be carefull!!!! this is needed in the _assemble_function      
  if (gridf>=_equation_systems.GetNumberOfGridNotRefined()) {   //_gridr
    _assemble_system_function(_equation_systems, gridf-1, gridn-1u, assemble_matrix);
  }
  
  bool matrix_reuse=true;
  if(assemble_matrix){
    if (gridf>=_equation_systems.GetNumberOfGridNotRefined()) {  //_gridr
      if (!_LinSolver[gridf-1]->_CC_flag) {
	_LinSolver[gridf-1]->_CC_flag=1;
	_LinSolver[gridf-1]->_CC->matrix_PtAP(*_LinSolver[gridf]->_PP,*_LinSolver[gridf]->_KK,!matrix_reuse);
	//_LinSolver[gridf-1]->_CC->matrix_ABC(*_LinSolver[gridf]->_RR,*_LinSolver[gridf]->_KK,*_LinSolver[gridf]->_PP,!matrix_reuse);
      } 
      else{
	_LinSolver[gridf-1]->_CC->matrix_PtAP(*_LinSolver[gridf]->_PP,*_LinSolver[gridf]->_KK,matrix_reuse);
	//_LinSolver[gridf-1]->_CC->matrix_ABC(*_LinSolver[gridf]->_RR,*_LinSolver[gridf]->_KK,*_LinSolver[gridf]->_PP,matrix_reuse);
      }
      _LinSolver[gridf-1u]->_KK->matrix_add(1.,*_LinSolver[gridf-1u]->_CC,"different_nonzero_pattern");
    } 
    else { //Projection of the Matrix on the lower level
      if (non_linear_iteration==0 && ( full_cycle*(gridf==gridn-1u) || !full_cycle )) {
	_LinSolver[gridf-1]->_KK->matrix_PtAP(*_LinSolver[gridf]->_PP,*_LinSolver[gridf]->_KK,!matrix_reuse);
	//_LinSolver[gridf-1]->_KK->matrix_ABC(*_LinSolver[gridf]->_RR,*_LinSolver[gridf]->_KK,*_LinSolver[gridf]->_PP,!matrix_reuse);
      }
      else{ 
	_LinSolver[gridf-1]->_KK->matrix_PtAP(*_LinSolver[gridf]->_PP,*_LinSolver[gridf]->_KK,matrix_reuse);
	//_LinSolver[gridf-1]->_KK->matrix_ABC(*_LinSolver[gridf]->_RR,*_LinSolver[gridf]->_KK,*_LinSolver[gridf]->_PP,matrix_reuse);
      }    
    }
  }
      
  _LinSolver[gridf-1u]->_RESC->matrix_mult_transpose(*_LinSolver[gridf]->_RES, *_LinSolver[gridf]->_PP);
  *_LinSolver[gridf-1u]->_RES += *_LinSolver[gridf-1u]->_RESC;
  
//   _LinSolver[gridf-1u]->_RESC->matrix_mult(*_LinSolver[gridf]->_RES, *_LinSolver[gridf]->_RR);
//   *_LinSolver[gridf-1u]->_RES += *_LinSolver[gridf-1u]->_RESC;
  
}

// *******************************************************
void LinearImplicitSystem::Prolongator(const unsigned &gridf) {
  
  _LinSolver[gridf]->_EPSC->matrix_mult(*_LinSolver[gridf-1]->_EPS,*_LinSolver[gridf]->_PP);
  _LinSolver[gridf]->UpdateResidual();
  _LinSolver[gridf]->SumEpsCToEps();

  
}


void LinearImplicitSystem::ProlongatorSol(unsigned gridf) {

  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    unsigned SolIndex=_SolSystemPdeIndex[k];
    unsigned Typeindex=_equation_systems.SolType[SolIndex];
    _equation_systems._solution[gridf]->_Sol[SolIndex]->matrix_mult(*_equation_systems._solution[gridf-1]->_Sol[SolIndex],*_equation_systems._solution[gridf]->_ProjMat[Typeindex]);
    _equation_systems._solution[gridf]->_Sol[SolIndex]->close(); 
    
  }

  
}


//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the FE matrix to finer grids 
//---------------------------------------------------------------------------------------------------

void LinearImplicitSystem::BuildProlongatorMatrix(unsigned gridf, const char pdename[]) {

  unsigned ipde = _sys_number;   //_equation_systems.GetPdeIndex(pdename);
      
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
  
/*  _LinSolver[gridf]->_RR = SparseMatrix::build().release();
  _LinSolver[gridf]->_RR->init(nc,nf,nc_loc,nf_loc,27,27); */   
   
//   SparseMatrix *RRt;
//   if( _TestIfPdeHasDisplacement[ipde]){
//     RRt = SparseMatrix::build().release();
//     RRt->init(nf,nc,nf_loc,nc_loc,27,27);
//   }
  
  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    unsigned SolIndex=_SolSystemPdeIndex[k];
//     bool TestDisp=0;
//     if( _TestIfPdeHasDisplacement[ipde] && _TestIfDisplacement[SolIndex] )   TestDisp=1;
    //TestDisp=0;
    
    // loop on the coarse grid 
    for(int isdom=iproc; isdom<iproc+1; isdom++) {
      for (int iel_mts=_equation_systems._msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom]; 
	   iel_mts < _equation_systems._msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = _equation_systems._msh[gridf-1]->IS_Mts2Gmt_elem[iel_mts];
	if(_equation_systems._msh[gridf-1]->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
    
	  short unsigned ielt=_equation_systems._msh[gridf-1]->el->GetElementType(iel);
	  
// 	  if( _TestIfPdeHasDisplacement[ipde]){
// 	    _equation_systems.type_elem[ielt][SolType[SolIndex]]->BuildRestrictionTranspose(*_LinSolver[ipde][gridf],*_LinSolver[ipde][gridf-1],iel,
// 									  RRt,SolIndex,k, TestDisp);
// 	  }
	  
	  _equation_systems.type_elem[ielt][_equation_systems.SolType[SolIndex]]->BuildProlongation(*_LinSolver[gridf],*_LinSolver[gridf-1],iel,
								 _LinSolver[gridf]->_PP,SolIndex,k);
	
	}
      }
    }
  }
  
  _LinSolver[gridf]->_PP->close();
  
//   if( _TestIfPdeHasDisplacement[ipde]){
//     RRt->close();
//     RRt->get_transpose( *_LinSolver[ipde][gridf]->_RR);
//     delete RRt;
//   }
//   else{
//  _LinSolver[gridf]->_PP->get_transpose( *_LinSolver[gridf]->_RR);
//   }
  
}











