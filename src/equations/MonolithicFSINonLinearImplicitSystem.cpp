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
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "ElemType.hpp"



// ------------------------------------------------------------
// NonLinearImplicitSystem implementation
MonolithicFSINonLinearImplicitSystem::MonolithicFSINonLinearImplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in) :
  NonLinearImplicitSystem (ml_probl, name_in, number_in)
{
}

MonolithicFSINonLinearImplicitSystem::~MonolithicFSINonLinearImplicitSystem() {
   this->clear(); 
}

void MonolithicFSINonLinearImplicitSystem::clear() {
}

void MonolithicFSINonLinearImplicitSystem::init() {
  Parent::init();
}

void MonolithicFSINonLinearImplicitSystem::Restrictor(const unsigned &gridf, const unsigned &gridn, 
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
	_LinSolver[gridf-1]->_CC->matrix_ABC(*_LinSolver[gridf]->_RR,*_LinSolver[gridf]->_KK,*_LinSolver[gridf]->_PP,!matrix_reuse);
      } 
      else{
	_LinSolver[gridf-1]->_CC->matrix_ABC(*_LinSolver[gridf]->_RR,*_LinSolver[gridf]->_KK,*_LinSolver[gridf]->_PP,matrix_reuse);
      }
      _LinSolver[gridf-1u]->_KK->matrix_add(1.,*_LinSolver[gridf-1u]->_CC,"different_nonzero_pattern");
    } 
    else { //Projection of the Matrix on the lower level
      if (non_linear_iteration==0 && ( full_cycle*(gridf==gridn-1u) || !full_cycle )) {
	_LinSolver[gridf-1]->_KK->matrix_ABC(*_LinSolver[gridf]->_RR,*_LinSolver[gridf]->_KK,*_LinSolver[gridf]->_PP,!matrix_reuse);
      }
      else{ 
	_LinSolver[gridf-1]->_KK->matrix_ABC(*_LinSolver[gridf]->_RR,*_LinSolver[gridf]->_KK,*_LinSolver[gridf]->_PP,matrix_reuse);
      }    
    }
  }
 
  _LinSolver[gridf-1u]->_RESC->matrix_mult(*_LinSolver[gridf]->_RES, *_LinSolver[gridf]->_RR);
  *_LinSolver[gridf-1u]->_RES += *_LinSolver[gridf-1u]->_RESC;
  
}




//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the FE matrix to finer grids 
//---------------------------------------------------------------------------------------------------

void MonolithicFSINonLinearImplicitSystem::BuildProlongatorMatrix(unsigned gridf, const char pdename[]) {

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
  
  _LinSolver[gridf]->_RR = SparseMatrix::build().release();
  _LinSolver[gridf]->_RR->init(nc,nf,nc_loc,nf_loc,27,27); 

  SparseMatrix *RRt;
  RRt = SparseMatrix::build().release();
  RRt->init(nf,nc,nf_loc,nc_loc,27,27);
  
  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    unsigned SolIndex=_SolSystemPdeIndex[k];
    bool TestDisp=0;
    if(_ml_sol->TestIfSolutionIsDisplacemenet(SolIndex) )   TestDisp=1;
        
    // loop on the coarse grid 
    for(int isdom=iproc; isdom<iproc+1; isdom++) {
      for (int iel_mts=_msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom]; 
	   iel_mts < _msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = _msh[gridf-1]->IS_Mts2Gmt_elem[iel_mts];
	if(_msh[gridf-1]->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
    
	  short unsigned ielt=_msh[gridf-1]->el->GetElementType(iel);
	  if(TestDisp){
	    _equation_systems._ml_msh->_type_elem[ielt][_ml_sol->GetSolutionType(SolIndex)]->BuildRestrictionTranspose(*_LinSolver[gridf],*_LinSolver[gridf-1],iel,
 													      RRt,SolIndex,k, TestDisp);
	  }
	  else{
	    _equation_systems._ml_msh->_type_elem[ielt][_ml_sol->GetSolutionType(SolIndex)]->BuildProlongation(*_LinSolver[gridf],*_LinSolver[gridf-1],iel,
												      RRt,SolIndex,k);
	  }
	  
	  _equation_systems._ml_msh->_type_elem[ielt][_ml_sol->GetSolutionType(SolIndex)]->BuildProlongation(*_LinSolver[gridf],*_LinSolver[gridf-1],iel,
												    _LinSolver[gridf]->_PP,SolIndex,k);
	}
      }
    }
  }
  
  _LinSolver[gridf]->_PP->close();
  
  RRt->close();
  RRt->get_transpose( *_LinSolver[gridf]->_RR);
  delete RRt;
  
}





