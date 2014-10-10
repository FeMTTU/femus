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


namespace femus {





// ------------------------------------------------------------
// NonLinearImplicitSystem implementation
MonolithicFSINonLinearImplicitSystem::MonolithicFSINonLinearImplicitSystem (MultiLevelProblem& ml_probl,
				const std::string& name_in,
				const unsigned int number_in,const MgSmoother & smoother_type) :
  NonLinearImplicitSystem (ml_probl, name_in, number_in, smoother_type)
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
    
  LinearEquationSolver* LinSolf=_LinSolver[gridf];
  LinearEquationSolver* LinSolc=_LinSolver[gridf-1];
  
  LinSolc->SetEpsZero();
  LinSolc->SetResZero();
  
  bool assemble_matrix = (linear_iteration == 0) ? true : false;  //Be carefull!!!! this is needed in the _assemble_function      
  if (gridf>=_gridr) {   //_gridr
    _assemble_system_function(_equation_systems, gridf-1, gridn-1u, assemble_matrix);
  }
  
  bool matrix_reuse=true;
  if(assemble_matrix){
    if (gridf>=_gridr) {  //_gridr
      if (!LinSolc->_CC_flag) {
	LinSolc->_CC_flag=1;
	LinSolc->_CC->matrix_ABC(*LinSolf->_RR,*LinSolf->_KK,*LinSolf->_PP,!matrix_reuse);
      } 
      else{
	LinSolc->_CC->matrix_ABC(*LinSolf->_RR,*LinSolf->_KK,*LinSolf->_PP,matrix_reuse);
      }
      LinSolc->_KK->matrix_add(1.,*LinSolc->_CC,"different_nonzero_pattern");
    } 
    else { //Projection of the Matrix on the lower level
      if (non_linear_iteration==0 && ( full_cycle*(gridf==gridn-1u) || !full_cycle )) {
	LinSolc->_KK->matrix_ABC(*LinSolf->_RR,*LinSolf->_KK,*LinSolf->_PP,!matrix_reuse);
      }
      else{ 
	LinSolc->_KK->matrix_ABC(*LinSolf->_RR,*LinSolf->_KK,*LinSolf->_PP,matrix_reuse);
      }    
    }
  }
 
  LinSolc->_RESC->matrix_mult(*LinSolf->_RES, *LinSolf->_RR);
  *LinSolc->_RES += *LinSolc->_RESC;
  
}
//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the FE matrix to finer grids 
//---------------------------------------------------------------------------------------------------

void MonolithicFSINonLinearImplicitSystem::BuildProlongatorMatrix(unsigned gridf) {

  if (gridf<1) {
    std::cout<<"Error! In function \"BuildProlongatorMatrix\" argument less then 1"<<std::endl;
    exit(0);
  }
  
  int iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  
  LinearEquationSolver* LinSolf=_LinSolver[gridf];
  LinearEquationSolver* LinSolc=_LinSolver[gridf-1];
  mesh* mshc = _msh[gridf-1];
  int nf= LinSolf->KKIndex[LinSolf->KKIndex.size()-1u];
  int nc= LinSolc->KKIndex[LinSolc->KKIndex.size()-1u];
  int nf_loc = LinSolf->KKoffset[LinSolf->KKIndex.size()-1][iproc]-LinSolf->KKoffset[0][iproc];
  int nc_loc = LinSolc->KKoffset[LinSolc->KKIndex.size()-1][iproc]-LinSolc->KKoffset[0][iproc];
  
  NumericVector *NNZ_d = NumericVector::build().release();
  NNZ_d->init(*LinSolf->_EPS);
  NNZ_d->zero();
    
  NumericVector *NNZ_o = NumericVector::build().release();
  NNZ_o->init(*LinSolf->_EPS);
  NNZ_o->zero();
  
  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    unsigned SolIndex = _SolSystemPdeIndex[k];
    unsigned  SolType = _ml_sol->GetSolutionType(SolIndex);
    // loop on the coarse grid 
    for(int isdom=iproc; isdom<iproc+1; isdom++) {
      for (int iel_mts=mshc->IS_Mts2Gmt_elem_offset[isdom]; iel_mts < mshc->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = mshc->IS_Mts2Gmt_elem[iel_mts];
	if(mshc->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
   	  short unsigned ielt=mshc->el->GetElementType(iel);
	  _equation_systems._ml_msh->_type_elem[ielt][SolType]->GetSparsityPatternSize(*LinSolf,*LinSolc,iel,NNZ_d, NNZ_o,SolIndex,k);
	}
      }
    }
  }
  
  NNZ_d->close();
  NNZ_o->close();
    
  unsigned offset = LinSolf->KKoffset[0][iproc];
  vector <int> nnz_d(nf_loc);
  vector <int> nnz_o(nf_loc);
  for(int i=0; i<nf_loc;i++){
    nnz_d[i]=static_cast <int> ((*NNZ_d)(offset+i));
    nnz_o[i]=static_cast <int> ((*NNZ_o)(offset+i));
  }
            
  delete NNZ_d;
  delete NNZ_o;
    
  LinSolf->_PP = SparseMatrix::build().release();
  LinSolf->_PP->init(nf,nc,nf_loc,nc_loc,nnz_d, nnz_o);
  
  LinSolf->_RR = SparseMatrix::build().release();
  LinSolf->_RR->init(nc,nf,nc_loc,nf_loc,nnz_d,nnz_o); 

  SparseMatrix *RRt;
  RRt = SparseMatrix::build().release();
  RRt->init(nf,nc,nf_loc,nc_loc,27,27);
  
  for (unsigned k=0; k<_SolSystemPdeIndex.size(); k++) {
    unsigned SolIndex=_SolSystemPdeIndex[k];
    unsigned  SolType = _ml_sol->GetSolutionType(SolIndex);
    bool TestDisp=0;
    if(_ml_sol->TestIfSolutionIsDisplacemenet(SolIndex) )   TestDisp=1;
        
    // loop on the coarse grid 
    for(int isdom=iproc; isdom<iproc+1; isdom++) {
      for (int iel_mts=mshc->IS_Mts2Gmt_elem_offset[isdom]; 
	   iel_mts < mshc->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = mshc->IS_Mts2Gmt_elem[iel_mts];
	if(mshc->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
    
	  short unsigned ielt=mshc->el->GetElementType(iel);
	  if(TestDisp){
	    _equation_systems._ml_msh->_type_elem[ielt][SolType]->BuildRestrictionTranspose(*LinSolf,*LinSolc,iel,RRt,SolIndex,k,TestDisp);
	  }
	  else{
	    _equation_systems._ml_msh->_type_elem[ielt][SolType]->BuildProlongation(*LinSolf,*LinSolc,iel, RRt,SolIndex,k);
	  }	  
	  _equation_systems._ml_msh->_type_elem[ielt][SolType]->BuildProlongation(*LinSolf,*LinSolc,iel, LinSolf->_PP,SolIndex,k);
	}
      }
    }
  }
  
  LinSolf->_PP->close();
  
  RRt->close();
  RRt->get_transpose( *LinSolf->_RR);
  delete RRt;
  
}



} //end namespace femus

