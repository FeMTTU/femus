/*=========================================================================

 Program: FEMUS
 Module: Writer
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include "Writer.hpp"

#include "mpi.h"

#include "MultiLevelProblem.hpp"
#include "SparseMatrix.hpp"
#include "ElemType.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"



namespace femus {



std::vector<SparseMatrix*> Writer::_ProlQitoQj[3][3];

Writer::Writer(MultiLevelSolution& ml_sol):
  _ml_sol(ml_sol)
{
  _gridn = ml_sol._ml_msh->GetNumberOfLevels();
  _gridr = ml_sol._ml_msh->GetNumberOfGridTotallyRefined();
  
  if(Writer::_ProlQitoQj[0][0].size() == 0)
  {
    BuildProlongatorMatrices();
  }
  _moving_mesh = 0;
}

Writer::~Writer()
{
  for (int igridn=0; igridn<_gridn; igridn++) {
    for (int itype=0; itype<3; itype++) {
      for (int jtype=0; jtype<3; jtype++) {
	if(_ProlQitoQj[itype][jtype][igridn])
	{
          delete _ProlQitoQj[itype][jtype][igridn];
	  _ProlQitoQj[itype][jtype][igridn] = NULL;
	}
      }
    }
  }
  
}


   std::auto_ptr<Writer> Writer::build(const WriterEnum format, MultiLevelSolution * ml_sol)  {
     
      switch (format) {
	case VTK: {
	  std::auto_ptr<Writer>   ap(new VTKWriter(*ml_sol)); return ap;
        }
	case GMV: {
	  std::auto_ptr<Writer>   ap(new GMVWriter(*ml_sol)); return ap;
        }
	case XDMF: {
	  std::auto_ptr<Writer>   ap(new XDMFWriter(*ml_sol)); return ap;
        }
	default: {
	 std::cout << "Format not supported" << std::endl; 
	 abort(); 
	}
	
      } //end switch

      
    }


void Writer::SetMovingMesh(std::vector<std::string>& movvars_in)
{
  _moving_mesh = 1;
  _moving_vars = movvars_in;
}


void Writer::BuildProlongatorMatrices() {

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      _ProlQitoQj[i][j].resize(_gridn);
    }
  }  
  
  for (unsigned igridn=0; igridn<_gridn; igridn++){
    Mesh* msh = _ml_sol._ml_msh->GetLevel(igridn);
    for(int itype=0;itype<3;itype++){
      int ni = msh->MetisOffset[itype][_nprocs];
      int ni_loc = msh->own_size[itype][_iproc];
      for (int jtype=0; jtype<3; jtype++) {
        int nj = msh->MetisOffset[jtype][_nprocs];
	int nj_loc = msh->own_size[itype][_iproc]; 	
	
	NumericVector *NNZ_d = NumericVector::build().release();
	NumericVector *NNZ_o = NumericVector::build().release();
	if(_nprocs==1) { // IF SERIAL
	  NNZ_d->init(ni,ni_loc,false,SERIAL);
	  NNZ_o->init(ni,ni_loc,false,SERIAL);
	} 
	else{
	  NNZ_d->init(ni,ni_loc,false,PARALLEL); 
	  NNZ_o->init(ni,ni_loc,false,PARALLEL); 
	}
	NNZ_d->zero();
    	NNZ_o->zero();
	
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
	  for (int iel_mts=msh->IS_Mts2Gmt_elem_offset[isdom];iel_mts < msh->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++){
	    unsigned iel = msh->IS_Mts2Gmt_elem[iel_mts];
	    short unsigned ielt=msh->el->GetElementType(iel);
            _ml_sol._ml_msh->_finiteElement[ielt][jtype]->GetSparsityPatternSize(*msh, iel, NNZ_d, NNZ_o, itype);	  
	  }
	}
	
	NNZ_d->close();
	NNZ_o->close();
    	
	unsigned offset = msh->MetisOffset[itype][_iproc];
	
	vector <int> nnz_d(ni_loc);
	vector <int> nnz_o(ni_loc);
	for(int i=0; i<ni_loc;i++){
	  nnz_d[i]=static_cast <int> ((*NNZ_d)(offset+i));
	  nnz_o[i]=static_cast <int> ((*NNZ_o)(offset+i));
	}
            
	delete NNZ_d;
	delete NNZ_o;
	
	_ProlQitoQj[itype][jtype][igridn] = SparseMatrix::build().release();
	_ProlQitoQj[itype][jtype][igridn]->init(ni,nj,msh->own_size[itype][_iproc],msh->own_size[jtype][_iproc],nnz_d,nnz_o);
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
	  for (int iel_mts=msh->IS_Mts2Gmt_elem_offset[isdom];iel_mts < msh->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++){
	    unsigned iel = msh->IS_Mts2Gmt_elem[iel_mts];
	    short unsigned ielt=msh->el->GetElementType(iel);
            _ml_sol._ml_msh->_finiteElement[ielt][jtype]->BuildProlongation(*msh, iel, _ProlQitoQj[itype][jtype][igridn], itype);	  
	  }
	}
	_ProlQitoQj[itype][jtype][igridn]->close();
      }
    }
  }
}




} //end namespace femus


