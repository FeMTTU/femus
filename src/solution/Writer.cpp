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

void Writer::SetMovingMesh(std::vector<std::string>& movvars_in)
{
  _moving_mesh = 1;
  _moving_vars = movvars_in;
}


void Writer::BuildProlongatorMatrices() {

  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      Writer::_ProlQitoQj[i][j].resize(_gridn);
    }
  }  
  
  for (unsigned igridn=0; igridn<_gridn; igridn++) {
    for(int itype=0;itype<3;itype++){
      int ni = _ml_sol._ml_msh->GetLevel(igridn)->MetisOffset[itype][_nprocs];
      bool *testnode=new bool [ni];
      for (int jtype=0; jtype<3; jtype++) {
        int nj = _ml_sol._ml_msh->GetLevel(igridn)->MetisOffset[jtype][_nprocs];
	memset(testnode,0,ni*sizeof(bool));
	Writer::_ProlQitoQj[itype][jtype][igridn] = SparseMatrix::build().release();
	Writer::_ProlQitoQj[itype][jtype][igridn]->init(ni,nj,_ml_sol._ml_msh->GetLevel(igridn)->own_size[itype][_iproc],_ml_sol._ml_msh->GetLevel(igridn)->own_size[jtype][_iproc],27,27);

	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
	  for (int iel_mts=_ml_sol._ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom]; 
	       iel_mts < _ml_sol._ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	    unsigned iel = _ml_sol._ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem[iel_mts];
	    short unsigned ielt=_ml_sol._ml_msh->GetLevel(igridn)->el->GetElementType(iel);
            _ml_sol._ml_msh->_type_elem[ielt][jtype]->ProlQitoQj(*_ml_sol._ml_msh->GetLevel(igridn),iel,Writer::_ProlQitoQj[itype][jtype][igridn],testnode,itype);	  
	  }
	}
	Writer::_ProlQitoQj[itype][jtype][igridn]->close();
      }
      delete [] testnode;
    }
  }
}




} //end namespace femus


