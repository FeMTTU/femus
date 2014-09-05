/*=========================================================================

 Program: FEMUS
 Module: Solution
 Authors: Eugenio Aulisa
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include <ctime>
#include <fstream>
#include <algorithm>
#include "Solution.hpp"
#include "ElemType.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"


namespace femus {

using std::cout;
using std::endl;

/**
 *  Contructor 
 **/
// ------------------------------------------------------------------
Solution::Solution(mesh *other_msh){    
 _msh = other_msh;
  for(int i=0;i<5;i++){
    _ProjMatFlag[i]=0;
  }
}

/**
 * Destructor
 **/
// ------------------------------------------------------------------
Solution::~Solution() {
  for (unsigned i=0; i<_SolName.size(); i++) {
    delete [] _SolName[i];
  }
}

/**
 * Add a new variable called 'name'
 */
void Solution::AddSolution( const char name[], const FEFamily fefamily, const FEOrder order, 
			    const unsigned& tmorder, const bool &Pde_type) {

  unsigned n=_Sol.size();

  _SolType.resize(n+1u);
  _SolName.resize(n+1u);
  _SolTmOrder.resize(n+1u);
  _family.resize(n+1u);
  _order.resize(n+1u);
  
  _Sol.resize(n+1u);
  _Res.resize(n+1u);
  _Eps.resize(n+1u);
  
  _Bdc.resize(n+1u);
  _ResEpsBdcFlag.resize(n+1u);
 
  _ResEpsBdcFlag[n]=Pde_type;
  _family[n] = fefamily;
  _order[n] = order;
  _SolType[n] = order - ((fefamily==LAGRANGE)?1:0) + fefamily*3;
  _SolTmOrder[n]=tmorder;
  _SolOld.resize(n+1u);
  _SolName[n]=new char [8];
  strcpy(_SolName[n],name);  
  
}

/**
 * Get the solution index for the variable called name
 **/
//-------------------------------------------------------------------
unsigned Solution::GetIndex(const char name[]) const {
  unsigned index=0;
  while (strcmp(_SolName[index],name)) {
    index++;
    if (index==_Res.size()) {
      cout<<"error! invalid name entry GetIndex(...)"<<endl;
      exit(0);
    }
  }
  return index;
}

/**
 * Allocate memory for the variable called name
 **/
// ------------------------------------------------------------------
void Solution::ResizeSolutionVector(const char name[]) {

  unsigned i=GetIndex(name);
  
  _Sol[i] = NumericVector::build().release();
  if(n_processors()==1) { // IF SERIAL
    _Sol[i]->init(_msh->MetisOffset[_SolType[i]][n_processors()],_msh->own_size[_SolType[i]][processor_id()],false,SERIAL);
  } 
  else { // IF PARALLEL
    if(_SolType[i]<3) {
      if(_msh->ghost_size[_SolType[i]][processor_id()]!=0) { 
	_Sol[i]->init(_msh->MetisOffset[_SolType[i]][n_processors()],_msh->own_size[_SolType[i]][processor_id()],
		      _msh->ghost_nd_mts[_SolType[i]][processor_id()],false,GHOSTED);
      } 
      else { 
	std::vector <int> fake_ghost(1,_msh->own_size[_SolType[i]][processor_id()]);
	_Sol[i]->init(_msh->MetisOffset[_SolType[i]][n_processors()],_msh->own_size[_SolType[i]][processor_id()],
		      fake_ghost,false,GHOSTED);
      }
    }
    else { //discontinuous pressure has no ghost nodes
      _Sol[i]->init(_msh->MetisOffset[_SolType[i]][n_processors()],_msh->own_size[_SolType[i]][processor_id()],false,PARALLEL); 
    }
  }
  
  if (_SolTmOrder[i]==2) { // only if the variable is time dependent
    _SolOld[i] = NumericVector::build().release();
    _SolOld[i]->init(*_Sol[i]);
  }
    
  if(_ResEpsBdcFlag[i]) { //only if the variable is a Pde type
    
    _Res[i] = NumericVector::build().release();
    _Res[i]->init(*_Sol[i]);

    _Eps[i] = NumericVector::build().release();
    _Eps[i]->init(*_Sol[i]);
    
    _Bdc[i] = NumericVector::build().release();
    _Bdc[i]->init(*_Sol[i]);
  }
}

/**
 * Deallocate memory for all variables
 **/
// ------------------------------------------------------------------
void Solution::FreeSolutionVectors() {
  for (unsigned i=0; i<_Sol.size(); i++) {
    if(_ResEpsBdcFlag[i]){
      delete _Res[i];
      delete _Eps[i];
      delete _Bdc[i];
    }
    delete _Sol[i];
    if (_SolTmOrder[i]==2) {
      delete _SolOld[i];
    }
  }
  for (unsigned i=0; i<5; i++) {
    if (_ProjMatFlag[i]) {
      delete _ProjMat[i];
    }
  }
}

/**
 * Initialize the coarse variable coordinates X,Y,Z 
 **/
// ------------------------------------------------------------------
void Solution::SetCoarseCoordinates( vector < vector < double> > &vt){
  unsigned indexSol;
  
  indexSol=GetIndex("X");
  *_Sol[indexSol]=vt[0];
  
  indexSol=GetIndex("Y");
  *_Sol[indexSol]=vt[1];
  
  indexSol=GetIndex("Z");
  *_Sol[indexSol]=vt[2];
  
}

/**
 * Update _Sol, _Res and _Eps based on EPS and RES 
 **/
//--------------------------------------------------------------------------------
void Solution::SumEpsToSol(const vector <unsigned> &_SolPdeIndex,  NumericVector* _EPS,  NumericVector* _RES, 
			   const vector <vector <unsigned> > &KKoffset){
 
  PetscScalar zero=0.;
  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    unsigned soltype =  _SolType[indexSol];

    int loc_size   = _Eps[indexSol]->local_size();
    int loc_offset_EPS = KKoffset[k][processor_id()];// - KKoffset[0][processor_id()]; //?????????

    int glob_offset_eps = _msh->MetisOffset[soltype][processor_id()];

    vector <int> index(_msh->own_size[soltype][processor_id()]);
    for(int i=0; i<_msh->own_size[soltype][processor_id()]; i++) {
      index[i]=loc_offset_EPS+i;
    }
    vector <double> valueEPS(_msh->own_size[soltype][processor_id()]);
    _EPS->get(index,valueEPS);
    vector <double> valueRES(_msh->own_size[soltype][processor_id()]);
    _RES->get(index,valueRES);

    for(int i=0; i<_msh->own_size[soltype][processor_id()]; i++) {
      _Eps[indexSol]->set(i+glob_offset_eps,valueEPS[i]);
      if ((*_Bdc[indexSol])(i+glob_offset_eps)>1.1) _Res[indexSol]->set(i+glob_offset_eps,valueRES[i]);
      else _Res[indexSol]->set(i+glob_offset_eps,zero);
    }
    _Res[indexSol]->close();
    _Eps[indexSol]->close();
  }

  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    _Sol[indexSol]->add(*_Eps[indexSol]);
    _Sol[indexSol]->close();
  }
    
}

void Solution::FlagAMRRegionBasedOnRes(const vector <unsigned> &SolIndex){
    
  vector <double> ResMax(SolIndex.size());
  vector <unsigned> SolType(SolIndex.size());
  vector <unsigned> SolEndInd(SolIndex.size());
  vector <unsigned> nve(SolIndex.size());  
  for (unsigned k=0; k<SolIndex.size(); k++) {
    ResMax[k] = _Res[SolIndex[k]]->linfty_norm ();
    SolType[k] = _SolType[SolIndex[k]];
    SolEndInd[k]   = _msh->GetEndIndex(SolType[k]);
    std::cout<< "ResMax of "   << _SolName[SolIndex[k]] << " = " << ResMax[k]<<std::endl;
    std::cout<< "orderInd of " << _SolName[SolIndex[k]] << " = " << SolType[k]<<std::endl;
    std::cout<< "endInd of "   << _SolName[SolIndex[k]] << " = " << SolEndInd[k]<<std::endl;
  }
  
  Solution* AMR = _msh->_coordinate;
  unsigned  AMRIndex= AMR->GetIndex("AMR");
  unsigned  AMRType = AMR->_SolType[AMRIndex];
  unsigned  AMREndInd = _msh->GetEndIndex(AMRType);
  
  
  unsigned nel= _msh->GetElementNumber();
  
  for (int iel=_msh->IS_Mts2Gmt_elem_offset[_iproc]; iel < _msh->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {

    unsigned kel = _msh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=_msh->el->GetElementType(kel);
    
    for (unsigned k=0; k<SolIndex.size(); k++) {
      nve[k]=_msh->el->GetElementDofNumber(kel,SolEndInd[k]);
      for(unsigned i=0; i<nve[k]; i++) {
	unsigned inode=(SolType[k]<3)?(_msh->el->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
	unsigned inode_metis=_msh->GetMetisDof(inode,SolType[k]);
	double soli = (*_Res[SolIndex[k]])(inode_metis);
      } 
    }
  }
  
}


// ------------------------------------------------------------------
void Solution::UpdateSolution() {
  for (unsigned i=0; i<_Sol.size(); i++) {
    // Copy the old vector
    if (_SolTmOrder[i]==2) {
      *(_SolOld[i]) = *(_Sol[i]);
    }
  }
}


} //end namespace femus


