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
 * Add new varible called name
 **/
// ------------------------------------------------------------------
void Solution::AddSolution( const char name[], const char order[],
                                const unsigned& tmorder, const bool &Pde_type) {
  unsigned n=_Sol.size();

  _SolType.resize(n+1u);
  _SolName.resize(n+1u);
  _SolTmOrder.resize(n+1u);
  

  _Sol.resize(n+1u);
  _Res.resize(n+1u);
  _Eps.resize(n+1u);
  
  _Bdc.resize(n+1u);
  _ResEpsBdcFlag.resize(n+1u);
  _ResEpsBdcFlag[n]=Pde_type;
  
  
  _SolTmOrder[n]=tmorder;
  _SolOld.resize(n+1u);

  if (!strcmp(order,"linear")) {
    _SolType[n]=0;
  } else if (!strcmp(order,"quadratic")) {
    _SolType[n]=1;
  } else if (!strcmp(order,"biquadratic")) {
    _SolType[n]=2;
  } else if (!strcmp(order,"constant")) {
    _SolType[n]=3;
  } else if (!strcmp(order,"disc_linear")) {
    _SolType[n]=4;
  } else {
    cout<<"error! invalid order entry in AddSolution(...)"<<endl;
    exit(0);
  }
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
  if(_msh->_nprocs==1) { // IF SERIAL
    _Sol[i]->init(_msh->MetisOffset[_SolType[i]][_msh->_nprocs],_msh->own_size[_SolType[i]][_msh->_iproc],false,SERIAL);
  } 
  else { // IF PARALLEL
    if(_SolType[i]<3) {
      if(_msh->ghost_size[_SolType[i]][_msh->_iproc]!=0) { 
	_Sol[i]->init(_msh->MetisOffset[_SolType[i]][_msh->_nprocs],_msh->own_size[_SolType[i]][_msh->_iproc],
		      _msh->ghost_nd_mts[_SolType[i]][_msh->_iproc],false,GHOSTED);
      } 
      else { 
	std::vector <int> fake_ghost(1,_msh->own_size[_SolType[i]][_msh->_iproc]);
	_Sol[i]->init(_msh->MetisOffset[_SolType[i]][_msh->_nprocs],_msh->own_size[_SolType[i]][_msh->_iproc],
		      fake_ghost,false,GHOSTED);
      }
    }
    else { //discontinuous pressure has no ghost nodes
      _Sol[i]->init(_msh->MetisOffset[_SolType[i]][_msh->_nprocs],_msh->own_size[_SolType[i]][_msh->_iproc],false,PARALLEL); 
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
    int loc_offset_EPS = KKoffset[k][_msh->_iproc];// - KKoffset[0][_msh->_iproc]; //?????????

    int glob_offset_eps = _msh->MetisOffset[soltype][_msh->_iproc];

    vector <int> index(_msh->own_size[soltype][_msh->_iproc]);
    for(int i=0; i<_msh->own_size[soltype][_msh->_iproc]; i++) {
      index[i]=loc_offset_EPS+i;
    }
    vector <double> valueEPS(_msh->own_size[soltype][_msh->_iproc]);
    _EPS->get(index,valueEPS);
    vector <double> valueRES(_msh->own_size[soltype][_msh->_iproc]);
    _RES->get(index,valueRES);

    for(int i=0; i<_msh->own_size[soltype][_msh->_iproc]; i++) {
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
    
  /*
//   PetscVector* EPSp=static_cast<PetscVector*> (_EPS);  //TODO
//   Vec EPS=EPSp->vec(); //TODO
//   PetscVector* RESp=static_cast<PetscVector*> (_RES);  //TODO
//   Vec RES=RESp->vec(); //TODO
// 
// 
//   PetscScalar* R;
//   PetscScalar* E;
   PetscScalar zero=0.;
//   int ierr;
//   PetscScalar value;
//   
//   Vec RESloc;
//   Vec EPSloc;
//   
//   if(_msh->_nprocs==1) {
//     ierr = VecGetArray(RES,&R);
//     CHKERRQ(ierr);
//     ierr = VecGetArray(EPS,&E);
//     CHKERRQ(ierr);
//   } 
//   else {
//     ierr=VecGhostGetLocalForm(RES,&RESloc);
//     CHKERRQ(ierr);
//     ierr=VecGhostGetLocalForm(EPS,&EPSloc);
//     CHKERRQ(ierr);
//     ierr = VecGetArray(RESloc,&R);
//     CHKERRQ(ierr);
//     ierr = VecGetArray(EPSloc,&E);
//     CHKERRQ(ierr);
//   }
  
  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    unsigned soltype =  _SolType[indexSol];

    int loc_size   = _Eps[indexSol]->local_size();
    int loc_offset_EPS = KKoffset[k][_msh->_iproc];// - KKoffset[0][_msh->_iproc];

    int glob_offset_eps = _msh->MetisOffset[soltype][_msh->_iproc];

    vector <int> index(_msh->own_size[soltype][_msh->_iproc]);
    vector <double> valueEPS(_msh->own_size[soltype][_msh->_iproc]);
    vector <double> valueRES(_msh->own_size[soltype][_msh->_iproc]);
    for(int i=0; i<_msh->own_size[soltype][_msh->_iproc]; i++) {
      index[i]=loc_offset_EPS+i;
    }

   _EPS->get(index,valueEPS);
   _RES->get(index,valueRES);

    for(int i=0; i<_msh->own_size[soltype][_msh->_iproc]; i++) {
      _Eps[indexSol]->set(i+glob_offset_eps,valueEPS[i]);
      if ((*_Bdc[indexSol])(i+glob_offset_eps)>1.1) _Res[indexSol]->set(i+glob_offset_eps,valueRES[i]);
      else _Res[indexSol]->set(i+glob_offset_eps,zero);
    }
    
     
     
     
//      for(int i=0; i<_msh->own_size[soltype][_msh->_iproc]; i++) {
//        _Eps[indexSol]->set(i+glob_offset_eps,E[loc_offset_EPS+i]);
//        if ((*_Bdc[indexSol])(i+glob_offset_eps)>1.1) _Res[indexSol]->set(i+glob_offset_eps,R[loc_offset_EPS+i]);
//        else _Res[indexSol]->set(i+glob_offset_eps,zero);
//     }
    _Res[indexSol]->close();
    _Eps[indexSol]->close();
  }

//   if(_msh->_nprocs==1) {
//     ierr = VecRestoreArray(RES,&R);
//     CHKERRQ(ierr);
//     ierr = VecRestoreArray(EPS,&E);
//     CHKERRQ(ierr);
//   } else {
//     ierr = VecRestoreArray(RESloc,&R);
//     CHKERRQ(ierr);
//     ierr = VecRestoreArray(EPSloc,&E);
//     CHKERRQ(ierr);
//     ierr=VecGhostRestoreLocalForm(RES,&RESloc);
//     CHKERRQ(ierr);
//     ierr=VecGhostRestoreLocalForm(EPS,&EPSloc);
//     CHKERRQ(ierr);
//   }

  for (unsigned k=0; k<_SolPdeIndex.size(); k++) {
    unsigned indexSol=_SolPdeIndex[k];
    _Sol[indexSol]->add(*_Eps[indexSol]);
    _Sol[indexSol]->close();
  }
  
  return 1;*/
}

/**
 * Set _SolOld=_Sol
 **/
// ------------------------------------------------------------------
void Solution::UpdateSolution() {
  for (unsigned i=0; i<_Sol.size(); i++) {
    // Copy the old vector
    if (_SolTmOrder[i]==2) {
      *(_SolOld[i]) = *(_Sol[i]);
    }
  }
}

/**
 * Flag the elements to be refined
 **/
//-------------------------------------------------------------------
void Solution::SetElementRefinement(const unsigned &test) {
 
  _msh->el->InitRefinedToZero();
  
  unsigned grid = _msh->GetGridNumber();
  unsigned nel = _msh->GetElementNumber();
   
  if (test==1) { //refine all next grid elements
    for (unsigned iel=0; iel<nel; iel++) {
      _msh->el->SetRefinedElementIndex(iel,1);
      _msh->el->AddToRefinedElementNumber(1);
      short unsigned elt=_msh->el->GetElementType(iel);
      _msh->el->AddToRefinedElementNumber(1,elt);
    }
  } 
  else if (test==2) { //refine based on the function SetRefinementFlag defined in the main;
    //serial loop
    std::vector<double> X_local;
    std::vector<double> Y_local;
    std::vector<double> Z_local;
    _Sol[0]->localize_to_one(X_local,0);
    _Sol[1]->localize_to_one(Y_local,0);
    _Sol[2]->localize_to_one(Z_local,0);
  
    for (unsigned iel=0; iel<nel; iel+=1) {
      unsigned nve=_msh->el->GetElementDofNumber(iel,0);
      double vtx=0.,vty=0.,vtz=0.;
      for ( unsigned i=0; i<nve; i++) {
        unsigned inode=_msh->el->GetElementVertexIndex(iel,i)-1u;
	unsigned inode_Metis=_msh->GetMetisDof(inode,2);
	vtx+=X_local[inode_Metis];  
	vty+=Y_local[inode_Metis]; 
	vtz+=Z_local[inode_Metis]; 
      }
      vtx/=nve;
      vty/=nve;
      vtz/=nve;
      if (_msh->_SetRefinementFlag(vtx,vty,vtz,_msh->el->GetElementGroup(iel),grid)) {
        _msh->el->SetRefinedElementIndex(iel,1);
        _msh->el->AddToRefinedElementNumber(1);
        short unsigned elt=_msh->el->GetElementType(iel);
        _msh->el->AddToRefinedElementNumber(1,elt);
      }
    }
  } 
  else if (test==3) { //refine all next grid even elements
    for (unsigned iel=0; iel<nel; iel+=2) {
      _msh->el->SetRefinedElementIndex(iel,1);
      _msh->el->AddToRefinedElementNumber(1);
      short unsigned elt=_msh->el->GetElementType(iel);
      _msh->el->AddToRefinedElementNumber(1,elt);
    }
  }
  _msh->el->AllocateChildrenElement(_msh->_ref_index);
}



} //end namespace femus


