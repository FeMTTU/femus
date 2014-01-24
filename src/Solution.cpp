//C++ include
#include <ctime>
#include <fstream>
#include <algorithm>
#include "Solution.hpp"
#include "ElemType.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"
#include "SparseRectangularMatrix.hpp"
using std::cout;
using std::endl;

const unsigned Solution::_END_IND[5]= {0,1,3,4,5};

// ------------ Contructor ------------------------------------------
Solution::Solution(mesh *other_msh){    
 _msh = other_msh;
  for(int i=0;i<5;i++){
    _ProjMatFlag[i]=0;
  }
}

// ------------- Destructor -----------------------------------------
Solution::~Solution() {
  for (unsigned i=0; i<_SolName.size(); i++) {
    delete [] _SolName[i];
  }
}

// ------------- Add a new solution ---------------------------------
void Solution::AddSolutionVector( const char name[], const char order[],
                                const unsigned& tmorder, const bool &PDE_type) {
  unsigned n=_Sol.size();

  _SolType.resize(n+1u);
  _SolName.resize(n+1u);
  _SolTmOrder.resize(n+1u);
  

  _Sol.resize(n+1u);
  _Res.resize(n+1u);
  _Eps.resize(n+1u);
  
  _Bdc.resize(n+1u);
  _ResEpsBdcFlag.resize(n+1u);
  _ResEpsBdcFlag[n]=PDE_type;
  
  
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
    cout<<"error! invalid order entry in AddSolutionVector(...)"<<endl;
    exit(0);
  }
  _SolName[n]=new char [8];
  strcpy(_SolName[n],name);
}

//-------------- Get the solution index base on its name ------------
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

// ------------- Allocate the solution memory -----------------------
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
    
  if(_ResEpsBdcFlag[i]) { //only if the variable is a PDE type
    
    _Res[i] = NumericVector::build().release();
    _Res[i]->init(*_Sol[i]);

    _Eps[i] = NumericVector::build().release();
    _Eps[i]->init(*_Sol[i]);
    
    _Bdc[i] = NumericVector::build().release();
    _Bdc[i]->init(*_Sol[i]);
  }
}

// -------------- Deallocate solution memory ------------------------
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

// ------------- Initialize the coordinate solution X,Y,Z -----------
void Solution::SetCoarseCoordinates( vector < vector < double> > &vt){
  unsigned indexSol;
  
  indexSol=GetIndex("X");
  *_Sol[indexSol]=vt[0];
  
  indexSol=GetIndex("Y");
  *_Sol[indexSol]=vt[1];
  
  indexSol=GetIndex("Z");
  *_Sol[indexSol]=vt[2];
  
}

//--------------------------------------------------------------------------------
int Solution::SumEpsToSol(const vector <unsigned> &MGIndex, const Vec &EPS, const Vec &RES, const vector <vector <unsigned> > &KKoffset ) {

  PetscScalar* R;
  PetscScalar* E;
  PetscScalar zero=0.;
  int ierr;
  PetscScalar value;
  
  Vec RESloc;
  Vec EPSloc;
  
  if(_msh->_nprocs==1) {
    ierr = VecGetArray(RES,&R);
    CHKERRQ(ierr);
    ierr = VecGetArray(EPS,&E);
    CHKERRQ(ierr);
  } 
  else {
    ierr=VecGhostGetLocalForm(RES,&RESloc);
    CHKERRQ(ierr);
    ierr=VecGhostGetLocalForm(EPS,&EPSloc);
    CHKERRQ(ierr);
    ierr = VecGetArray(RESloc,&R);
    CHKERRQ(ierr);
    ierr = VecGetArray(EPSloc,&E);
    CHKERRQ(ierr);
  }
  
  for (unsigned k=0; k<MGIndex.size(); k++) {
    unsigned indexSol=MGIndex[k];
    unsigned soltype =  _SolType[indexSol];

     int loc_size   = _Eps[indexSol]->local_size();
     int loc_offset_EPS = KKoffset[k][_msh->_iproc] - KKoffset[0][_msh->_iproc];

     int glob_offset_eps = _msh->MetisOffset[soltype][_msh->_iproc];

     for(int i=0; i<_msh->own_size[soltype][_msh->_iproc]; i++) {
       _Eps[indexSol]->set(i+glob_offset_eps,E[loc_offset_EPS+i]);
       if ((*_Bdc[indexSol])(i+glob_offset_eps)>1.1) _Res[indexSol]->set(i+glob_offset_eps,R[loc_offset_EPS+i]);
       else _Res[indexSol]->set(i+glob_offset_eps,zero);
    }
    _Res[indexSol]->close();
    _Eps[indexSol]->close();
  }

  if(_msh->_nprocs==1) {
    ierr = VecRestoreArray(RES,&R);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(EPS,&E);
    CHKERRQ(ierr);
  } else {
    ierr = VecRestoreArray(RESloc,&R);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(EPSloc,&E);
    CHKERRQ(ierr);
    ierr=VecGhostRestoreLocalForm(RES,&RESloc);
    CHKERRQ(ierr);
    ierr=VecGhostRestoreLocalForm(EPS,&EPSloc);
    CHKERRQ(ierr);
  }

  for (unsigned k=0; k<MGIndex.size(); k++) {
    unsigned indexSol=MGIndex[k];
    _Sol[indexSol]->add(*_Eps[indexSol]);
    _Sol[indexSol]->close();
  }
  
  return ierr;
}

// ---------------  Set _SolOld=_Sol ----------------------------------
void Solution::UpdateSolution() {
  for (unsigned i=0; i<_Sol.size(); i++) {
    // Copy the old vector
    if (_SolTmOrder[i]==2) {
      *(_SolOld[i]) = *(_Sol[i]);
    }
  }
}

/**
 * This function flags the elements that will be refined
 **/
//--------------------------------------------------------------------------------
void Solution::set_elr(const unsigned &test) {
 
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
  } else if (test==2) { //refine based of the position of the mid point;
    //serial loop
    for (unsigned iel=0; iel<nel; iel+=1) {
      unsigned nve=_msh->el->GetElementDofNumber(iel,0);
      double vtx=0.,vty=0.,vtz=0.;
      for ( unsigned i=0; i<nve; i++) {
        unsigned inode=_msh->el->GetElementVertexIndex(iel,i)-1u;
	unsigned inode_Metis=_msh->GetMetisDof(inode,2);
// 	vtx+=(*_Sol[0])(inode_Metis);  
// 	vty+=(*_Sol[1])(inode_Metis);
// 	vtz+=(*_Sol[2])(inode_Metis);
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
  
//     //parallel loop
//     for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
//       for(int iel=IS_Mts2Gmt_elem_offset[isdom]; iel < IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
//         unsigned kel = IS_Mts2Gmt_elem[iel];
//         unsigned nve=el->GetElementDofNumber(kel,0);
//         double vtx=0.,vty=0.,vtz=0.;
//         for ( unsigned i=0; i<nve; i++) {
//           unsigned inode=el->GetElementVertexIndex(kel,i)-1u;
// 	  unsigned inode_Metis=GetMetisDof(inode,2);
// 	  vtx+=(*_Sol[0])(inode_Metis);  
// 	  vty+=(*_Sol[1])(inode_Metis);
// 	  vtz+=(*_Sol[2])(inode_Metis);
//         }
//         vtx/=nve;
//         vty/=nve;
//         vtz/=nve;
//         if (_SetRefinementFlag(vtx,vty,vtz,el->GetElementGroup(kel),grid)) {
//           el->SetRefinedElementIndex(kel,1);
//           el->AddToRefinedElementNumber(1);
//           short unsigned elt=el->GetElementType(kel);
//           el->AddToRefinedElementNumber(1,elt);
//        }
//      } 
//    }
   
  } else if (test==3) { //refine all next grid even elements
    for (unsigned iel=0; iel<nel; iel+=2) {
      _msh->el->SetRefinedElementIndex(iel,1);
      _msh->el->AddToRefinedElementNumber(1);
      short unsigned elt=_msh->el->GetElementType(iel);
      _msh->el->AddToRefinedElementNumber(1,elt);
    }
  }
  _msh->el->AllocateChildrenElement(_msh->_ref_index);
}

