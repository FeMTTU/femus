/*=========================================================================

 Program: FEMUS
 Module: Elem
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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <assert.h>

#include "Elem.hpp"
#include "GeomElTypeEnum.hpp"

namespace femus {

using std::cout;
using std::endl;

/**
 * This constructor allocates the memory for the \textit{coarsest elem}
 **/
elem::elem(const unsigned &other_nel) {

  _level = 0;

  _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
  _nel = other_nel;

  _elementType = new unsigned short [ _nel ];
  _elementGroup = new unsigned short [ _nel ];
  _elementMaterial = new unsigned short [ _nel ];
  _elr = new unsigned [ _nel ];
  _nelf = 0;

  _kvert = new unsigned * [ _nel ];
  _kel = new int *[ _nel ];

  _kvertSize = _nel*NVE[0][2];
  _kelSize = _nel*NFC[0][1];

  _kvertMemory=new unsigned [_kvertSize];

  _kelMemory=new int [_kelSize];
  for (unsigned i=0; i<_kelSize; i++)
    _kelMemory[i]=-1;

  unsigned *pt_u = _kvertMemory;
  int *pt_i = _kelMemory;

  for (unsigned i=0; i<_nel; i++) {
    _kvert[i] = pt_u;
    pt_u += NVE[0][2];
    _kel[i] = pt_i;
    pt_i += NFC[0][1];
  }
  _nodeRegionFlag = false;
  _childElemFlag = false;

  _isFatherElementRefined = new bool [_nel];
  memset( _isFatherElementRefined, 0, _nel*sizeof(bool) );

  _kvtel = NULL;
  _kvtelMemory = NULL;
  _nve = NULL;
}

/**
 * This constructor allocates the memory for the \textit{finer elem}
 * starting from the paramenters of the \textit{coarser elem}
 **/
elem::elem(const elem *elc, const unsigned refindex) {

  _level = elc->_level + 1;

  _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
  _nel = elc->GetRefinedElementNumber()*refindex; //refined
  _nel += elc->GetElementNumber() - elc->GetRefinedElementNumber(); // + non-refined;

  _elementType = new unsigned short [_nel];
  //_elementGroup = new unsigned short [_nel];
  //_elementMaterial = new unsigned short [_nel];
  _elr = new unsigned [_nel];

  _isFatherElementRefined = new bool [_nel];

  memset( _isFatherElementRefined, 0, _nel*sizeof(bool) );
  _nelf = 0;


  _kvert = new unsigned * [_nel];
  _kel = new int * [_nel];

  _kvertSize = 0;
  _kelSize = 0;
  for (unsigned i = 0; i < N_GEOM_ELS; i++) {
    _kvertSize += elc->GetRefinedElementNumber(i) * refindex * NVE[i][2];
    _kelSize += elc->GetRefinedElementNumber(i) * refindex * NFC[i][1];
  }

  for (unsigned iel = 0; iel < elc->GetElementNumber(); iel++ ){
     if(!elc->GetRefinedElementIndex(iel) ){
       unsigned type = elc->GetElementType(iel);
       _kvertSize += NVE[type][2];
       _kelSize += NFC[type][1];
    }
  }

  _kvertMemory = new unsigned [ _kvertSize ];
  _kelMemory = new int [ _kelSize ];
  for (unsigned i=0; i < _kelSize; i++)
    _kelMemory[i] = -1;

  int *pt_i = _kelMemory;
  unsigned *pt_u = _kvertMemory;
  unsigned jel = 0;
  for (unsigned iel = 0; iel<elc->GetElementNumber(); iel++) {
    short unsigned elemt = elc->GetElementType(iel);
    int increment = 1;
    if ( elc->GetRefinedElementIndex(iel) ) {
      increment = NRE[elemt];
    }
    for (unsigned j = 0; j < increment; j++) {
      _kvert[jel+j] = pt_u;
      pt_u += elc->GetElementDofNumber(iel,2);

      _kel[jel+j] = pt_i;
      pt_i += elc->GetElementFaceNumber(iel);
    }
    jel += increment;
  }
  _nodeRegionFlag = false;
  _childElemFlag = false;

  _kvtel = NULL;
  _kvtelMemory = NULL;
  _nve = NULL;
}

void elem::ReorderMeshElements( const std::vector < unsigned > &elementMapping , elem *elc){
  //  REORDERING OF  ELT, ELG, ELMAT
  short unsigned *tempElt;
  bool *tempElRef;


  tempElt = _elementType;
  tempElRef = _isFatherElementRefined;

  _elementType = new short unsigned [_nel];
  _isFatherElementRefined = new bool [_nel];


  for(unsigned iel = 0; iel < _nel; iel++){
    _elementType[iel]   = tempElt[ elementMapping[iel] ];
    _isFatherElementRefined[iel] = tempElRef[ elementMapping[iel] ];
  }

  delete [] tempElt;
  delete [] tempElRef;

  if( _level == 0){
    short unsigned *tempElg;
    short unsigned *tempElmat;
    tempElg = _elementGroup;
    tempElmat = _elementMaterial;
    _elementGroup = new short unsigned [_nel];
    _elementMaterial = new short unsigned [_nel];
    for(unsigned iel = 0; iel < _nel; iel++){
      _elementGroup[iel]   = tempElg[ elementMapping[iel] ];
      _elementMaterial[iel] = tempElmat[ elementMapping[iel] ];
    }
    delete [] tempElg;
    delete [] tempElmat;
  }



  //  REORDERING OF KEL
  int **tempKel;
  int *tempKelMemory;

  tempKel = _kel;
  tempKelMemory = _kelMemory;

  _kel = new int * [_nel];
  _kelMemory = new int [ _kelSize ];

  int *ptKel= _kelMemory;

  for(unsigned iel=0; iel<_nel; iel++){
    _kel[iel] = ptKel;
    ptKel +=  NFC[_elementType[iel]][1];
  }

  for(unsigned iel=0; iel<_nel; iel++){
    for(unsigned iface=0; iface<NFC[_elementType[iel]][1]; iface++){
      _kel[iel][iface] = tempKel[elementMapping[iel]][iface];
    }
  }
  delete [] tempKelMemory;
  delete [] tempKel;

  //  REORDERING OF KVERT (ROWS)


  unsigned **tempKvert;
  unsigned *tempKvertMemory;

  tempKvert = _kvert;
  tempKvertMemory = _kvertMemory;

  _kvert = new unsigned * [_nel];
  _kvertMemory = new unsigned [ _kvertSize ];

  unsigned *ptKvert= _kvertMemory;

  for(unsigned iel=0; iel<_nel; iel++){
    _kvert[iel] = ptKvert;
    ptKvert +=  NVE[_elementType[iel]][2];
  }

  for(unsigned iel=0; iel<_nel; iel++){
    for(unsigned inode=0; inode<NVE[_elementType[iel]][2]; inode++){
      _kvert[iel][inode] = tempKvert[elementMapping[iel]][inode];
    }
  }


  delete [] tempKvert;
  delete [] tempKvertMemory;


  if(elc){
    std::vector < unsigned > InverseElementMapping(_nel);
    for(unsigned iel=0; iel<_nel; iel++){
      InverseElementMapping[ elementMapping[ iel ] ] = iel;
    }
    unsigned *pt = elc->_childElemMemory;
    for(int i=0; i < elc->_childElemSize; i++){
      unsigned iel = InverseElementMapping[*pt];
      *pt = iel;
      pt++;
    }
  }

}

void elem::ReorderMeshNodes( const std::vector < unsigned > &nodeMapping){
  for(unsigned iel=0; iel<_nel; iel++){
    for(unsigned inode=0; inode<NVE[_elementType[iel]][2]; inode++){
      _kvert[iel][inode] =  nodeMapping[ _kvert[iel][inode] -1u] + 1u;
    }
  }
}


elem::~elem() {
    delete [] _kvertMemory;
    delete [] _kvert;
    delete [] _kelMemory;
    delete [] _kel;
    delete [] _elementType;
    delete [] _isFatherElementRefined;
    delete [] _elr;

    delete [] _kvtelMemory;
    delete [] _kvtel;
    delete [] _nve;
    _kvtel = NULL;
    _kvtelMemory = NULL;
    _nve = NULL;

    if(_nodeRegionFlag) delete [] _nodeRegion;
    if(_childElemFlag){
      delete [] _childElemMemory;
      delete [] _childElem;
    }
  }

void elem::deleteParallelizedQuantities(){
  delete [] _elementGroup;
  delete [] _elementMaterial;
}

/**
 * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
 **/
unsigned elem::GetElementDofNumber(const unsigned &iel,const unsigned &type) const {
  return NVE[_elementType[iel]][type];
}

unsigned elem::GetElementFaceDofNumber(const unsigned &iel, const unsigned jface, const unsigned &type) const {
  assert(type<3);
  return NFACENODES[_elementType[iel]][jface][type];
}

const unsigned elem::GetElementFaceType(const unsigned &kel, const unsigned &jface) const{
  unsigned kelt = GetElementType(kel);
  const unsigned FELT[6][2]= {{3,3},{4,4},{3,4},{5,5},{5,5},{6,6}};
  const unsigned felt = FELT[kelt][jface >= GetElementFaceNumber(kel,0)];
  return felt;
}


/**
 * Return the local to global dof
 **/
unsigned elem::GetMeshDof(const unsigned iel,const unsigned &inode,const unsigned &SolType)const {
  unsigned Dof=(SolType<3)?GetElementVertexIndex(iel,inode)-1u:(_nel*inode)+iel;
  return Dof;
}


/**
 * Return the local->global node address
 **/
const unsigned*  elem::GetElementVertexAddress(const unsigned &iel,const unsigned &inode)const {
  return &_kvert[iel][inode];
}

/**
 * Set the local->global node number
 **/
void elem::SetElementVertexIndex(const unsigned &iel,const unsigned &inode, const unsigned &value) {
  _kvert[iel][inode]=value;
}

/**
 * Return the local->global face node number
 **/
unsigned elem::GetFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &inode)const {
  return _kvert[iel][ig[_elementType[iel]][iface][inode]];
}

/**
 * Return the local(edge/face)->local(surface/volume) node number
 **/
unsigned elem::GetLocalFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &iedgenode) const {
  return ig[_elementType[iel]][iface][iedgenode];
}

/**
 * Return the total node number
 **/
unsigned elem::GetNodeNumber()const {
  return _nvt;
}

/**
 * Set the total node number
 **/
void elem::SetNodeNumber(const unsigned &value) {
  _nvt=value;
}

/**
 * Return the total number of the element to refine
 **/
unsigned elem::GetRefinedElementNumber(const char* name) const {
  if (!strcmp(name,"All")) {
    return _nelr;
  }
  unsigned i;
  i=GetIndex(name);
  return _nelrt[i];
}
unsigned  elem::GetRefinedElementNumber(short unsigned ielt)const {
  return _nelrt[ielt];
}

/**
 * Add value to the total number of the refined element
 **/
void elem::AddToRefinedElementNumber(const unsigned &value, const char name[]) {
  if (!strcmp(name,"All")) {
    _nelr+=value;
    return;
  }
  unsigned i;
  i=this->GetIndex(name);
  _nelrt[i]+=value;
}
void elem::AddToRefinedElementNumber(const unsigned &value, short unsigned ielt) {
  _nelrt[ielt]+=value;
}


unsigned elem::GetRefinedElementIndex(const unsigned &iel) const {
  return _elr[iel];
}
void elem::SetRefinedElementIndex(const unsigned &iel, const unsigned &value) {
  _elr[iel]=value;
}
void elem::InitRefinedToZero() {
  _nelr=_nelrt[0]=_nelrt[1]=_nelrt[2]=_nelrt[3]=_nelrt[4]=_nelrt[5]=0;
  for (unsigned i=0; i<_nel; i++) _elr[i]=0;
}

/**
 * Return the total number of the element
 **/
unsigned elem::GetElementNumber(const char* name) const {
  if (!strcmp(name,"All")) {
    return _nel;
  }
  unsigned i;
  i=this->GetIndex(name);
  return _nelt[i];
}

/**
 * Add value to the total number of the element
 **/
void elem::AddToElementNumber(const unsigned &value, const char name[]) {
  unsigned i;
  i=GetIndex(name);
  _nelt[i]+=value;
}
void elem::AddToElementNumber(const unsigned &value, short unsigned ielt) {
  _nelt[ielt]+=value;
}
unsigned elem::GetElementFaceNumber(const unsigned &iel, const unsigned &type)const {
  return NFC[ _elementType[iel] ][type];
}
unsigned elem::GetElementSquareFaceNumber(const unsigned &iel)const {
  return NFC[ _elementType[iel] ][0];
}
unsigned elem::GetElementTriangleFaceNumber(const unsigned &iel)const {
  return NFC[ _elementType[iel] ][1]-NFC[ _elementType[iel] ][0];
}

/**
 * Return the global adiacent-to-face element number
 **/
int elem::GetFaceElementIndex(const unsigned &iel,const unsigned &iface) const {
  return _kel[iel][iface];
}

int elem::GetBoundaryIndex(const unsigned &iel,const unsigned &iface) const {
  return -(_kel[iel][iface]+1);
}


/**
 * Set the global adiacent-to-face element number
 **/
void elem::SetFaceElementIndex(const unsigned &iel,const unsigned &iface, const int &value) {
  _kel[iel][iface]=value;
}


/**
 * Return element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
 **/
short unsigned elem::GetElementType(const unsigned &iel) const {
  return _elementType[iel];
}

/**
 * Set element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
 **/
void elem::SetElementType(const unsigned &iel, const short unsigned &value) {
  _elementType[iel]=value;
}




/**
 * Return if the coarse element father has been refined
 **/
bool elem::IsFatherRefined(const unsigned &iel) const {
  return _isFatherElementRefined[iel];
}

/**
 * Set the coarse element father
 **/
void elem::SetIfFatherIsRefined(const unsigned &iel, const bool &refined) {
  _isFatherElementRefined[iel] = refined;
}


/**
 * Set the coarse element father
 **/

void elem::SetNumberElementFather(const unsigned &value) {
  _nelf=value;
}

/**
 * Return element group
 **/
short unsigned elem::GetElementGroup(const unsigned &iel) const {
  return _elementGroup[iel];
}

/**
 * Set element group
 **/
void elem::SetElementGroup(const unsigned &iel, const short unsigned &value) {
  _elementGroup[iel]=value;
}

/**
 * Set element Material
**/
void elem::SetElementMaterial(const unsigned &iel, const short unsigned &value) {
  _elementMaterial[iel]=value;
}

/**
 * Return element material
 **/
short unsigned elem::GetElementMaterial(const unsigned &iel) const {
  return _elementMaterial[iel];
}

/**
 * Return element group number
 **/
unsigned elem::GetElementGroupNumber() const {
  return _ngroup;
}

/**
 * Set element group
 **/
void elem::SetElementGroupNumber(const unsigned &value) {
  _ngroup=value;
}


/**
 * Set the memory storage and initialize nve and kvtel (node->element vectors)
 **/
void elem::AllocateVertexElementMemory() {
  unsigned counter=(_nelt[0]*NVE[0][0]+_nelt[1]*NVE[1][0]+_nelt[2]*NVE[2][0]+
                    _nelt[3]*NVE[3][0]+_nelt[4]*NVE[4][0]+_nelt[5]*NVE[5][0]);


  if( _kvtel != NULL) delete [] _kvtel;
  if( _kvtelMemory != NULL) delete [] _kvtelMemory;
  if( _nve != NULL ) delete [] _nve;

  _kvtel=new unsigned * [_nvt];
  _nve= new unsigned[_nvt];

  for (unsigned inode=0; inode<_nvt; inode++) {
    _nve[inode]=0;
  }
  for (unsigned iel=0; iel<_nel; iel++) {
    for (unsigned inode=0; inode < GetElementDofNumber(iel,0); inode++) {
      _nve[ GetElementVertexIndex(iel,inode)-1u ]++;
    }
  }

  _kvtelMemory=new unsigned[counter];
  unsigned *pt= _kvtelMemory;
  for (unsigned inode=0; inode<_nvt; inode++) {
    _kvtel[inode]=pt;
    pt+=_nve[inode];
  }
  memset(_kvtelMemory, 0, counter*sizeof(unsigned));
}

/**
 * Return the number of elements which have given node
 **/
unsigned elem::GetVertexElementNumber(const unsigned &inode)const {
  return _nve[inode];
}

/**
 * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
 **/
unsigned elem::GetVertexElementIndex(const unsigned &inode,const unsigned &jnode)const {
  return _kvtel[inode][jnode];
}

/**
 * Return the element address for the given i-node in the j-position with 0<=j<nve(i)
 **/
const unsigned * elem::GetVertexElementAddress(const unsigned &inode,const unsigned &jnode)const {
  return &_kvtel[inode][jnode];
}

/**
 * Set the element index for the given i-node in the j-position with 0<=j<nve(i)
 **/
void elem::SetVertexElementIndex(const unsigned &inode,const unsigned &jnode, const unsigned &value) {
  _kvtel[inode][jnode]=value;
}

/**
 * return the index 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
 **/
unsigned elem::GetIndex(const char name[]) const {
  unsigned index=0;
  if (!strcmp(name,"Hex")) {
    index=0;
  } else if (!strcmp(name,"Tet")) {
    index=1;
  } else if (!strcmp(name,"Wedge")) {
    index=2;
  } else if (!strcmp(name,"Quad")) {
    index=3;
  } else if (!strcmp(name,"Triangle")) {
    index=4;
  } else if (!strcmp(name,"Line")) {
    index=5;
  } else {
    cout<<"error! invalid Element Shape in elem::GetIndex(...)"<<endl;
    exit(0);
  }
  return index;
}


void elem::AllocateNodeRegion() {
  _nodeRegionFlag=1;
  _nodeRegion=new bool [_nvt];
  for (int i=0; i<_nvt; i++) _nodeRegion[i]=0;
  // 0 means Fluid - 1 means Solid  ==> Solid wins on Fluid on Interface nodes
}

bool elem::GetNodeRegion(const unsigned &jnode) const {
  return _nodeRegion[jnode];
}

void  elem::SetNodeRegion(const unsigned &jnode, const bool &value) {
  _nodeRegion[jnode]=value;
}

void elem::AllocateChildrenElement(const unsigned &refindex){
  if(_childElemFlag){
    delete [] _childElemMemory;
    delete [] _childElem;
  }

  _childElemSize = _nelr*refindex+(_nel-_nelr);
  _childElemMemory=new unsigned [_childElemSize];
  _childElem = new unsigned* [_nel];

  unsigned *ptr=_childElemMemory;
  for(int i=0;i<_nel;i++){
    _childElem[i]=ptr;
    if(_elr[i]==1) ptr+=refindex;
    else ptr+=1;
  }
  _childElemFlag=true;
  return;
}

void elem::SetChildElement(const unsigned &iel,const unsigned &json, const unsigned &value){
  _childElem[iel][json]=value;
  return;
}

unsigned elem::GetChildElement(const unsigned &iel,const unsigned &json) const{
  return _childElem[iel][json];
}


} //end namespace femus


