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
  
  _childElemFlag = false;

//   _isFatherElementRefined = new bool [_nel];
//   memset( _isFatherElementRefined, 0, _nel*sizeof(bool) );

  _kvtel = NULL;
  _kvtelMemory = NULL;
  _nve = NULL;

  _nelr = _nelrt[0] = _nelrt[1] = _nelrt[2] = _nelrt[3] = _nelrt[4] = _nelrt[5] = 0;
}

/**
 * This constructor allocates the memory for the \textit{finer elem}
 * starting from the paramenters of the \textit{coarser elem}
 **/
elem::elem(const elem *elc, const unsigned refindex, const std::vector < double > &coarseAmrVector, const std::vector < double > &coarseElementType) {

  _level = elc->_level + 1;

  _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
  _nel = elc->GetRefinedElementNumber()*refindex; //refined
  _nel += elc->GetElementNumber() - elc->GetRefinedElementNumber(); // + non-refined;

  _elementType = new unsigned short [_nel];

  _isFatherElementRefined = new bool [_nel];

  memset( _isFatherElementRefined, 0, _nel*sizeof(bool) );
  
  _kvert = new unsigned * [_nel];
  _kel = new int * [_nel];

  _kvertSize = 0;
  _kelSize = 0;
  for (unsigned i = 0; i < N_GEOM_ELS; i++) {
    _kvertSize += elc->GetRefinedElementTypeNumber(i) * refindex * NVE[i][2];
    _kelSize += elc->GetRefinedElementTypeNumber(i) * refindex * NFC[i][1];

  }

  for (unsigned iel = 0; iel < elc->GetElementNumber(); iel++ ){
    if( static_cast < short unsigned > ( coarseAmrVector[iel] + 0.25 ) == 0){
       //unsigned type = elc->GetElementType(iel);
       unsigned type = static_cast < short unsigned > ( coarseElementType[iel] + 0.25 );
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
    short unsigned elemt = static_cast < short unsigned > ( coarseElementType[iel] + 0.25 );
    int increment = 1;
    if( static_cast < short unsigned > ( coarseAmrVector[iel] + 0.25 ) == 1){
      increment = NRE[elemt];
    }
    for (unsigned j = 0; j < increment; j++) {
      _kvert[jel+j] = pt_u;
      pt_u += NVE[elemt][2];
      _kel[jel+j] = pt_i;
      pt_i += NFC[ elemt ][1];
    }
    jel += increment;
  }
  
  _childElemFlag = false;

  _kvtel = NULL;
  _kvtelMemory = NULL;
  _nve = NULL;

  _nelr = _nelrt[0] = _nelrt[1] = _nelrt[2] = _nelrt[3] = _nelrt[4] = _nelrt[5] = 0;
}

void elem::ReorderMeshElements( const std::vector < unsigned > &elementMapping , elem *elc){
  //  REORDERING OF  ELT, ELG, ELMAT

  short unsigned *tempElt;
  tempElt = _elementType;
  _elementType = new short unsigned [_nel];
  for(unsigned iel = 0; iel < _nel; iel++){
    _elementType[iel]   = tempElt[ elementMapping[iel] ];
  }
  delete [] tempElt;
  

  if(_level !=0 ){
    bool *tempElRef;
    tempElRef = _isFatherElementRefined;
    _isFatherElementRefined = new bool [_nel];
    for(unsigned iel = 0; iel < _nel; iel++){
      _isFatherElementRefined[iel] = tempElRef[ elementMapping[iel] ];
    } 
    delete [] tempElRef; 
  }
  

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
    //delete [] _isFatherElementRefined;

    delete [] _kvtelMemory;
    delete [] _kvtel;
    delete [] _nve;
    
    _kvtel = NULL;
    _kvtelMemory = NULL;
    _nve = NULL;
    
    if(_childElemFlag){
      delete [] _childElemMemory;
      delete [] _childElem;
    }
  }

void elem::DeleteGroupAndMaterial(){
  delete [] _elementGroup;
  delete [] _elementMaterial;
}

void elem::DeleteElementType(){
  delete [] _elementType;
}

void elem::DeleteElementFather(){
  delete [] _isFatherElementRefined;
}

/**
 * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
 **/
unsigned elem::GetElementDofNumber(const unsigned &iel,const unsigned &type) const {
  return NVE[_elementType[iel]][type];
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
// unsigned elem::GetLocalFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &iedgenode) const {
//   return ig[_elementType[iel]][iface][iedgenode];
// }

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
bool elem::GetIfFatherElementIsRefined(const unsigned &iel) const {
  return _isFatherElementRefined[iel];
}

/**
 * Set the coarse element father
 **/
void elem::SetIfFatherElementIsRefined(const unsigned &iel, const bool &refined) {
  _isFatherElementRefined[iel] = refined;
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


void elem::AllocateChildrenElement(const unsigned &refindex, const std::vector < double > &localizedAmrVector, 
				   const unsigned &offsetStart, const unsigned &offsetEnd){
  if(_childElemFlag){
    delete [] _childElemMemory;
    delete [] _childElem;
  }
  
  _offsetElementStart = offsetStart;
  
  unsigned localNel = offsetEnd - offsetStart;
  unsigned localNelr = 0;
  
  for(unsigned iel = offsetStart; iel < offsetEnd; iel++){
    if ( static_cast < short unsigned > (localizedAmrVector[iel] + 0.25) == 1 ){
      localNelr++;
    }
  }

  _childElemSize = localNelr*refindex + (localNel - localNelr);
  _childElemMemory=new unsigned [_childElemSize];
  
  memset( _childElemMemory, 0, _childElemSize * sizeof(unsigned) );
  
  _childElem = new unsigned* [localNel];

  unsigned *ptr=_childElemMemory;
  for(int i=0; i<localNel; i++){
    _childElem[i]=ptr;
    if( static_cast < short unsigned > ( localizedAmrVector[offsetStart + i] + 0.25) == 1 ) ptr+=refindex;
    else ptr+=1;
  }
  _childElemFlag = true;
  return;
}

void elem::SetChildElement(const unsigned &iel,const unsigned &json, const unsigned &value){
  _childElem[ iel - _offsetElementStart][json]=value;
  return;
}

unsigned elem::GetChildElement(const unsigned &iel,const unsigned &json) const{
  return _childElem[iel - _offsetElementStart ][json];
}

const unsigned elem::GetNVE(const unsigned &elementType, const unsigned &doftype) const{    
  return NVE[elementType][doftype];
}

const unsigned elem::GetNFACENODES(const unsigned &elementType, const unsigned &jface, const unsigned &dof) const{
  return NFACENODES[elementType][jface][dof];
}

const unsigned elem::GetNFC(const unsigned &elementType, const unsigned &type) const{
  return NFC[elementType][type];
}

const unsigned elem::GetIG(const unsigned &elementType, const unsigned &iface, const unsigned &jnode) const{
  return ig[elementType][iface][jnode];
}


} //end namespace femus


