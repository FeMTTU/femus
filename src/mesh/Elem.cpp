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

  _elementDof = new unsigned * [_nel + 1];
  _elementNearFace = new int *[_nel];

  _elementDofSize = _nel*NVE[0][2];
  _elementNearFaceSize = _nel*NFC[0][1];

  _elementDofMemory=new unsigned [_elementDofSize];

  _elementNearFaceMemory=new int [_elementNearFaceSize];
  for (unsigned i=0; i<_elementNearFaceSize; i++)
    _elementNearFaceMemory[i]=-1;

  unsigned *pt_u = _elementDofMemory;
  int *pt_i = _elementNearFaceMemory;

  for (unsigned i=0; i<_nel; i++) {
    _elementDof[i] = pt_u;
    pt_u += NVE[0][2];
    _elementNearFace[i] = pt_i;
    pt_i += NFC[0][1];
  }
  _elementDof[_nel] = pt_u;

  _childElemFlag = false;

  _nelr = _nelrt[0] = _nelrt[1] = _nelrt[2] = _nelrt[3] = _nelrt[4] = _nelrt[5] = 0;
}


void elem::ElementDofSharpAllocation( ){

  _elementDofSize = _nelt[0]*NVE[0][2]+_nelt[1]*NVE[1][2]+
               	    _nelt[2]*NVE[2][2]+_nelt[3]*NVE[3][2]+
		    _nelt[4]*NVE[4][2]+_nelt[5]*NVE[5][2];

  unsigned **tempElementDof = _elementDof;
  unsigned *tempElementDofMemory = _elementDofMemory;

  _elementDof = new unsigned* [_nel+1];
  _elementDofMemory = new unsigned [_elementDofSize];

  unsigned *ptElemDofMem = _elementDofMemory;
  for(unsigned iel = 0; iel<_nel; iel++){
    _elementDof[iel] = ptElemDofMem;
    unsigned ielType = GetElementType(iel);
    for(unsigned j = 0; j<NVE[ielType][2];j++ ){
      _elementDof[iel][j] = tempElementDof[iel][j];
    }
    ptElemDofMem += NVE[ielType][2];
  }
  _elementDof[_nel] = tempElementDofMemory;

  delete [] tempElementDofMemory;
  delete [] tempElementDof;
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

  _fatherElementIsRefined = new bool [_nel];

  memset( _fatherElementIsRefined, 0, _nel*sizeof(bool) );

  _elementDof = new unsigned * [_nel+1];
  _elementNearFace = new int * [_nel];

  _elementDofSize = 0;
  _elementNearFaceSize = 0;
  for (unsigned i = 0; i < N_GEOM_ELS; i++) {
    _elementDofSize += elc->GetRefinedElementTypeNumber(i) * refindex * NVE[i][2];
    _elementNearFaceSize += elc->GetRefinedElementTypeNumber(i) * refindex * NFC[i][1];

  }

  for (unsigned iel = 0; iel < elc->GetElementNumber(); iel++ ){
    if( static_cast < short unsigned > ( coarseAmrVector[iel] + 0.25 ) == 0){
       unsigned type = static_cast < short unsigned > ( coarseElementType[iel] + 0.25 );
       _elementDofSize += NVE[type][2];
       _elementNearFaceSize += NFC[type][1];
    }
  }

  _elementDofMemory = new unsigned [ _elementDofSize ];
  _elementNearFaceMemory = new int [ _elementNearFaceSize ];
  for (unsigned i=0; i < _elementNearFaceSize; i++)
    _elementNearFaceMemory[i] = -1;

  int *pt_i = _elementNearFaceMemory;
  unsigned *pt_u = _elementDofMemory;
  unsigned jel = 0;
  for (unsigned iel = 0; iel<elc->GetElementNumber(); iel++) {
    short unsigned elemt = static_cast < short unsigned > ( coarseElementType[iel] + 0.25 );
    int increment = 1;
    if( static_cast < short unsigned > ( coarseAmrVector[iel] + 0.25 ) == 1){
      increment = NRE[elemt];
    }
    for (unsigned j = 0; j < increment; j++) {
      _elementDof[jel+j] = pt_u;
      pt_u += NVE[elemt][2];
      _elementNearFace[jel+j] = pt_i;
      pt_i += NFC[ elemt ][1];
    }
    jel += increment;
  }
  _elementDof[_nel] = pt_u;


  _childElemFlag = false;

  _nelr = _nelrt[0] = _nelrt[1] = _nelrt[2] = _nelrt[3] = _nelrt[4] = _nelrt[5] = 0;
}

void elem::ReorderMeshElements( const std::vector < unsigned > &elementMapping , elem *elc){
  //  REORDERING OF  ELT, ELG, ELMAT

  short unsigned *tempElt;
  tempElt = _elementType;
  _elementType = new short unsigned [_nel];
  for(unsigned iel = 0; iel < _nel; iel++){
    _elementType[ elementMapping [iel] ]  = tempElt[iel] ;
  }
  delete [] tempElt;


  if(_level !=0 ){
    bool *tempElRef;
    tempElRef = _fatherElementIsRefined;
    _fatherElementIsRefined = new bool [_nel];
    for(unsigned iel = 0; iel < _nel; iel++){
      _fatherElementIsRefined[elementMapping [iel]] = tempElRef[iel];
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
      _elementGroup[elementMapping [iel]]   = tempElg[ iel ];
      _elementMaterial[elementMapping [iel]] = tempElmat[ iel] ;
    }
    delete [] tempElg;
    delete [] tempElmat;
  }

  //  REORDERING OF KEL
  int **tempKel;
  int *tempKelMemory;

  tempKel = _elementNearFace;
  tempKelMemory = _elementNearFaceMemory;

  _elementNearFace = new int * [_nel];
  _elementNearFaceMemory = new int [ _elementNearFaceSize ];

  int *ptKel= _elementNearFaceMemory;

  for(unsigned iel=0; iel<_nel; iel++){
    _elementNearFace[iel] = ptKel;
    ptKel +=  NFC[_elementType[iel]][1];
  }

  for(unsigned iel=0; iel<_nel; iel++){
    for(unsigned iface=0; iface<NFC[_elementType[iel]][1]; iface++){
      _elementNearFace[elementMapping [iel]][iface] = tempKel[iel][iface];
    }
  }
  delete [] tempKelMemory;
  delete [] tempKel;

  //  REORDERING OF ElementDof (ROWS)

  unsigned **tempElementDof = _elementDof;
  unsigned *tempElementDofMemory = _elementDofMemory;

  tempElementDof = _elementDof;
  tempElementDofMemory = _elementDofMemory;

  _elementDof = new unsigned* [_nel + 1];
  _elementDofMemory = new unsigned [ _elementDofSize ];

  unsigned *ptElementDof = _elementDofMemory;

  for(unsigned iel=0; iel<_nel; iel++){
    _elementDof[iel] = ptElementDof;
    ptElementDof +=  NVE[_elementType[iel]][2];
  }
  _elementDof[_nel] = ptElementDof;

  for(unsigned iel=0; iel<_nel; iel++){
    for(unsigned inode=0; inode<NVE[_elementType[iel]][2]; inode++){
      _elementDof[elementMapping [iel]][inode] = tempElementDof[iel][inode];
    }
  }

  delete [] tempElementDof;
  delete [] tempElementDofMemory;

  if(elc){
    for(unsigned i = 0; i < elc->_childElemMemorySize; i++){
      elc->_childElemMemory[i] =  elementMapping[ elc->_childElemMemory[i]];
    }
  }
}


void elem::ReorderMeshNodes( const std::vector < unsigned > &nodeMapping){
  for(unsigned i = 0; i < _elementDofSize; i++){
    _elementDofMemory[i] =  nodeMapping[ _elementDofMemory[i] ];
  }
}


elem::~elem() {
    delete [] _elementDofMemory;
    delete [] _elementDof;
    delete [] _elementNearFaceMemory;
    delete [] _elementNearFace;

    if(_childElemFlag){
      delete [] _childElemMemory;
      delete [] _childElem;

      delete [] _childElemDofMemory;
      delete [] _childElemDofMemoryPointer;
      delete [] _childElemDof;
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
  delete [] _fatherElementIsRefined;
}

void elem::DeleteElementNearVertex(){
  delete [] _elementNearVertexMemory;
  delete [] _elementNearVertex;
  delete [] _elementNearVertexNumber;
}

/**
 * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
 **/
unsigned elem::GetElementDofNumber(const unsigned &iel,const unsigned &type) const {
  return NVE[_elementType[iel]][type];
}


/**
 * Set the local->global node number
 **/
void elem::SetElementDofIndex(const unsigned &iel,const unsigned &inode, const unsigned &value) {
  _elementDof[iel][inode]=value;
}

/**
 * Return the local->global face node number
 **/
unsigned elem::GetFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &inode)const {
  return _elementDof[iel][ig[_elementType[iel]][iface][inode]];
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
  return _elementNearFace[iel][iface];
}

int elem::GetBoundaryIndex(const unsigned &iel,const unsigned &iface) const {
  return -(_elementNearFace[iel][iface]+1);
}


/**
 * Set the global adiacent-to-face element number
 **/
void elem::SetFaceElementIndex(const unsigned &iel,const unsigned &iface, const int &value) {
  _elementNearFace[iel][iface]=value;
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
  return _fatherElementIsRefined[iel];
}

/**
 * Set the coarse element father
 **/
void elem::SetIfFatherElementIsRefined(const unsigned &iel, const bool &refined) {
  _fatherElementIsRefined[iel] = refined;
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

void elem::BuildLocalElementNearVertex(){
  for (unsigned iel = _elementOffset; iel < _elementOffsetP1; iel++) {
    for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
      unsigned inode = GetElementDofIndex(iel,i);
      unsigned inodesize = 1;
      if( _localElementNearVertexMap.find(inode) != _localElementNearVertexMap.end() ){
        inodesize = _localElementNearVertexMap[inode].size();
	inodesize++;
      }
      else{
        _localElementNearVertexMap[inode].reserve(_elementNearVertexNumber[inode]);
      }
      _localElementNearVertexMap[inode].resize(inodesize);
      _localElementNearVertexMap[inode][inodesize-1] = iel;
    }
  }
}



void elem::BuildElementNearVertex() {
  unsigned counter=(_nelt[0]*NVE[0][0]+_nelt[1]*NVE[1][0]+_nelt[2]*NVE[2][0]+
                    _nelt[3]*NVE[3][0]+_nelt[4]*NVE[4][0]+_nelt[5]*NVE[5][0]);

  _elementNearVertex=new unsigned * [_nvt];
  _elementNearVertexNumber= new unsigned[_nvt];

  for (unsigned inode=0; inode<_nvt; inode++) {
    _elementNearVertexNumber[inode]=0;
  }
  for (unsigned iel=0; iel<_nel; iel++) {
    for (unsigned inode=0; inode < GetElementDofNumber(iel,0); inode++) {
      _elementNearVertexNumber[ GetElementDofIndex(iel,inode)]++;
    }
  }

  _elementNearVertexMemory = new unsigned[counter];
  unsigned *pt = _elementNearVertexMemory;
  for (unsigned inode=0; inode<_nvt; inode++) {
    _elementNearVertex[inode] = pt;
    pt+=_elementNearVertexNumber[inode];
    for(unsigned j = 0; j<_elementNearVertexNumber[inode]; j++){
      _elementNearVertex[inode][j] = _nel;
    }
  }

  for (unsigned iel = 0; iel < _nel; iel++) {
    for (unsigned inode = 0; inode < GetElementDofNumber(iel,0); inode++) {
      unsigned irow = GetElementDofIndex(iel,inode);
      unsigned j = 0;
      while ( _nel != _elementNearVertex[irow][j] ) j++;
       _elementNearVertex[irow][j]=iel;
    }
  }
}




/**
 * Return the number of elements which have given node
 **/
unsigned elem::GetElementNearVertexNumber(const unsigned &inode)const {
  return _elementNearVertexNumber[inode];
}

/**
 * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
 **/
unsigned elem::GetElementNearVertex(const unsigned &inode,const unsigned &j)const {
  return _elementNearVertex[inode][j];
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


void elem::AllocateChildrenElement(const unsigned &refindex, Mesh* msh){
  if(_childElemFlag){
    delete [] _childElemMemory;
    delete [] _childElem;

    delete [] _childElemDofMemory;
    delete [] _childElemDofMemoryPointer;
    delete [] _childElemDof;
  }

  unsigned localNel = _elementOffsetP1 - _elementOffset;
  unsigned localNelr = 0;


  _childElemDofMemorySize = 0;
  for(unsigned iel = _elementOffset; iel < _elementOffsetP1; iel++){
    unsigned elementType = msh->GetElementType(iel);
    if ( msh->GetRefinedElementIndex(iel) == 1){
      localNelr++;
      _childElemDofMemorySize += refindex * NVE[elementType][2];
    }
    else{
      _childElemDofMemorySize += NVE[elementType][2];
    }
  }

  _childElemMemorySize = localNelr*refindex + (localNel - localNelr);

  _childElemMemory = new unsigned [_childElemMemorySize];
  memset( _childElemMemory, 0, _childElemMemorySize * sizeof(unsigned) );
  _childElem = new unsigned* [localNel];

  _childElemDofMemory = new unsigned [_childElemDofMemorySize];
  _childElemDofMemoryPointer = new unsigned *[_childElemMemorySize];
  _childElemDof = new unsigned **[localNel];


  unsigned *ptr = _childElemMemory;

  unsigned *ptr1  = _childElemDofMemory;
  unsigned **pptr = _childElemDofMemoryPointer;
  unsigned counter = 0;

  for(unsigned i=0; i < localNel; i++){
    _childElem[i] = ptr;

    unsigned elementType = msh->GetElementType(_elementOffset + i);

    _childElemDof[i] = pptr;
    if ( msh->GetRefinedElementIndex(_elementOffset + i) == 1){
      ptr += refindex;

      pptr += refindex;
      for(unsigned j=0; j<refindex; j++){
	_childElemDofMemoryPointer[counter] = ptr1;
	counter++;
	ptr1 += NVE[elementType][2];
      }
    }
    else {
      ptr += 1;

      pptr += 1;
      _childElemDofMemoryPointer[counter] = ptr1;
      counter++;
      ptr1 += NVE[elementType][2];

    }
  }

  _childElemFlag = true;
  return;
}

void elem::SetChildElementDof(const unsigned &refIndex, Mesh *msh, const elem* elf){
  for(unsigned iel = _elementOffset; iel < _elementOffsetP1; iel++){
    unsigned elementType = msh->GetElementType(iel);
    unsigned endIndex = ( msh->GetRefinedElementIndex(iel) == 1)?refIndex : 1u;
    for(unsigned j=0; j<endIndex; j++){
      unsigned ielf = GetChildElement(iel,j);
      for(unsigned k=0; k<elf->GetElementDofNumber(ielf,2); k++){
	_childElemDof[iel-_elementOffset][j][k] = elf->GetElementDofIndex(ielf,k);
      }
    }
  }
}

void elem::SetChildElement(const unsigned &iel,const unsigned &json, const unsigned &value){
  _childElem[ iel - _elementOffset][json]=value;
  return;
}

unsigned elem::GetChildElement(const unsigned &iel,const unsigned &json) const{
  return _childElem[iel - _elementOffset ][json];
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


