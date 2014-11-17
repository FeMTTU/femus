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

#include "Elem.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <assert.h>


namespace femus {

using std::cout;
using std::endl;

/**
 * This constructor allocates the memory for the \textit{coarsest elem}
 **/
elem::elem(const unsigned &other_nel) {
  nelt[0]=nelt[1]=nelt[2]=nelt[3]=nelt[4]=nelt[5]=0;
  nel=other_nel;

  elt=new unsigned short [nel];
  elg=new unsigned short [nel];
  elmat=new unsigned short [nel];
  elr=new unsigned [nel];
  elf=new unsigned [nel];
  memset(elf,0,nel*sizeof(unsigned));
  nelf=0;

  kvert=new unsigned * [nel];
  kel=new int *[nel];

  kvert_memory=new unsigned [nel*NVE[0][3]];
  kel_memory=new int [nel*NFC[0][1]];
  for (unsigned i=0; i<nel*NFC[0][1]; i++)
    kel_memory[i]=-1;

  unsigned *pt_u=kvert_memory;
  int *pt_i=kel_memory;

  for (unsigned i=0; i<nel; i++) {
    kvert[i]=pt_u;
    pt_u+=NVE[0][3];
    kel[i]=pt_i;
    pt_i+=NFC[0][1];
  }
  _node_region_flag=false;
  _child_elem_flag=false;
}

/**
 * This constructor allocates the memory for the \textit{finer elem}
 * starting from the paramenters of the \textit{coarser elem}
 **/
elem::elem(const elem *elc, const unsigned refindex) {
  nelt[0]=nelt[1]=nelt[2]=nelt[3]=nelt[4]=nelt[5]=0;
//   nel=elc->GetRefinedElementNumber()*REF_INDEX;
  nel=elc->GetRefinedElementNumber()*refindex;

  elt=new unsigned short [nel];
  elg=new unsigned short [nel];
  elmat=new unsigned short [nel];
  elr=new unsigned [nel];
  elf=new unsigned [nel];
  memset(elf,0,nel*sizeof(unsigned));
  nelf=0;


  kvert=new unsigned*[nel];
  kel=new int *[nel];

  unsigned kvert_size=0;
  unsigned kel_size=0;
  for (unsigned i=0; i<6; i++) {
    kvert_size+=elc->GetRefinedElementNumber(i)*NVE[i][3];
    kel_size+=elc->GetRefinedElementNumber(i)*NFC[i][1];
  }
//   kvert_size*=REF_INDEX;
//   kel_size*=REF_INDEX;
  kvert_size*=refindex;
  kel_size*=refindex;

  kvert_memory=new unsigned [kvert_size];
  kel_memory=new int [kel_size];
  for (unsigned i=0; i<kel_size; i++)
    kel_memory[i]=0;

  int *pt_i=kel_memory;
  unsigned *pt_u=kvert_memory;
  unsigned jel=0;
  for (unsigned iel=0; iel<elc->GetElementNumber(); iel++) {
    if ( elc->GetRefinedElementIndex(iel) ) {
      short unsigned elemt=elc->GetElementType(iel);
      for (unsigned j=0; j<NRE[elemt]; j++) {
        kvert[jel+j]=pt_u;
        pt_u+=elc->GetElementDofNumber(iel);

        kel[jel+j]=pt_i;
        pt_i+=elc->GetElementFaceNumber(iel);
      }
      jel+=NRE[elemt];
    }
  }
  _node_region_flag=false;
  _child_elem_flag=false;
}


elem::~elem() {
    delete [] kvert_memory;
    delete [] kvert;
    delete [] kel_memory;
    delete [] kel;
    delete [] elt;
    delete [] elf;
    delete [] elg;
    delete [] elmat;
    delete [] elr;
    delete [] kvtel_memory;
    delete [] kvtel;
    delete [] nve;
    if(_node_region_flag) delete [] _node_region;
    if(_child_elem_flag){
      delete [] _child_elem_memory;
      delete [] _child_elem;
    }
  }
  
/**
 * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
 **/
unsigned elem::GetElementDofNumber(const unsigned &iel,const unsigned &type) const {
  return NVE[elt[iel]][type];
}

unsigned elem::GetElementFaceDofNumber(const unsigned &iel, const unsigned jface, const unsigned &type) const {
  assert(type<3);
  return NFACENODES[elt[iel]][jface][type];
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
  unsigned Dof=(SolType<3)?GetElementVertexIndex(iel,inode)-1u:(nel*inode)+iel;
  return Dof;
}


/**
 * Return the local->global node address
 **/
const unsigned*  elem::GetElementVertexAddress(const unsigned &iel,const unsigned &inode)const {
  return &kvert[iel][inode];
}

/**
 * Set the local->global node number
 **/
void elem::SetElementVertexIndex(const unsigned &iel,const unsigned &inode, const unsigned &value) {
  kvert[iel][inode]=value;
}

/**
 * Return the local->global face node number
 **/
unsigned elem::GetFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &inode)const {
  return kvert[iel][ig[elt[iel]][iface][inode]];
}

/**
 * Return the local(edge/face)->local(surface/volume) node number
 **/
unsigned elem::GetLocalFaceVertexIndex(const unsigned &iel, const unsigned &iface, const unsigned &iedgenode) const {
  return ig[elt[iel]][iface][iedgenode];
}

/**
 * Return the total node number
 **/
unsigned elem::GetNodeNumber()const {
  return nvt;
}

/**
 * Return the total vertex number
 **/
unsigned elem::GetVertexNodeNumber()const {
  return nv0;
}

/**
 * Return the total vertex+midpoint number
 **/
unsigned elem::GetMidpointNodeNumber()const {
  return nv1;
}

/**
 * Return the total vertex+midpoint+facepoint number
 **/
unsigned elem::GetCentralNodeNumber()const {
  return nv2;
}

/**
 * Set the total node number
 **/
void elem::SetNodeNumber(const unsigned &value) {
  nvt=value;
}

/**
 * Set the total vertex number
 **/
void elem::SetVertexNodeNumber(const unsigned &value) {
  nv0=value;
}

/**
 * Set the total vertex+midpoint number
 **/
void elem::SetMidpointNodeNumber(const unsigned &value) {
  nv1=value;
}

/**
 * Set the total vertex+midpoint+facepoint number
 **/
void elem::SetCentralNodeNumber(const unsigned &value) {
  nv2=value;
}

/**
 * Return the total number of the element to refine
 **/
unsigned elem::GetRefinedElementNumber(const char* name) const {
  if (!strcmp(name,"All")) {
    return nelr;
  }
  unsigned i;
  i=GetIndex(name);
  return nelrt[i];
}
unsigned  elem::GetRefinedElementNumber(short unsigned ielt)const {
  return nelrt[ielt];
}

/**
 * Add value to the total number of the refined element
 **/
void elem::AddToRefinedElementNumber(const unsigned &value, const char name[]) {
  if (!strcmp(name,"All")) {
    nelr+=value;
    return;
  }
  unsigned i;
  i=this->GetIndex(name);
  nelrt[i]+=value;
}
void elem::AddToRefinedElementNumber(const unsigned &value, short unsigned ielt) {
  nelrt[ielt]+=value;
}


unsigned elem::GetRefinedElementIndex(const unsigned &iel) const {
  return elr[iel];
}
void elem::SetRefinedElementIndex(const unsigned &iel, const unsigned &value) {
  elr[iel]=value;
}
void elem::InitRefinedToZero() {
  nelr=nelrt[0]=nelrt[1]=nelrt[2]=nelrt[3]=nelrt[4]=nelrt[5]=0;
  for (unsigned i=0; i<nel; i++) elr[i]=0;
}

/**
 * Return the total number of the element
 **/
unsigned elem::GetElementNumber(const char* name) const {
  if (!strcmp(name,"All")) {
    return nel;
  }
  unsigned i;
  i=this->GetIndex(name);
  return nelt[i];
}

/**
 * Add value to the total number of the element
 **/
void elem::AddToElementNumber(const unsigned &value, const char name[]) {
  unsigned i;
  i=GetIndex(name);
  nelt[i]+=value;
}
void elem::AddToElementNumber(const unsigned &value, short unsigned ielt) {
  nelt[ielt]+=value;
}
unsigned elem::GetElementFaceNumber(const unsigned &iel, const unsigned &type)const {
  return NFC[ elt[iel] ][type];
}
unsigned elem::GetElementSquareFaceNumber(const unsigned &iel)const {
  return NFC[ elt[iel] ][0];
}
unsigned elem::GetElementTriangleFaceNumber(const unsigned &iel)const {
  return NFC[ elt[iel] ][1]-NFC[ elt[iel] ][0];
}

/**
 * Return the global adiacent-to-face element number
 **/
int elem::GetFaceElementIndex(const unsigned &iel,const unsigned &iface) const {
  return kel[iel][iface];
}

/**
 * Set the global adiacent-to-face element number
 **/
void elem::SetFaceElementIndex(const unsigned &iel,const unsigned &iface, const int &value) {
  kel[iel][iface]=value;
}


/**
 * Return element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
 **/
short unsigned elem::GetElementType(const unsigned &iel) const {
  return elt[iel];
}

/**
 * Set element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
 **/
void elem::SetElementType(const unsigned &iel, const short unsigned &value) {
  elt[iel]=value;
}


/**
 * Return the coarse element father
 **/
unsigned elem::GetElementFather(const unsigned &iel) const {
  return elf[iel];
}

/**
 * Set the coarse element father
 **/
void elem::SetElementFather(const unsigned &iel, const unsigned &value) {
  elf[iel]=value;
}


/**
 * Set the coarse element father
 **/

void elem::SetNumberElementFather(const unsigned &value) {
  nelf=value;
}

/**
 * Return element group
 **/
short unsigned elem::GetElementGroup(const unsigned &iel) const {
  return elg[iel];
}

/**
 * Set element group
 **/
void elem::SetElementGroup(const unsigned &iel, const short unsigned &value) {
  elg[iel]=value;
}

/**
 * Set element Material
**/
void elem::SetElementMaterial(const unsigned &iel, const short unsigned &value) {
  elmat[iel]=value;
}

/**
 * Return element material
 **/
short unsigned elem::GetElementMaterial(const unsigned &iel) const {
  return elmat[iel];
}

/**
 * Return element group number
 **/
unsigned elem::GetElementGroupNumber() const {
  return ngroup;
}

/**
 * Set element group
 **/
void elem::SetElementGroupNumber(const unsigned &value) {
  ngroup=value;
}


/**
 * Set the memory storage and initialize nve and kvtel (node->element vectors)
 **/
void elem::AllocateVertexElementMemory() {
  unsigned counter=(nelt[0]*NVE[0][0]+nelt[1]*NVE[1][0]+nelt[2]*NVE[2][0]+
                    nelt[3]*NVE[3][0]+nelt[4]*NVE[4][0]+nelt[5]*NVE[5][0]);

  kvtel=new unsigned * [nv0];
  nve= new unsigned[nv0];

  for (unsigned inode=0; inode<nv0; inode++) {
    nve[inode]=0;
  }
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode < GetElementDofNumber(iel,0); inode++) {
      nve[ GetElementVertexIndex(iel,inode)-1u ]++;
    }
  }

  kvtel_memory=new unsigned[counter];
  unsigned *pt= kvtel_memory;
  for (unsigned inode=0; inode<nv0; inode++) {
    kvtel[inode]=pt;
    pt+=nve[inode];
  }
  memset(kvtel_memory, 0, counter*sizeof(unsigned));
}

/**
 * Return the number of elements which have given node
 **/
unsigned elem::GetVertexElementNumber(const unsigned &inode)const {
  return nve[inode];
}

/**
 * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
 **/
unsigned elem::GetVertexElementIndex(const unsigned &inode,const unsigned &jnode)const {
  return kvtel[inode][jnode];
}

/**
 * Return the element address for the given i-node in the j-position with 0<=j<nve(i)
 **/
const unsigned * elem::GetVertexElementAddress(const unsigned &inode,const unsigned &jnode)const {
  return &kvtel[inode][jnode];
}

/**
 * Set the element index for the given i-node in the j-position with 0<=j<nve(i)
 **/
void elem::SetVertexElementIndex(const unsigned &inode,const unsigned &jnode, const unsigned &value) {
  kvtel[inode][jnode]=value;
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
  _node_region_flag=1;
  _node_region=new bool [nvt];
  for (int i=0; i<nvt; i++) _node_region[i]=0;  
  // 0 means Fluid - 1 means Solid  ==> Solid wins on Fluid on Interface nodes
}

bool elem::GetNodeRegion(const unsigned &jnode) const {
  return _node_region[jnode];
}

void  elem::SetNodeRegion(const unsigned &jnode, const bool &value) {
  _node_region[jnode]=value;
}

void elem::AllocateChildrenElement(const unsigned &refindex){
  if(_child_elem_flag){
    delete [] _child_elem_memory;
    delete [] _child_elem;
  }
  
  _child_elem_memory=new unsigned [nelr*refindex];
  _child_elem = new unsigned* [nel];
  
  unsigned *ptr=_child_elem_memory;
  for(int i=0;i<nel;i++){
    _child_elem[i]=ptr;
    if(elr[i]==1) ptr+=refindex;
  }
  _child_elem_flag=true;
  return;
}

void elem::SetChildElement(const unsigned &iel,const unsigned &json, const unsigned &value){
  _child_elem[iel][json]=value;
  return;
}

unsigned elem::GetChildElement(const unsigned &iel,const unsigned &json) const{
  return _child_elem[iel][json];
}


} //end namespace femus


