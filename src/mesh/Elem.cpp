/*=========================================================================

 Program: FEMUS
 Module: Elem
 Authors: Eugenio Aulisa, Sara Calandrini, Giacomo Capodaglio

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
  elem::elem(const unsigned& other_nel) {
    _coarseElem = NULL;

    _level = 0;

    _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
    _nel = other_nel;

    _elementType.resize(_nel);
    _elementGroup.resize(_nel);
    _elementMaterial.resize(_nel);
    _elementLevel.resize(_nel, _level);

    _elementDof.resize(_nel, NVE[0][2]);
    _elementNearFace.resize(_nel, NFC[0][1]);

    _nelr = _nelrt[0] = _nelrt[1] = _nelrt[2] = _nelrt[3] = _nelrt[4] = _nelrt[5] = 0;
  }

  void elem::SharpMemoryAllocation() {

    MyVector <unsigned> rowSizeElDof(_nel);
    MyVector <unsigned> rowSizeElNearFace(_nel);
    for(unsigned iel = 0; iel < _nel; iel++) {
      unsigned ielType = GetElementType(iel);
      rowSizeElDof[iel] = NVE[ielType][2];
      rowSizeElNearFace[iel] = NFC[ielType][1];
    }

    MyMatrix <unsigned> tmpElDof(rowSizeElDof);
    MyMatrix <int> tmpElNearFace(rowSizeElNearFace);

    rowSizeElDof.clear();
    rowSizeElNearFace.clear();

    for(unsigned i = tmpElDof.begin(); i < tmpElDof.end(); i++) {
      for(unsigned j = tmpElDof.begin(i); j < tmpElDof.end(i); j++) {
        tmpElDof[i][j] = _elementDof[i][j];
      }
    }

    for(unsigned i = tmpElNearFace.begin(); i < tmpElNearFace.end(); i++) {
      for(unsigned j = tmpElNearFace.begin(i); j < tmpElNearFace.end(i); j++) {
        tmpElNearFace[i][j] = _elementNearFace[i][j];
      }
    }

    _elementDof = tmpElDof;
    _elementNearFace = tmpElNearFace;

    tmpElDof.clear();
    tmpElNearFace.clear();

  }
  /**
   * This constructor allocates the memory for the \textit{finer elem}
   * starting from the parameters of the \textit{coarser elem}
   **/
  elem::elem(elem* elc, const unsigned refindex, const std::vector < double >& coarseAmrVector) {
    _coarseElem = elc;

    _level = elc->_level + 1;

    _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
    _nel = elc->GetRefinedElementNumber() * refindex; //refined
    _nel += elc->GetElementNumber() - elc->GetRefinedElementNumber(); // + non-refined;

    _elementType.resize(_nel);
    _elementGroup.resize(_nel);
    _elementMaterial.resize(_nel);
    _elementLevel.resize(_nel, _level);

    //**************************
    MyVector <unsigned> rowSizeElDof(_nel);
    MyVector <unsigned> rowSizeElNearFace(_nel);
    unsigned jel = 0;
    for(unsigned isdom = 0; isdom < elc->_nprocs; isdom++) {
      elc->_elementType.localize(isdom);
      for(unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
        short unsigned elType = elc->_elementType[iel];
        int increment = 1;
        if(static_cast < short unsigned >(coarseAmrVector[iel] + 0.25) == 1) {
          increment = NRE[elType];
        }
        for(unsigned j = 0; j < increment; j++) {
          rowSizeElDof[jel + j] += NVE[elType][2];
          rowSizeElNearFace[jel + j] += NFC[elType][1];
        }
        jel += increment;
      }
      elc->_elementType.clearLocalized();
    }
    _elementDof = MyMatrix <unsigned> (rowSizeElDof);
    _elementNearFace = MyMatrix <int> (rowSizeElNearFace, -1);

    rowSizeElDof.clear();
    rowSizeElNearFace.clear();

    _nelr = _nelrt[0] = _nelrt[1] = _nelrt[2] = _nelrt[3] = _nelrt[4] = _nelrt[5] = 0;

  }

  void elem::ReorderMeshElements(const std::vector < unsigned >& elementMapping) {

    //BEGIN reordering _elementType
    MyVector <short unsigned> tempElementType;
    tempElementType = _elementType;
    for(unsigned iel = 0; iel < _nel; iel++) {
      _elementType[ elementMapping [iel] ]  = tempElementType[iel] ;
    }
    tempElementType.clear();
    //END reordering _elementType

    //BEGIN reordering _elementLevel
    if(_level != 0) {
      MyVector <short unsigned> tempElementLevel;
      tempElementLevel = _elementLevel;
      for(unsigned iel = 0; iel < _nel; iel++) {
        _elementLevel[elementMapping [iel]] = tempElementLevel[iel];
      }
      tempElementLevel.clear();
    }
    //END reordering _elementCanBeRefined

    //BEGIN reordering _elementGroup and _elementMaterial

    MyVector <short unsigned> tempElementGroup;
    MyVector <short unsigned> tempElementMaterial;
    tempElementGroup = _elementGroup;
    tempElementMaterial = _elementMaterial;

    for(unsigned iel = 0; iel < _nel; iel++) {
      _elementGroup[elementMapping[iel]]   = tempElementGroup[iel];
      _elementMaterial[elementMapping[iel]] = tempElementMaterial[iel] ;
    }
    tempElementGroup.clear();
    tempElementMaterial.clear();

    //END reordering _elementGroup and _elementMaterial


    //BEGIN reordering _elementDof (rows)
    MyVector <unsigned> rowSize(_nel, 0);
    for(unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[elementMapping[i]] = _elementDof.size(i);
    }

    MyMatrix <unsigned> tmpElDof = _elementDof;
    _elementDof = MyMatrix <unsigned>(rowSize, 0);
    for(unsigned i = tmpElDof.begin(); i < tmpElDof.end(); i++) {
      for(unsigned j = tmpElDof.begin(i); j < tmpElDof.end(i); j++) {
        _elementDof[elementMapping [i]][j] = tmpElDof[i][j];
      }
    }
    tmpElDof.clear();
    //END reordering OF _elementDof

    //BEGIN reordering _elementNearFace (rows)
    for(unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[elementMapping[i]] = _elementNearFace.size(i);
    }
    MyMatrix <int> tmpElNearFace = _elementNearFace;
    _elementNearFace = MyMatrix <int>(rowSize, 0);
    for(unsigned i = tmpElNearFace.begin(); i < tmpElNearFace.end(); i++) {
      for(unsigned j = tmpElNearFace.begin(i); j < tmpElNearFace.end(i); j++) {
        _elementNearFace[elementMapping [i]][j] = tmpElNearFace[i][j];
      }
    }
    tmpElNearFace.clear();
    //END reordering _elementNearFace

    //BEGIN reordering _childElementDof (columns) on coarse level
    if(_level != 0) {
      for(unsigned i = _coarseElem->_childElem.begin(); i < _coarseElem->_childElem.end(); i++) {
        for(unsigned j = _coarseElem->_childElem.begin(i); j < _coarseElem->_childElem.end(i); j++) {
          _coarseElem->_childElem[i][j] =  elementMapping[ _coarseElem->_childElem[i][j]];
        }
      }
    }
    //END reordering _childElementDof
  }


  void elem::ReorderMeshNodes(const std::vector < unsigned >& nodeMapping) {
    for(unsigned i = _elementDof.begin(); i < _elementDof.end(); i++) {
      for(unsigned j = _elementDof.begin(i); j < _elementDof.end(i); j++) {
        _elementDof[i][j] =  nodeMapping[_elementDof[i][j]];
      }
    }
  }

  elem::~elem() {}

  void elem::DeleteElementNearVertex() {
    _elementNearVertex.clear();
  }

  /**
   * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
   **/
  unsigned elem::GetElementDofNumber(const unsigned& iel, const unsigned& type) {
    return NVE[_elementType[iel]][type];
  }

  /**
   * Set the local->global node number
   **/
  void elem::SetElementDofIndex(const unsigned& iel, const unsigned& inode, const unsigned& value) {
    _elementDof[iel][inode] = value;
  }

  /**
   * Return the local->global face node index
   **/
  unsigned elem::GetFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& inode) {
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
  void elem::SetNodeNumber(const unsigned& value) {
    _nvt = value;
  }

  /**
   * Return the total number of elements
   **/
  unsigned elem::GetElementNumber(const char* name) const {
    if(!strcmp(name, "All")) {
      return _nel;
    }
    unsigned i;
    i = this->GetIndex(name);
    return _nelt[i];
  }

  /**
   * Add value to the total number of the element
   **/
  void elem::AddToElementNumber(const unsigned& value, const char name[]) {
    unsigned i;
    i = GetIndex(name);
    _nelt[i] += value;
  }
  void elem::AddToElementNumber(const unsigned& value, short unsigned ielt) {
    _nelt[ielt] += value;
  }
  unsigned elem::GetElementFaceNumber(const unsigned& iel, const unsigned& type) {
    return NFC[ _elementType[iel] ][type];
  }

  /**
   * Return the global adiacent-to-face element number
   **/
  int elem::GetFaceElementIndex(const unsigned& iel, const unsigned& iface) {
    return _elementNearFace[iel][iface];
  }

  int elem::GetBoundaryIndex(const unsigned& iel, const unsigned& iface) {
    return  -(GetFaceElementIndex(iel, iface) + 1);
  }

  /**
   * Set the global adiacent-to-face element number
   **/
  void elem::SetFaceElementIndex(const unsigned& iel, const unsigned& iface, const int& value) {
    _elementNearFace[iel][iface] = value;
  }

  /**
   * Return element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  short unsigned elem::GetElementType(const unsigned& iel) {
    return _elementType[iel];
  }

  /**
   * Set element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  void elem::SetElementType(const unsigned& iel, const short unsigned& value) {
    _elementType[iel] = value;
  }

  /**
   * Return element group
   **/
  short unsigned elem::GetElementGroup(const unsigned& iel) {
    return _elementGroup[iel];
  }

  /**
   * Set element group
   **/
  void elem::SetElementGroup(const unsigned& iel, const short unsigned& value) {
    _elementGroup[iel] = value;
  }

  /**
   * Set element Material
  **/
  void elem::SetElementMaterial(const unsigned& iel, const short unsigned& value) {
    _elementMaterial[iel] = value;
  }

  /**
   * Return element material
   **/
  short unsigned elem::GetElementMaterial(const unsigned& iel) {
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
  void elem::SetElementGroupNumber(const unsigned& value) {
    _ngroup = value;
  }

  /**
   * Set the memory storage and initialize nve and kvtel (node->element vectors)
   **/

  void elem::BuildElementNearElement() {
    MyVector < unsigned > rowSize(_elementOffset, 1);
    for(unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
      std::map< unsigned, bool> elements;
      for(unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
        unsigned inode = GetElementDofIndex(iel, i);
        for(unsigned j = _elementNearVertex.begin(inode); j < _elementNearVertex.end(inode); j++) {
          elements[_elementNearVertex[inode][j]] = true;
        }
      }
      rowSize[iel] = elements.size();
    }
    _elementNearElement = MyMatrix <unsigned> (rowSize, UINT_MAX);
    for(unsigned iel = _elementNearElement.begin(); iel < _elementNearElement.end(); iel++) {
      std::map< unsigned, bool> elements;
      _elementNearElement[iel][0] = iel;
      for(unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
        unsigned inode = GetElementDofIndex(iel, i);
        for(unsigned j = _elementNearVertex.begin(inode); j < _elementNearVertex.end(inode); j++) {
          if(_elementNearVertex[inode][j] != iel) {
            elements[_elementNearVertex[inode][j]] = true;
          }
        }
      }
      unsigned j = 1;
      for(std::map<unsigned, bool>::iterator it = elements.begin(); it != elements.end(); it++, j++) {
        _elementNearElement[iel][j] = it->first;
      }
    }
  }

  void elem::BuildElementNearVertex() {
    MyVector <unsigned> rowSize(_nvt, 0);
    for(unsigned iel = 0; iel < _nel; iel++) {
      for(unsigned inode = 0; inode < GetElementDofNumber(iel, 0); inode++) {
        rowSize[GetElementDofIndex(iel, inode)]++;
      }
    }
    _elementNearVertex = MyMatrix <unsigned> (rowSize, _nel);
    for(unsigned iel = 0; iel < _nel; iel++) {
      for(unsigned inode = 0; inode < GetElementDofNumber(iel, 0); inode++) {
        unsigned irow = GetElementDofIndex(iel, inode);
        unsigned j = 0;
        while(_nel != _elementNearVertex[irow][j]) j++;
        _elementNearVertex[irow][j] = iel;
      }
    }
  }

  /**
   * Return the number of elements which have a given vertex
   **/
  unsigned elem::GetElementNearVertexNumber(const unsigned& inode) {
    return _elementNearVertex.size(inode);
  }

  /**
   * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
   **/
  unsigned elem::GetElementNearVertex(const unsigned& inode, const unsigned& j) {
    return _elementNearVertex[inode][j];
  }

  /**
   * return the index 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  unsigned elem::GetIndex(const char name[]) const {
    unsigned index = 0;
    if(!strcmp(name, "Hex")) {
      index = 0;
    }
    else if(!strcmp(name, "Tet")) {
      index = 1;
    }
    else if(!strcmp(name, "Wedge")) {
      index = 2;
    }
    else if(!strcmp(name, "Quad")) {
      index = 3;
    }
    else if(!strcmp(name, "Triangle")) {
      index = 4;
    }
    else if(!strcmp(name, "Line")) {
      index = 5;
    }
    else {
      cout << "error! invalid Element Shape in elem::GetIndex(...)" << endl;
      exit(0);
    }
    return index;
  }

  void elem::AllocateChildrenElement(const unsigned& refindex, Mesh* msh) {
    MyVector <unsigned> rowSize(_elementOffset, 0);
    for(unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[i] = (msh->GetRefinedElementIndex(i) == 1) ? refindex : 1;
    }
    _childElem = MyMatrix <unsigned> (rowSize, 0);

    for(unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      unsigned elementType = msh->GetElementType(i);
      rowSize[i] = (msh->GetRefinedElementIndex(i) == 1) ? refindex * NVE[elementType][2] : NVE[elementType][2];
    }
    _childElemDof = MyMatrix <unsigned> (rowSize, 0);
  }

  void elem::SetChildElementDof(elem* elf) {
    for(unsigned iel = _childElemDof.begin(); iel < _childElemDof.end(); iel++) {
      unsigned nDofs = GetElementDofNumber(iel, 2);
      for(unsigned j = _childElemDof.begin(iel); j < _childElemDof.end(iel); j++) {
        unsigned ielf = GetChildElement(iel, j / nDofs);
        _childElemDof[iel][j] = elf->GetElementDofIndex(ielf, j % nDofs);
      }
    }
  }

  unsigned elem::GetChildElementDof(const unsigned& iel, const unsigned& i0, const unsigned i1)  {
    return _childElemDof[iel][i0 * GetElementDofNumber(iel, 2) + i1];
  }

  void elem::SetChildElement(const unsigned& iel, const unsigned& json, const unsigned& value) {
    _childElem[iel][json] = value;
  }

  unsigned elem::GetChildElement(const unsigned& iel, const unsigned& json) {
    return _childElem[iel][json];
  }

  const unsigned elem::GetNVE(const unsigned& elementType, const unsigned& doftype) const {
    return NVE[elementType][doftype];
  }

  const unsigned elem::GetNFACENODES(const unsigned& elementType, const unsigned& jface, const unsigned& dof) const {
    return NFACENODES[elementType][jface][dof];
  }

  const unsigned elem::GetNFC(const unsigned& elementType, const unsigned& type) const {
    return NFC[elementType][type];
  }

  const unsigned elem::GetIG(const unsigned& elementType, const unsigned& iface, const unsigned& jnode) const {
    return ig[elementType][iface][jnode];
  }

  void elem::ScatterElementDof() {
    _elementDof.scatter(_elementOffset);
  }

  void elem::LocalizeElementDof(const unsigned& jproc) {
    _elementDof.localize(jproc);
  }

  unsigned elem::GetElementDofIndex(const unsigned& iel, const unsigned& inode) {
    return _elementDof[iel][inode];
  };

  void elem::FreeLocalizedElementDof() {
    _elementDof.clearLocalized();
  }

  void elem::ScatterElementNearFace() {
    _elementNearFace.scatter(_elementOffset);
  };

  void elem::LocalizeElementNearFace(const unsigned& jproc) {
    _elementNearFace.localize(jproc);
  }

  void elem::FreeLocalizedElementNearFace() {
    _elementNearFace.clearLocalized();
  }

  void elem::SetLevelInterfaceElement() {

    std::cout << "Mesh level" << _level << std::endl;

    _levelInterfaceElement.resize(_level + 1);
    _levelInterfaceLocalDofs.resize(_level + 1);
    for(unsigned ilevel = 0; ilevel <= _level; ilevel++) {
      unsigned counter = 0;
      for(unsigned i = _elementLevel.begin(); i < _elementLevel.end(); i++) {
        if(ilevel == _elementLevel[i]) {
          for(unsigned j = _elementNearFace.begin(i); j < _elementNearFace.end(i); j++) {
            if(-1 == _elementNearFace[i][j]) {
              counter++;
              break;
            }
          }
        }
      }
      _levelInterfaceElement[ilevel] = MyVector <unsigned> (counter);
      counter = 0;
      for(unsigned i = _elementLevel.begin(); i < _elementLevel.end(); i++) {
        if(ilevel == _elementLevel[i]) {
          for(unsigned j = _elementNearFace.begin(i); j < _elementNearFace.end(i); j++) {
            if(-1 == _elementNearFace[i][j]) {
              _levelInterfaceElement[ilevel][counter] = i;
              counter++;
              break;
            }
          }
        }
      }
      _levelInterfaceElement[ilevel].gather();
      std::cout << "level =" << ilevel << std::endl;
      std::cout << _levelInterfaceElement[ilevel] << std::endl;

      std::vector< unsigned > offset = _levelInterfaceElement[ilevel].getOffset();
      MyVector<unsigned> rowSize(offset, 0);
      for(unsigned i = _levelInterfaceElement[ilevel].begin(); i < _levelInterfaceElement[ilevel].end(); i++) {
        unsigned iel =  _levelInterfaceElement[ilevel][i];
	std::map<unsigned,bool> lDofs;
        for(unsigned jface = _elementNearFace.begin(iel); jface < _elementNearFace.end(iel); jface++) {
//           if(-1 == _elementNearFace[i][jface]) {
// 	    for(unsigned k=0;k < GetNFACENODES(GetElementType(iel), jface, 2);k++){
// 	      unsigned index = GetIG(GetElementType(iel), jface, k);
// 	      lDofs[index] = true;
// 	    }
//           }
        }
        rowSize[i] = lDofs.size();
      }
      std::cout << rowSize <<std::endl;
//       _levelInterfaceLocalDofs[ilevel] = MyMatrix <unsigned>(rowSize,0);
//       std::cout << _levelInterfaceLocalDofs[ilevel] <<std::endl;
    }
  }

} //end namespace femus


