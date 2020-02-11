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

namespace femus
{

  using std::cout;
  using std::endl;

  /**
   * This constructor allocates the memory for the \textit{coarsest elem}
   **/
  elem::elem(const unsigned& other_nel)
  {
    _coarseElem = NULL;

    _level = 0;

    _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
    _nel = other_nel;

    _elementType.resize(_nel);
    _elementGroup.resize(_nel);
    _elementMaterial.resize(_nel);
    _elementLevel.resize(_nel, _level);

    _elementDof.resize(_nel, NVE[0][2], UINT_MAX);
    _elementNearFace.resize(_nel, NFC[0][1], -1);

}

  void elem::ShrinkToFit()
  {

    _elementDof.shrinkToFit(UINT_MAX);

    MyVector <unsigned> rowSize(_nel);
    for (unsigned iel = 0; iel < _nel; iel++) {
      unsigned ielType = GetElementType(iel);
      rowSize[iel] = NFC[ielType][1];
    }
    _elementNearFace.shrinkToFit(rowSize);

  }
  /**
   * This constructor allocates the memory for the \textit{finer elem}
   * starting from the parameters of the \textit{coarser elem}
   **/
  elem::elem(elem* elc, const unsigned refindex, const std::vector < double >& coarseAmrVector)
  {
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
    for (unsigned isdom = 0; isdom < elc->_nprocs; isdom++) {
      elc->_elementType.broadcast(isdom);
      for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
        short unsigned elType = elc->_elementType[iel];
        int increment = 1;
        if (static_cast < short unsigned >(coarseAmrVector[iel] + 0.25) == 1) {
          increment = NRE[elType];
        }
        for (unsigned j = 0; j < increment; j++) {
          rowSizeElDof[jel + j] += NVE[elType][2];
          rowSizeElNearFace[jel + j] += NFC[elType][1];
        }
        jel += increment;
      }
      elc->_elementType.clearBroadcast();
    }
    _elementDof = MyMatrix <unsigned> (rowSizeElDof);
    _elementNearFace = MyMatrix <int> (rowSizeElNearFace, -1);

    rowSizeElDof.clear();
    rowSizeElNearFace.clear();

  }

  void elem::ReorderMeshElements(const std::vector < unsigned >& elementMapping)
  {

    //BEGIN reordering _elementType
    MyVector <short unsigned> tempElementType;
    tempElementType = _elementType;
    for (unsigned iel = 0; iel < _nel; iel++) {
      _elementType[ elementMapping [iel] ]  = tempElementType[iel] ;
    }
    tempElementType.clear();
    //END reordering _elementType

    //BEGIN reordering _elementLevel
    if (_level != 0) {
      MyVector <short unsigned> tempElementLevel;
      tempElementLevel = _elementLevel;
      for (unsigned iel = 0; iel < _nel; iel++) {
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

    for (unsigned iel = 0; iel < _nel; iel++) {
      _elementGroup[elementMapping[iel]]   = tempElementGroup[iel];
      _elementMaterial[elementMapping[iel]] = tempElementMaterial[iel] ;
    }
    tempElementGroup.clear();
    tempElementMaterial.clear();

    //END reordering _elementGroup and _elementMaterial


    //BEGIN reordering _elementDof (rows)
    MyVector <unsigned> rowSize(_nel, 0);
    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[elementMapping[i]] = _elementDof.size(i);
    }

    MyMatrix <unsigned> tmpElDof = _elementDof;
    _elementDof = MyMatrix <unsigned>(rowSize, 0);
    for (unsigned i = tmpElDof.begin(); i < tmpElDof.end(); i++) {
      for (unsigned j = tmpElDof.begin(i); j < tmpElDof.end(i); j++) {
        _elementDof[elementMapping [i]][j] = tmpElDof[i][j];
      }
    }
    tmpElDof.clear();
    //END reordering OF _elementDof

    //BEGIN reordering _elementNearFace (rows)
    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[elementMapping[i]] = _elementNearFace.size(i);
    }
    MyMatrix <int> tmpElNearFace = _elementNearFace;
    _elementNearFace = MyMatrix <int>(rowSize, 0);
    for (unsigned i = tmpElNearFace.begin(); i < tmpElNearFace.end(); i++) {
      for (unsigned j = tmpElNearFace.begin(i); j < tmpElNearFace.end(i); j++) {
        _elementNearFace[elementMapping [i]][j] = tmpElNearFace[i][j];
      }
    }
    tmpElNearFace.clear();
    //END reordering _elementNearFace

    //BEGIN reordering _childElementDof (columns) on coarse level
    if (_level != 0) {
      for (unsigned i = _coarseElem->_childElem.begin(); i < _coarseElem->_childElem.end(); i++) {
        for (unsigned j = _coarseElem->_childElem.begin(i); j < _coarseElem->_childElem.end(i); j++) {
          _coarseElem->_childElem[i][j] =  elementMapping[ _coarseElem->_childElem[i][j]];
        }
      }
    }
    //END reordering _childElementDof
  }


  void elem::ReorderMeshNodes(const std::vector < unsigned >& nodeMapping)
  {
    for (unsigned i = _elementDof.begin(); i < _elementDof.end(); i++) {
      for (unsigned j = _elementDof.begin(i); j < _elementDof.end(i); j++) {
        _elementDof[i][j] =  nodeMapping[_elementDof[i][j]];
      }
    }
  }

  elem::~elem()
  {
  }

  void elem::DeleteElementNearVertex()
  {
    _elementNearVertex.clear();
  }

  /**
   * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
   **/
  unsigned elem::GetElementDofNumber(const unsigned& iel, const unsigned& type)
  {
    return NVE[_elementType[iel]][type];
  }

  /**
   * Set the local->global node number
   **/
  void elem::SetElementDofIndex(const unsigned& iel, const unsigned& inode, const unsigned& value)
  {
    _elementDof[iel][inode] = value;
  }

  /**
   * Return the local->global face node index
   **/
  unsigned elem::GetFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& inode)
  {
    return _elementDof[iel][ig[_elementType[iel]][iface][inode]];
  }

  /**
   * Return the total node number
   **/
  unsigned elem::GetNodeNumber()const
  {
    return _nvt;
  }

  /**
   * Set the total node number
   **/
  void elem::SetNodeNumber(const unsigned& value)
  {
    _nvt = value;
  }

  /**
   * Return the total number of elements
   **/
  unsigned elem::GetElementNumber(const char* name) const
  {
    if (!strcmp(name, "All")) {
      return _nel;
    }
    unsigned i;
    i = this->GetIndex(name);
    return _nelt[i];
  }

  /**
   * Add value to the total number of the element
   **/
  void elem::AddToElementNumber(const unsigned& value, const char name[])
  {
    unsigned i;
    i = GetIndex(name);
    _nelt[i] += value;
  }
  void elem::AddToElementNumber(const unsigned& value, short unsigned ielt)
  {
    _nelt[ielt] += value;
  }
  unsigned elem::GetElementFaceNumber(const unsigned& iel, const unsigned& type)
  {
    return NFC[ _elementType[iel] ][type];
  }

  /**
   * Return the global adiacent-to-face element number
   **/
  int elem::GetFaceElementIndex(const unsigned& iel, const unsigned& iface)
  {
    return _elementNearFace[iel][iface];
  }

  int elem::GetBoundaryIndex(const unsigned& iel, const unsigned& iface)
  {
    return  -(GetFaceElementIndex(iel, iface) + 1);
  }

  /**
   * Set the global adiacent-to-face element number
   **/
  void elem::SetFaceElementIndex(const unsigned& iel, const unsigned& iface, const int& value)
  {
    _elementNearFace[iel][iface] = value;
  }

  /**
   * Return element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  short unsigned elem::GetElementType(const unsigned& iel)
  {
    return _elementType[iel];
  }

  /**
   * Set element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  void elem::SetElementType(const unsigned& iel, const short unsigned& value)
  {
    _elementType[iel] = value;
  }

  /**
   * Return element group
   **/
  short unsigned elem::GetElementGroup(const unsigned& iel)
  {
    return _elementGroup[iel];
  }

  /**
   * Set element group
   **/
  void elem::SetElementGroup(const unsigned& iel, const short unsigned& value)
  {
    _elementGroup[iel] = value;
  }

  /**
   * Set element Material
  **/
  void elem::SetElementMaterial(const unsigned& iel, const short unsigned& value)
  {
    _elementMaterial[iel] = value;
  }

  /**
   * Return element material
   **/
  short unsigned elem::GetElementMaterial(const unsigned& iel)
  {
    return _elementMaterial[iel];
  }

  /**
   * Return element group number
   **/
  unsigned elem::GetElementGroupNumber() const
  {
    return _ngroup;
  }

  /**
   * Set element group
   **/
  void elem::SetElementGroupNumber(const unsigned& value)
  {
    _ngroup = value;
  }

  /**
   * Set the memory storage and initialize nve and kvtel (node->element vectors)
   **/

  void elem::BuildElementNearElement()
  {
    MyVector < unsigned > rowSize(_elementOffset, 1);
    for (unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
      std::map< unsigned, bool> elements;
      for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
        unsigned inode = GetElementDofIndex(iel, i);
        for (unsigned j = _elementNearVertex.begin(inode); j < _elementNearVertex.end(inode); j++) {
          elements[_elementNearVertex[inode][j]] = true;
        }
      }
      rowSize[iel] = elements.size();
    }
    _elementNearElement = MyMatrix <unsigned> (rowSize, UINT_MAX);
    for (unsigned iel = _elementNearElement.begin(); iel < _elementNearElement.end(); iel++) {
      std::map< unsigned, bool> elements;
      _elementNearElement[iel][0] = iel;
      for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
        unsigned inode = GetElementDofIndex(iel, i);
        for (unsigned j = _elementNearVertex.begin(inode); j < _elementNearVertex.end(inode); j++) {
          if (_elementNearVertex[inode][j] != iel) {
            elements[_elementNearVertex[inode][j]] = true;
          }
        }
      }
      unsigned j = 1;
      for (std::map<unsigned, bool>::iterator it = elements.begin(); it != elements.end(); it++, j++) {
        _elementNearElement[iel][j] = it->first;
      }
    }
  }

  void elem::BuildElementNearVertex()
  {
    MyVector <unsigned> rowSize(_nvt, 0);
    for (unsigned iel = 0; iel < _nel; iel++) {
      for (unsigned inode = 0; inode < GetElementDofNumber(iel, 0); inode++) {
        rowSize[GetElementDofIndex(iel, inode)]++;
      }
    }
    _elementNearVertex = MyMatrix <unsigned> (rowSize, _nel);
    for (unsigned iel = 0; iel < _nel; iel++) {
      for (unsigned inode = 0; inode < GetElementDofNumber(iel, 0); inode++) {
        unsigned irow = GetElementDofIndex(iel, inode);
        unsigned j = 0;
        while (_nel != _elementNearVertex[irow][j]) j++;
        _elementNearVertex[irow][j] = iel;
      }
    }
  }

  /**
   * Return the number of elements which have a given vertex
   **/
  unsigned elem::GetElementNearVertexNumber(const unsigned& inode)
  {
    return _elementNearVertex.size(inode);
  }

  /**
   * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
   **/
  unsigned elem::GetElementNearVertex(const unsigned& inode, const unsigned& j)
  {
    return _elementNearVertex[inode][j];
  }

  /**
   * return the index 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  unsigned elem::GetIndex(const char name[]) const
  {
    unsigned index = 0;
    if (!strcmp(name, "Hex")) {
      index = 0;
    }
    else if (!strcmp(name, "Tet")) {
      index = 1;
    }
    else if (!strcmp(name, "Wedge")) {
      index = 2;
    }
    else if (!strcmp(name, "Quad")) {
      index = 3;
    }
    else if (!strcmp(name, "Triangle")) {
      index = 4;
    }
    else if (!strcmp(name, "Line")) {
      index = 5;
    }
    else {
      cout << "error! invalid Element Shape in elem::GetIndex(...)" << endl;
      exit(0);
    }
    return index;
  }

  void elem::AllocateChildrenElement(const unsigned& refindex, Mesh* msh)
  {
    MyVector <unsigned> rowSize(_elementOffset, 0);
    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[i] = (msh->GetRefinedElementIndex(i) == 1) ? refindex : 1;
    }
    _childElem = MyMatrix <unsigned> (rowSize, 0);

    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      unsigned elementType = msh->GetElementType(i);
      rowSize[i] = (msh->GetRefinedElementIndex(i) == 1) ? refindex * NVE[elementType][2] : NVE[elementType][2];
    }
    _childElemDof = MyMatrix <unsigned> (rowSize, 0);
  }

  void elem::SetChildElementDof(elem* elf)
  {
    for (unsigned iel = _childElemDof.begin(); iel < _childElemDof.end(); iel++) {
      unsigned nDofs = GetElementDofNumber(iel, 2);
      for (unsigned j = _childElemDof.begin(iel); j < _childElemDof.end(iel); j++) {
        unsigned ielf = GetChildElement(iel, j / nDofs);
        _childElemDof[iel][j] = elf->GetElementDofIndex(ielf, j % nDofs);
      }
    }
  }

  unsigned elem::GetChildElementDof(const unsigned& iel, const unsigned& i0, const unsigned i1)
  {
    return _childElemDof[iel][i0 * GetElementDofNumber(iel, 2) + i1];
  }

  void elem::SetChildElement(const unsigned& iel, const unsigned& json, const unsigned& value)
  {
    _childElem[iel][json] = value;
  }

  unsigned elem::GetChildElement(const unsigned& iel, const unsigned& json)
  {
    return _childElem[iel][json];
  }

  const unsigned elem::GetNVE(const unsigned& elementType, const unsigned& doftype) const
  {
    return NVE[elementType][doftype];
  }

  const unsigned elem::GetNFACENODES(const unsigned& elementType, const unsigned& jface, const unsigned& dof) const
  {
    return NFACENODES[elementType][jface][dof];
  }

  const unsigned elem::GetNFC(const unsigned& elementType, const unsigned& type) const
  {
    return NFC[elementType][type];
  }

  const unsigned elem::GetIG(const unsigned& elementType, const unsigned& iface, const unsigned& jnode) const
  {
    return ig[elementType][iface][jnode];
  }

  void elem::ScatterElementDof()
  {
    _elementDof.scatter(_elementOffset);
  }

  void elem::LocalizeElementDof(const unsigned& jproc)
  {
    _elementDof.broadcast(jproc);
  }

  unsigned elem::GetElementDofIndex(const unsigned& iel, const unsigned& inode)
  {
    return _elementDof[iel][inode];
  };

  void elem::FreeLocalizedElementDof()
  {
    _elementDof.clearBroadcast();
  }

  void elem::ScatterElementNearFace()
  {
    _elementNearFace.scatter(_elementOffset);
  };

  void elem::LocalizeElementNearFace(const unsigned& jproc)
  {
    _elementNearFace.broadcast(jproc);
  }

  void elem::FreeLocalizedElementNearFace()
  {
    _elementNearFace.clearBroadcast();
  }

  void elem::GetAMRRestriction(Mesh *msh)
  {

    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > &restriction = msh->GetAmrRestrictionMap();
    restriction.resize(3);

    std::vector < std::map < unsigned, bool > > &interfaceSolidMark = msh->GetAmrSolidMark();
    interfaceSolidMark.resize(3);

    std::vector < MyVector<unsigned> > interfaceElement;
    std::vector < MyMatrix<unsigned> > interfaceLocalDof;
    std::vector < std::vector < MyMatrix<unsigned> > > interfaceDof;
    std::vector < std::vector < MyMatrix<unsigned> > > levelInterfaceSolidMark;
    std::vector < std::vector < MyMatrix< double > > > interfaceNodeCoordinates;

    interfaceElement.resize(_level + 1);
    interfaceLocalDof.resize(_level + 1);
    interfaceDof.resize(3);
    levelInterfaceSolidMark.resize(3);
    for (unsigned i = 0; i < 3; i++) {
      interfaceDof[i].resize(_level + 1);
      levelInterfaceSolidMark[i].resize(_level + 1);
    }
    interfaceNodeCoordinates.resize(_level + 1);
    unsigned dim = msh->GetDimension();

    for (unsigned ilevel = 0; ilevel <= _level; ilevel++) {
      //BEGIN interface element search
      interfaceElement[ilevel] = MyVector <unsigned> (_elementOwned);
      unsigned counter = 0;
      for (unsigned i = _elementLevel.begin(); i < _elementLevel.end(); i++) {
        if (ilevel == _elementLevel[i]) {
          for (unsigned j = _elementNearFace.begin(i); j < _elementNearFace.end(i); j++) {
            if (-1 == _elementNearFace[i][j]) {
              interfaceElement[ilevel][counter] = i;
              counter++;
              break;
            }
          }
        }
      }
      interfaceElement[ilevel].resize(counter);
      interfaceElement[ilevel].stack();
      //END interface element search

      //BEGIN interface node search
      std::vector< unsigned > offset = interfaceElement[ilevel].getOffset();
      interfaceLocalDof[ilevel] = MyMatrix <unsigned>(offset, NVE[0][2], UINT_MAX);
      for (unsigned i = interfaceElement[ilevel].begin(); i < interfaceElement[ilevel].end(); i++) {
        unsigned iel =  interfaceElement[ilevel][i];
        std::map <unsigned, bool> ldofs;
        for (unsigned jface = _elementNearFace.begin(iel); jface < _elementNearFace.end(iel); jface++) {
          if (-1 == _elementNearFace[iel][jface]) {
            for (unsigned k = 0; k < GetNFACENODES(GetElementType(iel), jface, 2); k++) {
              unsigned index = GetIG(GetElementType(iel), jface, k);
              ldofs[index] = true;
            }
          }
        }
        unsigned j = 0;
        for (std::map<unsigned, bool>::iterator it = ldofs.begin(); it != ldofs.end(); it++) {
          interfaceLocalDof[ilevel][i][j] = it->first;
          j++;
        }
      }
      interfaceLocalDof[ilevel].shrinkToFit(UINT_MAX);
      //END interface node search

      //BEGIN interface node dof global search, one for each soltype
      MyVector <unsigned> rowSize = interfaceLocalDof[ilevel].getRowSize();
      for (unsigned soltype = 0; soltype < 3; soltype++) {
        interfaceDof[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
        levelInterfaceSolidMark[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
        for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
          unsigned iel = interfaceElement[ilevel][i];
          unsigned counter = 0;
          for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
            unsigned jloc = interfaceLocalDof[ilevel][i][j];
            if (jloc < GetElementDofNumber(iel, soltype)) {
              unsigned jdof  = msh->GetSolutionDof(jloc, iel, soltype);
              interfaceDof[soltype][ilevel][i][counter] = jdof;
              unsigned jdof2  = msh->GetSolutionDof(jloc, iel, 2);
              levelInterfaceSolidMark[soltype][ilevel][i][counter] = msh->GetSolidMark(jdof2);
              counter++;
            }
            else {
              break;
            }
          }
        }
        interfaceDof[soltype][ilevel].shrinkToFit(UINT_MAX);
        levelInterfaceSolidMark[soltype][ilevel].shrinkToFit(UINT_MAX);
      }
      //END interface node dof global search, one for each soltype

      //BEGIN interface node coordinates
      interfaceNodeCoordinates[ilevel].resize(dim);
      for (unsigned k = 0; k < dim; k++) {
        interfaceNodeCoordinates[ilevel][k] = MyMatrix< double > (rowSize, 0.);
      }
      for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
        unsigned iel = interfaceElement[ilevel][i];
        for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
          unsigned jnode = interfaceLocalDof[ilevel][i][j];
          unsigned xDof  = msh->GetSolutionDof(jnode, iel, 2);
          for (unsigned k = 0; k < dim; k++) {
            interfaceNodeCoordinates[ilevel][k][i][j] = (*msh->_topology->_Sol[k])(xDof);
          }
        }
      }
      //END interface node coordinates search
    }

    for (unsigned soltype = 0; soltype < 3; soltype++) {
      for (int ilevel = 0; ilevel < _level; ilevel++) {
        for (int jlevel = ilevel + 1; jlevel <= _level; jlevel++) {
          for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
            interfaceDof[soltype][jlevel].broadcast(lproc);
            levelInterfaceSolidMark[soltype][jlevel].broadcast(lproc);
            for (unsigned d = 0; d < dim; d++) {
              interfaceNodeCoordinates[jlevel][d].broadcast(lproc);
            }
            std::map< unsigned, bool> candidateNodes;
            std::map< unsigned, bool> elementNodes;

            for (unsigned i = interfaceDof[soltype][ilevel].begin(); i < interfaceDof[soltype][ilevel].end(); i++) {

              candidateNodes.clear();

              std::vector < std::vector < std::vector <double > > > aP(3);
              bool aPIsInitialized = false;

              unsigned iel = interfaceElement[ilevel][i];
              short unsigned ielType = _elementType[iel];

              elementNodes.clear();
              for (unsigned j = 0; j < GetElementDofNumber(iel, soltype); j++) {
                unsigned jdof  = msh->GetSolutionDof(j, iel, soltype);
                elementNodes[jdof] = true;
              }

              std::vector < std::vector <double > > xv;
              msh->GetElementNodeCoordinates(xv, iel);
              unsigned ndofs = xv[0].size();

              double r;
              std::vector <double> xc;
              GetConvexHullSphere(xv, xc, r, 0.01);
              double r2 = r * r;

              std::vector < std::vector< double > > xe;
              GetBoundingBox(xv, xe, 0.01);


              for (unsigned k = interfaceDof[soltype][jlevel].begin(); k < interfaceDof[soltype][jlevel].end(); k++) {
                for (unsigned l = interfaceDof[soltype][jlevel].begin(k); l < interfaceDof[soltype][jlevel].end(k); l++) {
                  unsigned ldof = interfaceDof[soltype][jlevel][k][l];
                  if (candidateNodes.find(ldof) == candidateNodes.end() || candidateNodes[ldof] != false) {
                    double d2 = 0.;
                    std::vector<double> xl(dim);
                    for (int d = 0; d < dim; d++) {
                      xl[d] = interfaceNodeCoordinates[jlevel][d][k][l];
                      d2 += (xl[d] - xc[d]) * (xl[d] - xc[d]);
                    }
                    bool insideHull = true;
                    if (d2 > r2) {
                      insideHull = false;
                    }
                    for (unsigned d = 0; d < dim; d++) {
                      if (xl[d] < xe[d][0] || xl[d] > xe[d][1]) {
                        insideHull = false;
                      }
                    }
                    if (insideHull) {
                      if (elementNodes.find(ldof) == elementNodes.end()) {

                        if (!aPIsInitialized) {
                          aPIsInitialized = true;
                          std::vector < std::vector <double> > x1(dim);
                          for (unsigned jtype = 0; jtype < 3; jtype++) {
                            ProjectNodalToPolynomialCoefficients(aP[jtype], xv, ielType, jtype) ;
                          }
                        }

                        std::vector <double> xi;
                        GetClosestPointInReferenceElement(xv, xl, ielType, xi);
                        GetInverseMapping(2, ielType, aP, xl, xi);

                        bool insideDomain = CheckIfPointIsInsideReferenceDomain(xi, ielType, 0.0001);
                        if (insideDomain) {
                          for (unsigned j = interfaceDof[soltype][ilevel].begin(i); j < interfaceDof[soltype][ilevel].end(i); j++) {
                            unsigned jloc = interfaceLocalDof[ilevel][i][j];

                            basis* base = msh->GetBasis(ielType, soltype);
                            double value = base->eval_phi(jloc, xi);

                            if (fabs(value) >= 1.0e-10) {
                              unsigned jdof = interfaceDof[soltype][ilevel][i][j];
                              if (restriction[soltype][jdof].find(jdof) == restriction[soltype][jdof].end()) {
                                restriction[soltype][jdof][jdof] = 1.;
                                unsigned jdof2  = msh->GetSolutionDof(jloc, iel, 2);
                                interfaceSolidMark[soltype][jdof] = levelInterfaceSolidMark[soltype][ilevel][i][j];
                              }
                              restriction[soltype][jdof][ldof] = value;
                              restriction[soltype][ldof][ldof] = 10.;
                              interfaceSolidMark[soltype][ldof] = levelInterfaceSolidMark[soltype][jlevel][k][l];
                              candidateNodes[ldof] = true;
                            }
                          }
                        }
                        else {
                          candidateNodes[ldof] = false;
                        }
                      }
                      else {
                        candidateNodes[ldof] = false;
                      }
                    }
                  }
                }
              }
            }
            interfaceDof[soltype][jlevel].clearBroadcast();
            levelInterfaceSolidMark[soltype][jlevel].clearBroadcast();
            for (unsigned d = 0; d < dim; d++) {
              interfaceNodeCoordinates[jlevel][d].clearBroadcast();
            }
          }
        }
      }




      NumericVector* pvector;
      pvector = NumericVector::build().release();
      pvector->init(_nprocs, 1 , false, AUTOMATIC);

      unsigned counter = 1;
      while (counter != 0) {
        counter = 0;

        //BEGIN  saving the restriction object in parallel vectors and matrices

        MyVector <unsigned> rowSize(restriction[soltype].size(), 0);
        unsigned cnt1 = 0;
        for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
          rowSize[cnt1] = restriction[soltype][it1->first].size();
          cnt1++;
        }
        rowSize.stack();

        std::vector< unsigned > offset = rowSize.getOffset();

        MyVector <unsigned> masterNode(offset);
        MyMatrix <unsigned> slaveNodes(rowSize);
        MyMatrix <double> slaveNodesValues(rowSize);

        cnt1 = 0;
        for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
          masterNode[offset[_iproc] + cnt1] = it1->first;
          unsigned cnt2 = 0;
          for (std::map<unsigned, double> ::iterator it2 = restriction[soltype][it1->first].begin(); it2 != restriction[soltype][it1->first].end(); it2++) {
            slaveNodes[offset[_iproc] + cnt1][cnt2] = it2->first;
            slaveNodesValues[offset[_iproc] + cnt1][cnt2] = it2->second;
            cnt2++;
          }
          cnt1++;
        }
        //END  saving the restriction object in parallel vectors and matrices


        //BEGIN filling the restriction object with infos coming form the parallel vectors and matrices
        unsigned solutionOffset = msh->_dofOffset[soltype][_iproc];
        unsigned solutionOffsetp1 = msh->_dofOffset[soltype][_iproc + 1];
        for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
          masterNode.broadcast(lproc);
          slaveNodes.broadcast(lproc);
          slaveNodesValues.broadcast(lproc);
          for (unsigned i = slaveNodes.begin(); i < slaveNodes.end(); i++) {
            unsigned inode = masterNode[i];
            if (inode >= solutionOffset && inode < solutionOffsetp1 && // inode belongs to _iproc
                restriction[soltype].find(inode) == restriction[soltype].end()) { // but inode is not set as master node of _iproc
              counter++;
              for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { //copy information for lproc to _iproc
                unsigned jnode = slaveNodes[i][j];
                restriction[soltype][inode][jnode] = slaveNodesValues[i][j];
              }
            }
            else { // either inode does not belong to _iproc or it was already defined as a master for _iproc
              if (restriction[soltype].find(inode) != restriction[soltype].end()) { // inode is already defined as master node in _iproc (either it does or does not belong to _iproc)
                for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { // loop on all the columns of restriction[lproc][inode]
                  unsigned jnode = slaveNodes[i][j];
                  double value = slaveNodesValues[i][j];
                  if (inode != jnode || value > 5.) { // if off-diagonal or hanging node for lproc
                    restriction[soltype][inode][jnode] =  value;
                  }
                  if (restriction[soltype].find(jnode) == restriction[soltype].end()) { // if jnode is not yet a master node for _iproc
                    counter++;
                    for (unsigned k = masterNode.begin(); k < masterNode.end(); k++) {
                      if (masterNode[k] == jnode) { // and if jnode is also a master node for lproc
                        for (unsigned l = slaveNodes.begin(k); l < slaveNodes.end(k); l++) {
                          unsigned lnode = slaveNodes[k][l];
                          restriction[soltype][jnode][lnode] = slaveNodesValues[k][l]; //copy the rule of lptoc into _iproc
                        }
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
          masterNode.clearBroadcast();
          slaveNodes.clearBroadcast();
          slaveNodesValues.clearBroadcast();
        }
        //END filling the restriction object with infos coming form the parallel vectors and matrices


        pvector->set(_iproc, counter);
        pvector->close();
        counter = static_cast <unsigned>(floor(pvector->l1_norm() + 0.5));
      }
      delete pvector;

//       for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
//         unsigned inode = it1->first;
//         if (restriction[soltype][inode][inode] > 5.) {
//           if (restriction[soltype][inode].size() > 1) {
//             for (std::map<unsigned, std::map<unsigned, double> >::iterator it2 = restriction[soltype].begin(); it2 != restriction[soltype].end(); it2++) {
//               unsigned jnode = it2->first;
//               if (jnode != inode && restriction[soltype][jnode].find(inode) != restriction[soltype][jnode].end()) {
//                 double value =  restriction[soltype][jnode][inode];
//                 for (std::map<unsigned, double> ::iterator it3 = restriction[soltype][inode].begin(); it3 != restriction[soltype][inode].end(); it3++) {
//                   unsigned knode = it3->first;
//                   if (knode != inode) {
//                     restriction[soltype][jnode][knode] = it3->second * value;
//                   }
//                 }
//               }
//             }
//           }
//           restriction[soltype][inode].clear();
//           restriction[soltype][inode][inode] = 0.;
//         }
//       }



      std::vector<std::vector < unsigned > > genealogy;
      std::vector<std::vector < double > > heredity;
      std::vector< unsigned > index;
      std::map < unsigned,  std::map < unsigned, double  > >  restrictionCopy = restriction[soltype];

      for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restrictionCopy.begin(); it1 != restrictionCopy.end(); it1++) { // loop all over master, hanging and master+hanging nodes

        genealogy.resize(1);
        heredity.resize(1);
        index.resize(1);

        genealogy[0].resize(1);
        heredity[0].resize(1);

        unsigned inode = it1->first;

        if (restrictionCopy[inode][inode] < 5.) { // only if a real master node

          // initialize master node genealogy and heredity at level 0
          genealogy[0][0] = inode;
          heredity[0][0] = 1.;
          index[0] = 0;

          restriction[soltype][inode].clear();
          restriction[soltype][inode][inode] = 1.;

          unsigned level = 1;
          while ( level > 0 ) {

            // initialize master node genealogy and heredity at genemeric level

            unsigned father = genealogy[level - 1][index[level - 1]];

            genealogy.resize(level + 1);
            genealogy[level].reserve( restrictionCopy[ father ].size() - 1 );
            genealogy[level].resize(0);

            heredity.resize(level + 1);
            heredity[level].reserve( restrictionCopy[ father ].size() - 1 );
            heredity[level].resize(0);

            index.resize(level + 1);
            index[level] = 0;

            unsigned cnt  = 0;
            for (std::map <unsigned, double>::iterator it2 = restrictionCopy[ father ].begin(); it2 != restrictionCopy[ father ].end(); it2++) { // loop on all the father sons
              unsigned son = it2->first;
              bool alreadyFound = false;
              for (unsigned klevel = 0; klevel < level; klevel++) { // check if the son is in the previous genealogy
                for (unsigned k = 0; k < genealogy[klevel].size(); k++) {
                  if ( genealogy[klevel][k] == son ) alreadyFound = true;
                }
              }
              if (!alreadyFound ) { // if never found add the the restionction value in the master node line and zero the hanging node line
                genealogy[level].resize( genealogy[level].size() + 1);
                heredity[level].resize(heredity[level].size() + 1);

                genealogy[level][cnt] = son;
                heredity[level][cnt] = it2->second * heredity[level - 1][index[level - 1]];

                restriction[soltype][inode][son] += heredity[level][cnt];
                cnt++;

                restriction[soltype][son].clear();
                restriction[soltype][son][son] = 0.;
              }
            }

            if ( cnt > 0) {
              level++;
            }
            else {
              bool test = true;
              while ( test && level > 0) {
                index[level - 1]++;
                test = false;
                if ( index[level - 1] == genealogy[level - 1].size() ) {
                  level--;
                  test = true;
                }
              }
            }
          }
        }
	else{
	  restriction[soltype][inode].clear();
	  restriction[soltype][inode][inode] = 0.;
	}
      }

      MyVector <unsigned> InterfaceSolidMarkNode(interfaceSolidMark[soltype].size());
      MyVector <short unsigned> InterfaceSolidMarkValue(interfaceSolidMark[soltype].size());

      unsigned cnt = 0;
      for (std::map<unsigned, bool >::iterator it = interfaceSolidMark[soltype].begin(); it != interfaceSolidMark[soltype].end(); it++) {
        InterfaceSolidMarkNode[cnt] = it->first;
        InterfaceSolidMarkValue[cnt] = it->second;
        cnt++;
      }
      InterfaceSolidMarkNode.stack();
      InterfaceSolidMarkValue.stack();

      for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
        InterfaceSolidMarkNode.broadcast(lproc);
        InterfaceSolidMarkValue.broadcast(lproc);
        for (unsigned i = InterfaceSolidMarkNode.begin(); i < InterfaceSolidMarkNode.end(); i++) {
          unsigned jnode = InterfaceSolidMarkNode[i];
          if ( restriction[soltype].find(jnode) != restriction[soltype].end()) {
            interfaceSolidMark[soltype][jnode] = InterfaceSolidMarkValue[i];
          }
        }
        InterfaceSolidMarkNode.clearBroadcast();
        InterfaceSolidMarkValue.clearBroadcast();
      }
    }
  }


  void Mesh::GetElementNodeCoordinates(std::vector < std::vector <double > > &xv, const unsigned &iel, const unsigned &solType) {
    xv.resize(_dimension);
    unsigned ndofs = el->GetElementDofNumber(iel, solType);
    for(int d = 0; d < _dimension; d++) {
      xv[d].resize(ndofs);
    }
    for(unsigned j = 0; j < ndofs; j++) {
      unsigned xdof  = GetSolutionDof(j, iel, solType);
      for(int d = 0; d < _dimension; d++) {
        xv[d][j] = (*_topology->_Sol[d])(xdof);
      }
    }
  }



} //end namespace femus




