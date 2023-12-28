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

#include "Elem.hpp"
#include "Mesh.hpp"
#include "GeomElTypeEnum.hpp"
#include "NumericVector.hpp"

#include "Basis.hpp"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cassert>
#include <climits>



namespace femus
{


  /**
   * This constructor allocates the memory for the \textit{coarsest elem}
   **/
  elem::elem(const unsigned& other_nel, const unsigned dim_in)
  {
      
    _dim = dim_in;
      
    _coarseElem = NULL;

    _level = 0;

    InitializeNumberOfElementsPerGeomType();
    
    _nel = other_nel;
    

    ResizeElement_Level_Type_Group_Material(_nel, _level);
    
    
    _elementDof.resize(_nel, NVE[0][2], UINT_MAX);
    
    _elementNearFace.resize(_nel, NFC[0][1], -1);

}

  
  /**
   * This constructor allocates the memory for the \textit{finer elem}
   * starting from the parameters of the \textit{coarser elem}
   **/
  elem::elem(elem* elc, const unsigned dim_in, const unsigned refindex, const std::vector < double >& coarseAmrVector)
  {
      
    _dim = dim_in;
      
    _coarseElem = elc;

    _level = elc->_level + 1;

    
    InitializeNumberOfElementsPerGeomType();
    
    _nel = InitializeNumberOfElementsFromCoarseList(elc, refindex);
    

    ResizeElement_Level_Type_Group_Material(_nel, _level);
    

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


  

   unsigned int  elem::InitializeNumberOfElementsFromCoarseList(elem* elc, const unsigned refindex)  {
     
     unsigned int nelem;
     
    nelem = elc->GetRefinedElementNumber() * refindex; //refined
    nelem += elc->GetElementNumber() - elc->GetRefinedElementNumber(); // + non-refined;
    
    return nelem;
    
   }

  
  
  void elem::ShrinkToFitElementDof() {
    
    _elementDof.shrinkToFit(UINT_MAX);
    
  }
  
  void elem::ShrinkToFitElementNearFace() {
    

    MyVector <unsigned> rowSize(_nel);
    for (unsigned iel = 0; iel < _nel; iel++) {
      unsigned ielType = GetElementType(iel);
      rowSize[iel] = NFC[ielType][1];
    }
    _elementNearFace.shrinkToFit(rowSize);
    
  }
  
  
  void elem::ReorderElementNearFace_rows(const std::vector < unsigned >& elementMapping) {
    
    //BEGIN reordering _elementNearFace (rows)
    MyVector <unsigned> rowSize(_nel, 0);
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
    //END reordering _elementNearFace (rows)
  }
  
  
  void elem::ReorderElementDof_rows(const std::vector < unsigned >& elementMapping) {

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
    //END reordering OF _elementDof (rows)
    
  }

  
  void elem::ReorderChildElement_OnCoarseElem_columns(const std::vector < unsigned >& elementMapping) {
    
     //BEGIN reordering _childElementDof (columns) on coarse levels
    if (_level != 0) {
      for (unsigned i = _coarseElem->_childElem.begin(); i < _coarseElem->_childElem.end(); i++) {
        for (unsigned j = _coarseElem->_childElem.begin(i); j < _coarseElem->_childElem.end(i); j++) {
          _coarseElem->_childElem[i][j] =  elementMapping[ _coarseElem->_childElem[i][j]];
        }
      }
    }
    //END reordering _childElementDof (columns) on coarse levels
 
  } 
  
  void elem::ReorderMeshElement_Type_Level_Group_Material___NearFace_rows_ChildElem_columns(const std::vector < unsigned >& elementMapping)
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
    //END reordering  _elementLevel

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

    
    //Elem
    //BEGIN reordering _elementNearFace (rows)
    ReorderElementNearFace_rows(elementMapping);
    //END reordering _elementNearFace (rows)

    //Elem
    //BEGIN reordering _childElement (columns) on coarse levels
    ReorderChildElement_OnCoarseElem_columns(elementMapping);
    //END reordering _childElement (columns) on coarse levels
 
    
  }


  void elem::ReorderMeshElement_Dof_stuff(const std::vector < unsigned >& elementMapping)
  {
    
    //BEGIN reordering _elementDof (rows)
    ReorderElementDof_rows(elementMapping);
    //END reordering OF _elementDof (rows)

    
  }
  
  
  
  
  
  void elem::ReorderElementDof_columns_Using_node_mapping(const std::vector < unsigned >& nodeMapping)
  {
    for (unsigned i = _elementDof.begin(); i < _elementDof.end(); i++) {
      for (unsigned j = _elementDof.begin(i); j < _elementDof.end(i); j++) {
        _elementDof[i][j] =  nodeMapping[_elementDof[i][j]];
      }
    }
  }

  void elem::DeleteElementNearVertex()
  {
    _elementNearVertex.clear();
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
  unsigned elem::GetFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& inode) const
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
   * Return the total number of elements, where the shape can be specified as an input string
   **/
  unsigned elem::GetElementNumber(const std::string name) const
  {
    if (!strcmp(name.c_str(), "All")) {
      return _nel;
    }
    unsigned i;
    i = this->GetIndex(name);
    return _nelt[i];
  }

  /**
   * Add value to the total number of the element
   **/
  void elem::AddToElementNumber(const unsigned& value, const std::string name)
  {
    unsigned i;
    i = GetIndex(name);
    _nelt[i] += value;
  }
  
  void elem::AddToElementNumber(const unsigned& value, short unsigned ielt)
  {
    _nelt[ielt] += value;
  }
  


  /**
   * Return the global adiacent-to-face element number
   **/
  int elem::GetFaceElementIndex(const unsigned& iel, const unsigned& iface) const
  {
    return _elementNearFace[iel][iface];
  }

  int elem::GetBoundaryIndex(const unsigned& iel, const unsigned& iface) const
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
  const short unsigned elem::GetElementType(const unsigned& iel) const
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
  short unsigned elem::GetElementGroup(const unsigned& iel) const
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
  short unsigned elem::GetElementMaterial(const unsigned& iel) const
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

  
  void elem::BuildElementNearFace() {
      
    for(unsigned iel = 0; iel < GetElementNumber(); iel++) {
      for(unsigned iface = 0; iface < GetElementFaceNumber(iel); iface++) {
        if(GetFaceElementIndex(iel, iface) <= 0) {   /// @todo probably just == -1
          unsigned i1 = GetFaceVertexIndex(iel, iface, 0);
          unsigned i2 = GetFaceVertexIndex(iel, iface, 1);
          unsigned i3 = GetFaceVertexIndex(iel, iface, 2);

          for(unsigned j = 0; j < GetElementNearVertexNumber(i1); j++) {
            unsigned jel = GetElementNearVertex(i1, j);

            if(jel > iel) {
              for(unsigned jface = 0; jface < GetElementFaceNumber(jel); jface++) {
                if(GetFaceElementIndex(jel, jface) <= 0) {
                  unsigned j1 = GetFaceVertexIndex(jel, jface, 0);
                  unsigned j2 = GetFaceVertexIndex(jel, jface, 1);
                  unsigned j3 = GetFaceVertexIndex(jel, jface, 2);
                  unsigned j4 = GetFaceVertexIndex(jel, jface, 3);
                  
                  const bool faces_coincide_three_dim = ( GetDimension() == 3 &&
                                         (i1 == j1 || i1 == j2 || i1 == j3 ||  i1 == j4) &&
                                         (i2 == j1 || i2 == j2 || i2 == j3 ||  i2 == j4) &&
                                         (i3 == j1 || i3 == j2 || i3 == j3 ||  i3 == j4));
                  
                  const bool faces_coincide_two_dim = ( GetDimension() == 2 &&
                                       (i1 == j1 || i1 == j2) &&
                                       (i2 == j1 || i2 == j2));
                  
                  const bool faces_coincide_one_dim = ( GetDimension() == 1 && (i1 == j1));

                  if(faces_coincide_three_dim
                      ||
                     faces_coincide_two_dim 
                      ||
                     faces_coincide_one_dim
                    ) {
                    SetFaceElementIndex(iel, iface, jel + 1u);
                    SetFaceElementIndex(jel, jface, iel + 1u);
                  }
                }
              }
            }
          }
        }
      }
    }
    
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

  
  
  void elem::BuildElem_NearFace_NearElem_using_NearVertex() {
    
      
    BuildElementNearVertex();

       BuildElementNearFace();     //needs ElementNearVertex

       BuildElementNearElement();  //needs ElementNearVertex
    
    DeleteElementNearVertex();

  }
  
  
  void elem::ScatterElement_Level_Type_Group_Material___NearFace() {

    
    ScatterElement_Level_Type_Group_Material();
    
    ScatterElementNearFace();

    
  }
  
  
  
  /**
   * Return the number of elements which have a given vertex
   **/
  unsigned elem::GetElementNearVertexNumber(const unsigned& inode) const
  {
    return _elementNearVertex.size(inode);
  }

  /**
   * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
   **/
  unsigned elem::GetElementNearVertex(const unsigned& inode, const unsigned& j) const
  {
    return _elementNearVertex[inode][j];
  }

  /**
   * return the index 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  unsigned elem::GetIndex(const std::string name) const
  {
    unsigned index = 0;
    if (!strcmp(name.c_str(), geom_elems[HEX].c_str() )) {
      index = HEX;
    }
    else if (!strcmp(name.c_str(), geom_elems[TET].c_str() )) {
      index = TET;
    }
    else if (!strcmp(name.c_str(), geom_elems[WEDGE].c_str() )) {
      index = WEDGE;
    }
    else if (!strcmp(name.c_str(), geom_elems[QUAD].c_str() )) {
      index = QUAD;
    }
    else if (!strcmp(name.c_str(), geom_elems[TRI].c_str() )) {
      index = TRI;
    }
    else if (!strcmp(name.c_str(), geom_elems[LINE].c_str() ) ) {
      index = LINE;
    }
    else {
      std::cout << "error! invalid Element Shape in elem::GetIndex(...)" << std::endl;
      exit(0);
    }
    return index;
  }
  
  
  void elem::AllocateChildrenElement(const unsigned& refindex, const Mesh* msh) {
    
    MyVector <unsigned> rowSize(_elementOffset, 0);
    
    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[i] = (msh->GetRefinedElementIndex(i) == 1) ? refindex : 1; //for every element, establish if it is refined or not, and then return the number of children
    }
    _childElem = MyMatrix <unsigned> (rowSize, 0);
    
  }

  void elem::AllocateChildrenElementDof(const unsigned& refindex, const Mesh* msh) {

    MyVector <unsigned> rowSize(_elementOffset, 0);
    
    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      const unsigned elementType = GetElementType(i);
      rowSize[i] = (msh->GetRefinedElementIndex(i) == 1) ? refindex * NVE[elementType][2] : NVE[elementType][2];  //for every element, establish if it is refined or not, and return the number of DofCarriers of all its children
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

  unsigned elem::GetChildElementDof(const unsigned& iel, const unsigned& i0, const unsigned i1) const
  {
    return _childElemDof[iel][i0 * GetElementDofNumber(iel, 2) + i1];
  }

  void elem::SetChildElement(const unsigned& iel, const unsigned& json, const unsigned& value)
  {
    _childElem[iel][json] = value;
  }

  unsigned elem::GetChildElement(const unsigned& iel, const unsigned& json) const
  {
    return _childElem[iel][json];
  }







  void elem::ScatterElementDof()
  {
    _elementDof.scatter(_elementOffset);
  }

  void elem::LocalizeElementDof(const unsigned& jproc)
  {
    _elementDof.broadcast(jproc);
  }

  unsigned elem::GetElementDofIndex(const unsigned& iel, const unsigned& inode) const
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
  



} //end namespace femus




