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
  elem::elem( const unsigned& other_nel ) {

    _coarseElem = NULL;

    _level = 0;

    _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
    _nel = other_nel;

    _elementType.resize(_nel);
    _elementGroup.resize(_nel);
    _elementMaterial.resize(_nel);
    _elementLevel.resize(_nel,_level);
   
    _elementDof.resize(_nel,NVE[0][2]);
    _elementNearFace.resize(_nel,NFC[0][1]);
     
    _childElemFlag = false;

    _nelr = _nelrt[0] = _nelrt[1] = _nelrt[2] = _nelrt[3] = _nelrt[4] = _nelrt[5] = 0;
 
  }


  void elem::SharpMemoryAllocation() {

    MyVector <unsigned> rowSizeElDof(_nel);
    MyVector <unsigned> rowSizeElNearFace(_nel);
    for( unsigned iel = 0; iel < _nel; iel++ ) {
      unsigned ielType = GetElementType( iel );
      rowSizeElDof[iel] = NVE[ielType][2];
      rowSizeElNearFace[iel] = NFC[ielType][1];
    }
    
    MyMatrix <unsigned> tmpElDof(rowSizeElDof);
    MyMatrix <int> tmpElNearFace(rowSizeElNearFace);
    
    rowSizeElDof.clear();
    rowSizeElNearFace.clear();
        
    for(unsigned i=tmpElDof.begin();i<tmpElDof.end();i++){
      for(unsigned j=tmpElDof.begin(i);j<tmpElDof.end(i);j++){
	tmpElDof[i][j] = _elementDof[i][j];
      }
    }
    
    for(unsigned i=tmpElNearFace.begin();i<tmpElNearFace.end();i++){
      for(unsigned j=tmpElNearFace.begin(i);j<tmpElNearFace.end(i);j++){
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
  //elem::elem( elem* elc, const unsigned refindex, const std::vector < double >& coarseAmrVector, const std::vector < double >& coarseElementType ) {
  elem::elem( elem* elc, const unsigned refindex, const std::vector < double >& coarseAmrVector) {
    _coarseElem = elc;

    _level = elc->_level + 1;

    _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
    _nel = elc->GetRefinedElementNumber() * refindex; //refined
    _nel += elc->GetElementNumber() - elc->GetRefinedElementNumber(); // + non-refined;

    _elementType.resize(_nel);
    _elementGroup.resize(_nel);
    _elementMaterial.resize(_nel);
    _elementLevel.resize(_nel,_level);
        
    //**************************
    MyVector <unsigned> rowSizeElDof(_nel);
    MyVector <unsigned> rowSizeElNearFace(_nel);
    unsigned jel = 0;
    for(unsigned isdom = 0; isdom < elc->_nprocs; isdom++) {
      elc->_elementType.localize(isdom);
      for(unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
   	short unsigned elemt = elc->_elementType[iel];
	int increment = 1;
	if( static_cast < short unsigned >( coarseAmrVector[iel] + 0.25 ) == 1 ) {
	  increment = NRE[elemt];
	}
	for( unsigned j = 0; j < increment; j++ ) {
	  rowSizeElDof[jel + j] += NVE[elemt][2];
	  rowSizeElNearFace[jel + j] += NFC[ elemt ][1];
	}
	jel += increment;
      }
      elc->_elementType.clearLocalized();
    }
    _elementDof = MyMatrix <unsigned> (rowSizeElDof);
    _elementNearFace = MyMatrix <int> (rowSizeElNearFace,-1);
    
    rowSizeElDof.clear();
    rowSizeElNearFace.clear();
        
    _childElemFlag = false;

    _nelr = _nelrt[0] = _nelrt[1] = _nelrt[2] = _nelrt[3] = _nelrt[4] = _nelrt[5] = 0;
  
  }

  void elem::ReorderMeshElements( const std::vector < unsigned >& elementMapping ) {

    //BEGIN reordering _elementType
    MyVector <short unsigned> tempElementType;
    tempElementType = _elementType;
    //_elementType = new short unsigned [_nel];
    for( unsigned iel = 0; iel < _nel; iel++ ) {
      _elementType[ elementMapping [iel] ]  = tempElementType[iel] ;
    }
    tempElementType.clear();
    //END reordering _elementType

    //BEGIN reordering _elementLevel
    if( _level != 0 ) {      
      MyVector <short unsigned> tempElementLevel;
      tempElementLevel = _elementLevel;
      for( unsigned iel = 0; iel < _nel; iel++ ) {
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
      
      for( unsigned iel = 0; iel < _nel; iel++ ) {
        _elementGroup[elementMapping[iel]]   = tempElementGroup[iel];
        _elementMaterial[elementMapping[iel]] = tempElementMaterial[iel] ;
      }
      tempElementGroup.clear();
      tempElementMaterial.clear();
 
    //END reordering _elementGroup and _elementMaterial

      
    //BEGIN reordering _elementDof (rows)
    MyMatrix <unsigned> tmpElDof = _elementDof;
    for(unsigned i=_elementDof.begin();i<_elementDof.end();i++){
      for(unsigned j=_elementDof.begin(i);j<_elementDof.end(i);j++){
	_elementDof[elementMapping [i]][j] = tmpElDof[i][j];
      }
    }
    tmpElDof.clear();
    //END reordering OF _elementDof
    
    //BEGIN reordering _elementNearFace (rows)
    MyMatrix <int> tmpElNearFace = _elementNearFace;
    for(unsigned i=_elementNearFace.begin();i<_elementNearFace.end();i++){
      for(unsigned j=_elementNearFace.begin(i);j<_elementNearFace.end(i);j++){
	_elementNearFace[elementMapping [i]][j] = tmpElNearFace[i][j];
      }
    }
    tmpElNearFace.clear();
    //END reordering _elementNearFace
    
    //BEGIN reordering _childElementDof (columns) on coarse level
    if( _level != 0 ) {
      for( unsigned i = 0; i < _coarseElem->_childElemMemorySize; i++ ) {
        _coarseElem->_childElemMemory[i] =  elementMapping[ _coarseElem->_childElemMemory[i]];
      }
    }
    //END reordering _childElementDof
  }


  void elem::ReorderMeshNodes( const std::vector < unsigned >& nodeMapping ) {  
    for( unsigned i = _elementDof.begin(); i < _elementDof.end(); i++ ) {
      for( unsigned j = _elementDof.begin(i); j < _elementDof.end(i); j++ ) {
	_elementDof[i][j] =  nodeMapping[ _elementDof[i][j] ];
      }
    }   
  }

  elem::~elem() {
    if( _childElemFlag ) {
      delete [] _childElemMemory;
      delete [] _childElem;

      delete [] _childElemDofMemory;
      delete [] _childElemDofMemoryPointer;
      delete [] _childElemDof;
    }
  }

  void elem::DeleteElementNearVertex() {
    delete [] _elementNearVertexMemory;
    delete [] _elementNearVertex;
    delete [] _elementNearVertexNumber;
  }

  /**
   * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
   **/
  unsigned elem::GetElementDofNumber( const unsigned& iel, const unsigned& type ) {
    return NVE[_elementType[iel]][type];
  }


  /**
   * Set the local->global node number
   **/
  void elem::SetElementDofIndex( const unsigned& iel, const unsigned& inode, const unsigned& value ) {
    //_elementDof[iel][inode] = value;
    _elementDof[iel][inode] = value;
  }

  /**
   * Return the local->global face node index
   **/
  unsigned elem::GetFaceVertexIndex( const unsigned& iel, const unsigned& iface, const unsigned& inode ) {
    return _elementDof[iel][ig[_elementType[iel]][iface][inode]];
    //return _elementDof[iel][ig[_elementType[iel]][iface][inode]];
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
  void elem::SetNodeNumber( const unsigned& value ) {
    _nvt = value;
  }

  /**
   * Return the total number of elements
   **/
  unsigned elem::GetElementNumber( const char* name ) const {
    if( !strcmp( name, "All" ) ) {
      return _nel;
    }
    unsigned i;
    i = this->GetIndex( name );
    return _nelt[i];
  }

  /**
   * Add value to the total number of the element
   **/
  void elem::AddToElementNumber( const unsigned& value, const char name[] ) {
    unsigned i;
    i = GetIndex( name );
    _nelt[i] += value;
  }
  void elem::AddToElementNumber( const unsigned& value, short unsigned ielt ) {
    _nelt[ielt] += value;
  }
  unsigned elem::GetElementFaceNumber( const unsigned& iel, const unsigned& type ) {
    return NFC[ _elementType[iel] ][type];
  }

  /**
   * Return the global adiacent-to-face element number
   **/
  int elem::GetFaceElementIndex( const unsigned& iel, const unsigned& iface ) {
    return _elementNearFace[iel][iface];
  }

  int elem::GetBoundaryIndex( const unsigned& iel, const unsigned& iface ) {
    return  -( GetFaceElementIndex( iel, iface ) + 1 );
  }

  /**
   * Set the global adiacent-to-face element number
   **/
  void elem::SetFaceElementIndex( const unsigned& iel, const unsigned& iface, const int& value ) {
    //_elementNearFace[iel][iface] = value;
    _elementNearFace[iel][iface] = value;
  }

  /**
   * Return element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  short unsigned elem::GetElementType( const unsigned& iel ) {
    return _elementType[iel];
  }

  /**
   * Set element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  void elem::SetElementType( const unsigned& iel, const short unsigned& value ) {
    _elementType[iel] = value;
  }

  /**
   * Return element group
   **/
  short unsigned elem::GetElementGroup( const unsigned& iel ) {
    return _elementGroup[iel];
  }

  /**
   * Set element group
   **/
  void elem::SetElementGroup( const unsigned& iel, const short unsigned& value ) {
    _elementGroup[iel] = value;
  }

  /**
   * Set element Material
  **/
  void elem::SetElementMaterial( const unsigned& iel, const short unsigned& value ) {
    _elementMaterial[iel] = value;
  }

  /**
   * Return element material
   **/
  short unsigned elem::GetElementMaterial( const unsigned& iel ) {
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
  void elem::SetElementGroupNumber( const unsigned& value ) {
    _ngroup = value;
  }

  /**
   * Set the memory storage and initialize nve and kvtel (node->element vectors)
   **/

  void elem::BuildLocalElementNearVertex() {
    for( unsigned iel = _elementOffset[ _iproc ]; iel < _elementOffset[ _iproc + 1]; iel++ ) {
      for( unsigned i = 0; i < GetElementDofNumber( iel, 0 ); i++ ) {
        unsigned inode = GetElementDofIndex( iel, i );
        unsigned inodesize = 1;
        if( _localElementNearVertexMap.find( inode ) != _localElementNearVertexMap.end() ) {
          inodesize = _localElementNearVertexMap[inode].size();
          inodesize++;
        }
        else {
          _localElementNearVertexMap[inode].reserve( _elementNearVertexNumber[inode] );
        }
        _localElementNearVertexMap[inode].resize( inodesize );
        _localElementNearVertexMap[inode][inodesize - 1] = iel;
      }
    }
  }

  void elem::BuildElementNearVertex() {
    unsigned counter = ( _nelt[0] * NVE[0][0] + _nelt[1] * NVE[1][0] + _nelt[2] * NVE[2][0] +
                         _nelt[3] * NVE[3][0] + _nelt[4] * NVE[4][0] + _nelt[5] * NVE[5][0] );

    _elementNearVertex = new unsigned * [_nvt];
    _elementNearVertexNumber = new unsigned[_nvt];

    for( unsigned inode = 0; inode < _nvt; inode++ ) {
      _elementNearVertexNumber[inode] = 0;
    }
    for( unsigned iel = 0; iel < _nel; iel++ ) {
      for( unsigned inode = 0; inode < GetElementDofNumber( iel, 0 ); inode++ ) {
        _elementNearVertexNumber[ GetElementDofIndex( iel, inode )]++;
      }
    }

    _elementNearVertexMemory = new unsigned[counter];
    unsigned* pt = _elementNearVertexMemory;
    for( unsigned inode = 0; inode < _nvt; inode++ ) {
      _elementNearVertex[inode] = pt;
      pt += _elementNearVertexNumber[inode];
      for( unsigned j = 0; j < _elementNearVertexNumber[inode]; j++ ) {
        _elementNearVertex[inode][j] = _nel;
      }
    }

    for( unsigned iel = 0; iel < _nel; iel++ ) {
      for( unsigned inode = 0; inode < GetElementDofNumber( iel, 0 ); inode++ ) {
        unsigned irow = GetElementDofIndex( iel, inode );
        unsigned j = 0;
        while( _nel != _elementNearVertex[irow][j] ) j++;
        _elementNearVertex[irow][j] = iel;
      }
    }
  }

  /**
   * Return the number of elements which have given node
   **/
  unsigned elem::GetElementNearVertexNumber( const unsigned& inode )const {
    return _elementNearVertexNumber[inode];
  }

  /**
   * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
   **/
  unsigned elem::GetElementNearVertex( const unsigned& inode, const unsigned& j )const {
    return _elementNearVertex[inode][j];
  }

  /**
   * return the index 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  unsigned elem::GetIndex( const char name[] ) const {
    unsigned index = 0;
    if( !strcmp( name, "Hex" ) ) {
      index = 0;
    } else if( !strcmp( name, "Tet" ) ) {
      index = 1;
    } else if( !strcmp( name, "Wedge" ) ) {
      index = 2;
    } else if( !strcmp( name, "Quad" ) ) {
      index = 3;
    } else if( !strcmp( name, "Triangle" ) ) {
      index = 4;
    } else if( !strcmp( name, "Line" ) ) {
      index = 5;
    } else {
      cout << "error! invalid Element Shape in elem::GetIndex(...)" << endl;
      exit( 0 );
    }
    return index;
  }


  void elem::AllocateChildrenElement( const unsigned& refindex, Mesh* msh ) {
    if( _childElemFlag ) {
      delete [] _childElemMemory;
      delete [] _childElem;

      delete [] _childElemDofMemory;
      delete [] _childElemDofMemoryPointer;
      delete [] _childElemDof;
    }

    unsigned localNel = _elementOwned;
    unsigned localNelr = 0;


    _childElemDofMemorySize = 0;
    for( unsigned iel = _elementOffset[ _iproc ]; iel < _elementOffset[ _iproc + 1]; iel++ ) {
      unsigned elementType = msh->GetElementType( iel );
      if( msh->GetRefinedElementIndex( iel ) == 1 ) {
        localNelr++;
        _childElemDofMemorySize += refindex * NVE[elementType][2];
      }
      else {
        _childElemDofMemorySize += NVE[elementType][2];
      }
    }

    _childElemMemorySize = localNelr * refindex + ( localNel - localNelr );

    _childElemMemory = new unsigned [_childElemMemorySize];
    memset( _childElemMemory, 0, _childElemMemorySize * sizeof( unsigned ) );
    _childElem = new unsigned* [localNel];

    _childElemDofMemory = new unsigned [_childElemDofMemorySize];
    _childElemDofMemoryPointer = new unsigned *[_childElemMemorySize];
    _childElemDof = new unsigned** [localNel];


    unsigned* ptr = _childElemMemory;

    unsigned* ptr1  = _childElemDofMemory;
    unsigned** pptr = _childElemDofMemoryPointer;
    unsigned counter = 0;

    for( unsigned i = 0; i < localNel; i++ ) {
      _childElem[i] = ptr;

      unsigned elementType = msh->GetElementType( _elementOffset[ _iproc ] + i );

      _childElemDof[i] = pptr;
      if( msh->GetRefinedElementIndex( _elementOffset[ _iproc ] + i ) == 1 ) {
        ptr += refindex;

        pptr += refindex;
        for( unsigned j = 0; j < refindex; j++ ) {
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

  void elem::SetChildElementDof( const unsigned& refIndex, Mesh* msh, elem* elf ) {
    for( unsigned iel = _elementOffset[ _iproc ]; iel < _elementOffset[ _iproc + 1]; iel++ ) {
      unsigned elementType = msh->GetElementType( iel );
      unsigned endIndex = ( msh->GetRefinedElementIndex( iel ) == 1 ) ? refIndex : 1u;
      for( unsigned j = 0; j < endIndex; j++ ) {
        unsigned ielf = GetChildElement( iel, j );
        for( unsigned k = 0; k < elf->GetElementDofNumber( ielf, 2 ); k++ ) {
          _childElemDof[iel - _elementOffset[ _iproc ]][j][k] = elf->GetElementDofIndex( ielf, k );
        }
      }
    }
  }

  void elem::SetChildElement( const unsigned& iel, const unsigned& json, const unsigned& value ) {
    _childElem[ iel - _elementOffset[ _iproc ]][json] = value;
    return;
  }

  unsigned elem::GetChildElement( const unsigned& iel, const unsigned& json ) const {
    return _childElem[iel - _elementOffset[ _iproc ] ][json];
  }

  const unsigned elem::GetNVE( const unsigned& elementType, const unsigned& doftype ) const {
    return NVE[elementType][doftype];
  }

  const unsigned elem::GetNFACENODES( const unsigned& elementType, const unsigned& jface, const unsigned& dof ) const {
    return NFACENODES[elementType][jface][dof];
  }

  const unsigned elem::GetNFC( const unsigned& elementType, const unsigned& type ) const {
    return NFC[elementType][type];
  }

  const unsigned elem::GetIG( const unsigned& elementType, const unsigned& iface, const unsigned& jnode ) const {
    return ig[elementType][iface][jnode];
  }

  void elem::ScatterElementDof() {
    _elementDof.scatter(_elementOffset);
  }

  void elem::LocalizeElementDofFromOneToAll( const unsigned& jproc ) {
    _elementDof.localize(jproc);
  }

  void elem::LocalizeElementDofFromOneToOne( const unsigned& jproc, const unsigned& kproc ) {
    _elementDof.localize(jproc);
  }

  unsigned elem::GetElementDofIndex( const unsigned& iel, const unsigned& inode ) {  
    return _elementDof[iel][inode];
  };

  void elem::FreeLocalizedElementDof() {  
    _elementDof.clearLocalized();
  }

  void elem::ScatterElementNearFace() {
    _elementNearFace.scatter(_elementOffset);
  };

  void elem::LocalizeElementNearFaceFromOneToAll( const unsigned& jproc ) {
    _elementNearFace.localize(jproc);
  }

  void elem::FreeLocalizedElementNearFace() {
    _elementNearFace.clearLocalized();
  }

} //end namespace femus


