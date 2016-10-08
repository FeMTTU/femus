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
      
    /////////////////////////////
    _elementDof1.resize(_nel,NVE[0][2]);
    _elementNearFace1.resize(_nel,NFC[0][1]);
    ///////////////////////////////////
    
    _elementDof = new unsigned * [_nel + 1];
    _elementNearFace = new int *[_nel + 1];

    _elementDofMemorySize = _nel * NVE[0][2];
    _elementNearFaceMemorySize = _nel * NFC[0][1];

    _elementDofMemory = new unsigned [_elementDofMemorySize];

    _elementNearFaceMemory = new int [_elementNearFaceMemorySize];
    for( unsigned i = 0; i < _elementNearFaceMemorySize; i++ )
      _elementNearFaceMemory[i] = -1;

    unsigned* pt_u = _elementDofMemory;
    int* pt_i = _elementNearFaceMemory;

    for( unsigned i = 0; i < _nel; i++ ) {
      _elementDof[i] = pt_u;
      pt_u += NVE[0][2];
      _elementNearFace[i] = pt_i;
      pt_i += NFC[0][1];
    }
    _elementDof[_nel] = pt_u;
    _elementNearFace[_nel] = pt_i;
    
    _elementDofIsScattered = false;
    _elementDofIsLocalizedFromJproc = false;
    _localElementDofMemorySize = 0;

    _elementNearFaceIsScattered = false;
    _elementNearFaceIsLocalizedFromJproc = false;
    _localElementNearFaceMemorySize = 0;
    
    /////////////////////////////
    
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
	tmpElDof[i][j] = _elementDof1[i][j];
      }
    }
    
    for(unsigned i=tmpElNearFace.begin();i<tmpElNearFace.end();i++){
      for(unsigned j=tmpElNearFace.begin(i);j<tmpElNearFace.end(i);j++){
	tmpElNearFace[i][j] = _elementNearFace1[i][j];
      }
    }
    
    _elementDof1 = tmpElDof;
    _elementNearFace1 = tmpElNearFace;
    
    tmpElDof.clear();
    tmpElNearFace.clear();
    
    // ***************************************************************************
    
    _elementDofMemorySize = _nelt[0] * NVE[0][2] + _nelt[1] * NVE[1][2] +
                            _nelt[2] * NVE[2][2] + _nelt[3] * NVE[3][2] +
                            _nelt[4] * NVE[4][2] + _nelt[5] * NVE[5][2];
    
    unsigned** tempElementDof = _elementDof;

    _elementDof = new unsigned* [_nel + 1];
    _elementDofMemory = new unsigned [_elementDofMemorySize];

    unsigned* ptElemDofMem = _elementDofMemory;
    for( unsigned iel = 0; iel < _nel; iel++ ) {
      _elementDof[iel] = ptElemDofMem;
      unsigned ielType = GetElementType( iel );
      for( unsigned j = 0; j < NVE[ielType][2]; j++ ) {
        _elementDof[iel][j] = tempElementDof[iel][j];
      }
      ptElemDofMem += NVE[ielType][2];
    }
    _elementDof[_nel] = ptElemDofMem;

    delete [] tempElementDof[0];
    delete [] tempElementDof;


    // *****************************************************************


    _elementNearFaceMemorySize = _nelt[0] * NFC[0][1] + _nelt[1] * NFC[1][1] +
                                 _nelt[2] * NFC[2][1] + _nelt[3] * NFC[3][1] +
                                 _nelt[4] * NFC[4][1] + _nelt[5] * NFC[5][1];

    int** tempElementNearFace = _elementNearFace;

    _elementNearFace = new int* [_nel + 1];
    _elementNearFaceMemory = new int [_elementNearFaceMemorySize];

    int* ptElemNearFaceMem = _elementNearFaceMemory;
    for( unsigned iel = 0; iel < _nel; iel++ ) {
      _elementNearFace[iel] = ptElemNearFaceMem;
      unsigned ielType = GetElementType( iel );
      for( unsigned j = 0; j < NFC[ielType][1]; j++ ) {
        _elementNearFace[iel][j] = tempElementNearFace[iel][j];
      }
      ptElemNearFaceMem += NFC[ielType][1];
    }
    _elementNearFace[_nel] = ptElemNearFaceMem;

    delete [] tempElementNearFace[0];
    delete [] tempElementNearFace;

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
    _elementDof1 = MyMatrix <unsigned> (rowSizeElDof);
    _elementNearFace1 = MyMatrix <int> (rowSizeElNearFace,-1);
    
    rowSizeElDof.clear();
    rowSizeElNearFace.clear();
    
    //std::cout << _elementDof1 << std::endl;
    //std::cout << _elementNearFace1 <<std::endl;
    
    // *************************
    _elementDof = new unsigned * [_nel + 1];
    _elementNearFace = new int * [_nel + 1];

    _elementDofMemorySize = 0;
    _elementNearFaceMemorySize = 0;
    for( unsigned i = 0; i < N_GEOM_ELS; i++ ) {
      _elementDofMemorySize += elc->GetRefinedElementTypeNumber( i ) * refindex * NVE[i][2];
      _elementNearFaceMemorySize += elc->GetRefinedElementTypeNumber( i ) * refindex * NFC[i][1];
    }
    
    for(unsigned isdom = 0; isdom < elc->_nprocs; isdom++) {
      elc->_elementType.localize(isdom);
      for(unsigned iel = elc->_elementOffset[isdom]; iel < elc->_elementOffset[isdom + 1]; iel++) {
	if( static_cast < short unsigned >( coarseAmrVector[iel] + 0.25 ) == 0 ) {
	  short unsigned type = elc->_elementType[iel];
	  _elementDofMemorySize += NVE[type][2];
	  _elementNearFaceMemorySize += NFC[type][1];
	}
      }
      elc->_elementType.clearLocalized();
    }

    _elementDofMemory = new unsigned [ _elementDofMemorySize ];
    _elementNearFaceMemory = new int [ _elementNearFaceMemorySize ];
    for( unsigned i = 0; i < _elementNearFaceMemorySize; i++ )
      _elementNearFaceMemory[i] = -1;
    
    int* ptElementNearFaceMemory = _elementNearFaceMemory;
    unsigned* ptElementDofMemory = _elementDofMemory;
    jel = 0;
    
    for(unsigned isdom = 0; isdom < elc->_nprocs; isdom++) {
      elc->_elementType.localize(isdom);
      for(unsigned iel = elc->_elementOffset[isdom]; iel < elc->_elementOffset[isdom + 1]; iel++) {
   	short unsigned elemt = elc->_elementType[iel];
	int increment = 1;
	if( static_cast < short unsigned >( coarseAmrVector[iel] + 0.25 ) == 1 ) {
	  increment = NRE[elemt];
	}
	for( unsigned j = 0; j < increment; j++ ) {
	  _elementDof[jel + j] = ptElementDofMemory;
	  ptElementDofMemory += NVE[elemt][2];
	  _elementNearFace[jel + j] = ptElementNearFaceMemory;
	  ptElementNearFaceMemory += NFC[ elemt ][1];
	}
	jel += increment;
      }
      elc->_elementType.clearLocalized();
    }
    _elementDof[_nel] = ptElementDofMemory;
    _elementNearFace[_nel] = ptElementNearFaceMemory;

    _elementDofIsScattered = false;
    _elementDofIsLocalizedFromJproc = false;
    _localElementDofMemorySize = 0;

    _elementNearFaceIsScattered = false;
    _elementNearFaceIsLocalizedFromJproc = false;
    _localElementNearFaceMemorySize = 0;
    
    // ***********************************
        
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
    //if( _level == 0 ) {
      MyVector <short unsigned> tempElementGroup;
      MyVector <short unsigned> tempElementMaterial;
      tempElementGroup = _elementGroup;
      tempElementMaterial = _elementMaterial;
      //_elementGroup = new short unsigned [_nel];
      //_elementMaterial = new short unsigned [_nel];
      for( unsigned iel = 0; iel < _nel; iel++ ) {
        _elementGroup[elementMapping[iel]]   = tempElementGroup[iel];
        _elementMaterial[elementMapping[iel]] = tempElementMaterial[iel] ;
      }
      tempElementGroup.clear();
      tempElementMaterial.clear();
    //}
    //END reordering _elementGroup and _elementMaterial

      
    //***********************
      
    
    MyMatrix <unsigned> tmpElDof = _elementDof1;
    for(unsigned i=_elementDof1.begin();i<_elementDof1.end();i++){
      for(unsigned j=_elementDof1.begin(i);j<_elementDof1.end(i);j++){
	_elementDof1[elementMapping [i]][j] = tmpElDof[i][j];
      }
    }
    tmpElDof.clear();
    
    
    MyMatrix <int> tmpElNearFace = _elementNearFace1;
    for(unsigned i=_elementNearFace1.begin();i<_elementNearFace1.end();i++){
      for(unsigned j=_elementNearFace1.begin(i);j<_elementNearFace1.end(i);j++){
	_elementNearFace1[elementMapping [i]][j] = tmpElNearFace[i][j];
      }
    }
    tmpElNearFace.clear();
    
          
    //BEGIN reordering _elementNearFace (rows)
    int** tempElementNearFace;
    tempElementNearFace = _elementNearFace;

    _elementNearFace = new int * [_nel + 1];
    _elementNearFaceMemory = new int [ _elementNearFaceMemorySize ];

    int* ptKel = _elementNearFaceMemory;

    for( unsigned iel = 0; iel < _nel; iel++ ) {
      _elementNearFace[iel] = ptKel;
      ptKel +=  NFC[_elementType[iel]][1];
    }
    _elementNearFace[_nel] = ptKel;

    for( unsigned iel = 0; iel < _nel; iel++ ) {
      for( unsigned iface = 0; iface < NFC[ _elementType[ elementMapping [iel] ]][1]; iface++ ) {
        _elementNearFace[elementMapping [iel]][iface] = tempElementNearFace[iel][iface];
      }
    }
    delete [] tempElementNearFace[0];
    delete [] tempElementNearFace;
    //END reordering _elementNearFace

    //BEGIN reordering _elementDof (rows)
    unsigned** tempElementDof = _elementDof;

    tempElementDof = _elementDof;

    _elementDof = new unsigned* [_nel + 1];
    _elementDofMemory = new unsigned [ _elementDofMemorySize ];

    unsigned* ptElementDof = _elementDofMemory;

    for( unsigned iel = 0; iel < _nel; iel++ ) {
      _elementDof[iel] = ptElementDof;
      ptElementDof +=  NVE[_elementType[iel]][2];
    }
    _elementDof[_nel] = ptElementDof;

    for( unsigned iel = 0; iel < _nel; iel++ ) {
      for( unsigned inode = 0; inode < NVE[_elementType[elementMapping [iel]]][2]; inode++ ) {
        _elementDof[elementMapping [iel]][inode] = tempElementDof[iel][inode];
      }
    }

    delete [] tempElementDof[0];
    delete [] tempElementDof;
    //END reordering OF _elementDof

    //********************************************************
    
    //BEGIN reordering _childElementDof (columns) on coarse level
    if( _level != 0 ) {
      for( unsigned i = 0; i < _coarseElem->_childElemMemorySize; i++ ) {
        _coarseElem->_childElemMemory[i] =  elementMapping[ _coarseElem->_childElemMemory[i]];
      }
    }
    //END reordering _childElementDof
  }


  void elem::ReorderMeshNodes( const std::vector < unsigned >& nodeMapping ) {
    for( unsigned i = 0; i < _elementDofMemorySize; i++ ) {
      _elementDofMemory[i] =  nodeMapping[ _elementDofMemory[i] ];
    }
    
    for( unsigned i = _elementDof1.begin(); i < _elementDof1.end(); i++ ) {
      for( unsigned j = _elementDof1.begin(i); j < _elementDof1.end(i); j++ ) {
	_elementDof1[i][j] =  nodeMapping[ _elementDof1[i][j] ];
      }
    }
    
    
  }

  elem::~elem() {

    delete [] _localElementDof;
    delete [] _localElementDofMemory;

    if( _elementDof != NULL ) {
      delete [] _elementDof;
      delete [] _elementDofMemory;
      _elementDof = NULL;
    }

    delete [] _localElementNearFace;
    delete [] _localElementNearFaceMemory;

    if( _elementNearFace != NULL ) {
      delete [] _elementNearFace;
      delete [] _elementNearFaceMemory;
      _elementNearFace == NULL;
    }

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
    _elementDof[iel][inode] = value;
    _elementDof1[iel][inode] = value;
  }

  /**
   * Return the local->global face node index
   **/
  unsigned elem::GetFaceVertexIndex( const unsigned& iel, const unsigned& iface, const unsigned& inode ) {
    unsigned value = _elementDof1[iel][ig[_elementType[iel]][iface][inode]];
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
    return _elementNearFace1[iel][iface];
    
    if (_elementNearFaceIsScattered ) {
      return _localElementNearFace[iel - _elementOffset[ _iproc ]][iface];
    }
    else if(_elementNearFaceIsLocalizedFromJproc){
      return _elementNearFace[iel- _elementOffset[ _jprocElementNearFaceIsLocalizedFrom ]][iface];
    }
    else{
      return _elementNearFace[iel][iface];
    }
  }

  int elem::GetBoundaryIndex( const unsigned& iel, const unsigned& iface ) {
    return  -( GetFaceElementIndex( iel, iface ) + 1 );
  }


  /**
   * Set the global adiacent-to-face element number
   **/
  void elem::SetFaceElementIndex( const unsigned& iel, const unsigned& iface, const int& value ) {
    _elementNearFace[iel][iface] = value;
    _elementNearFace1[iel][iface] = value;
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

    _elementDof1.scatter(_elementOffset);
    
    if( _localElementDofMemorySize == 0 ) {
      _localElementDofMemorySize = _elementDof[_elementOffset[ _iproc + 1]] - _elementDof[_elementOffset[ _iproc ]];

      _elementDofOffset = _elementDof[_elementOffset[ _iproc ]] - _elementDof[0];

      _localElementDofMemory = new unsigned [_localElementDofMemorySize];

      for( unsigned i = 0; i < _localElementDofMemorySize; i++ ) {
        _localElementDofMemory[i] =  *( _elementDof[_elementOffset[ _iproc ]] + i );
      }

      _localElementDof = new unsigned *[_elementOwned + 1];
      unsigned* pt_u = _localElementDofMemory;
      for( unsigned i = 0; i < _elementOwned; i++ ) {
        _localElementDof[i] = pt_u;
        pt_u += _elementDof[_elementOffset[ _iproc ] + ( i + 1 )] - _elementDof[_elementOffset[ _iproc ] + i];
      }
      _localElementDof[_elementOwned] = pt_u;
    }
    _elementDofIsScattered = true;
    _elementDofIsLocalizedFromJproc = false;

    if( _elementDof != NULL ) {
      delete [] _elementDofMemory;
      delete [] _elementDof;
      _elementDof = NULL;
    }

  };

  void elem::LocalizeElementDofFromOneToAll( const unsigned& jproc ) {

    _elementDof1.localize(jproc);
    
    if( _elementDof != NULL ) {
      delete [] _elementDofMemory;
      delete [] _elementDof;
      _elementDof = NULL;
    }

    unsigned jLocalElementDofMemorySize = _localElementDofMemorySize;
    MPI_Bcast( &jLocalElementDofMemorySize, 1, MPI_UNSIGNED, jproc, MPI_COMM_WORLD );

    _elementDofMemory = new unsigned [jLocalElementDofMemorySize];

    if( _iproc == jproc ) {
      for( unsigned i = 0; i < jLocalElementDofMemorySize; i++ ) {
        _elementDofMemory[i] = _localElementDofMemory[i];
      }
    }
    MPI_Bcast( _elementDofMemory, jLocalElementDofMemorySize, MPI_UNSIGNED, jproc, MPI_COMM_WORLD );

    unsigned bufferSize = _elementOffset[ jproc + 1 ] - _elementOffset[ jproc ];
    unsigned* buffer = new unsigned [ bufferSize ];

    if( _iproc == jproc ) {
      for( unsigned i = 0; i < _elementOwned; i++ ) {
        buffer[i] = _localElementDof[i + 1] - _localElementDof[i];
      }
    }

    MPI_Bcast( buffer , bufferSize, MPI_UNSIGNED, jproc, MPI_COMM_WORLD );

    _elementDof = new unsigned * [bufferSize + 1];
    unsigned* pt_u = _elementDofMemory;
    for( unsigned iel = 0; iel < bufferSize; iel++ ) {
      _elementDof[iel] = pt_u;
      pt_u += buffer[iel];
    }
    _elementDof[bufferSize] = pt_u;

    delete [] buffer;

    _elementDofIsScattered = false;
    _elementDofIsLocalizedFromJproc = true;
    _jprocElementDofIsLocalizedFrom = jproc;
  }


  void elem::LocalizeElementDofFromOneToOne( const unsigned& jproc, const unsigned& kproc ) {

    _elementDof1.localize(jproc);
    
    if( (_iproc == jproc || _iproc == kproc ) && jproc != kproc) {
      if( _iproc == kproc && _elementDof != NULL ) {
        delete [] _elementDofMemory;
        delete [] _elementDof;
        _elementDof = NULL;
      }

      unsigned jLocalElementDofMemorySize;
      if( _iproc == jproc ) {
        jLocalElementDofMemorySize = _localElementDofMemorySize;
        MPI_Send( &jLocalElementDofMemorySize, 1, MPI_UNSIGNED, kproc, 0, MPI_COMM_WORLD );
      }
      else {
        MPI_Recv( &jLocalElementDofMemorySize, 1, MPI_UNSIGNED, jproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      }

      if( _iproc == jproc ) {
        MPI_Send( _localElementDofMemory, jLocalElementDofMemorySize, MPI_UNSIGNED, kproc, 1 , MPI_COMM_WORLD );
      }
      else {
        _elementDofMemory = new unsigned [jLocalElementDofMemorySize];
        MPI_Recv( _elementDofMemory, jLocalElementDofMemorySize, MPI_UNSIGNED, jproc, 1,  MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      }

      unsigned bufferSize = _elementOffset[ jproc + 1 ] - _elementOffset[ jproc ];
      unsigned* buffer = new unsigned [ bufferSize ];

      if( _iproc == jproc ) {
        for( unsigned i = 0; i < _elementOwned; i++ ) {
          buffer[i] = _localElementDof[i + 1] - _localElementDof[i];
        }
        MPI_Send( buffer , bufferSize, MPI_UNSIGNED, kproc, 2, MPI_COMM_WORLD );
      }
      else {
        MPI_Recv( buffer , bufferSize, MPI_UNSIGNED, jproc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

        _elementDof = new unsigned * [bufferSize + 1];
        unsigned* pt_u = _elementDofMemory;
        for( unsigned iel = 0; iel < bufferSize; iel++ ) {
          _elementDof[iel] = pt_u;
          pt_u += buffer[iel];
        }
        _elementDof[bufferSize] = pt_u;
        _elementDofIsScattered = false;
        _elementDofIsLocalizedFromJproc = true;
        _jprocElementDofIsLocalizedFrom = jproc;
      }

      delete [] buffer;
    }
  }


  unsigned elem::GetElementDofIndex( const unsigned& iel, const unsigned& inode ) {
    
    return _elementDof1[iel][inode];
    
    if( _elementDofIsScattered ) {
      return _localElementDof[ iel - _elementOffset[ _iproc ]][inode];
    }
    else if( _elementDofIsLocalizedFromJproc ) {
      return _elementDof[ iel - _elementOffset[ _jprocElementDofIsLocalizedFrom ]][inode];
    }
    else {
      return _elementDof[iel][inode];
    }
  };

  void elem::FreeLocalizedElementDof() {
    
    _elementDof1.clearLocalized();
    
    _elementDofIsScattered = true;
    _elementDofIsLocalizedFromJproc = false;
    if( _elementDof != NULL ) {
      delete [] _elementDofMemory;
      delete [] _elementDof;
      _elementDof = NULL;
    }
  }

  void elem::ScatterElementNearFace() {

    _elementNearFace1.scatter(_elementOffset);
    
    if( _localElementNearFaceMemorySize == 0 ) {
      _localElementNearFaceMemorySize = _elementNearFace[_elementOffset[ _iproc + 1 ]] - _elementNearFace[_elementOffset[ _iproc ]];

      _elementNearFaceOffset = _elementNearFace[_elementOffset[ _iproc ]] - _elementNearFace[0];

      _localElementNearFaceMemory = new int [_localElementNearFaceMemorySize];

      for( unsigned i = 0; i < _localElementNearFaceMemorySize; i++ ) {
        _localElementNearFaceMemory[i] =  *( _elementNearFace[_elementOffset[ _iproc ]] + i );
      }

      _localElementNearFace = new int *[_elementOwned + 1];
      int* pt_i = _localElementNearFaceMemory;
      for( unsigned i = 0; i < _elementOwned; i++ ) {
        _localElementNearFace[i] = pt_i;
        pt_i += _elementNearFace[_elementOffset[ _iproc ] + ( i + 1 )] - _elementNearFace[_elementOffset[ _iproc ] + i];
      }
      _localElementNearFace[_elementOwned] = pt_i;
    }
    _elementNearFaceIsScattered = true;
    _elementNearFaceIsLocalizedFromJproc = false;

    if( _elementNearFace != NULL ) {
      delete [] _elementNearFaceMemory;
      delete [] _elementNearFace;
      _elementNearFace = NULL;
    }

  };

  void elem::LocalizeElementNearFaceFromOneToAll( const unsigned& jproc ) {
   
    _elementNearFace1.localize(jproc);
    
    if( _elementNearFace != NULL ) {
      delete [] _elementNearFaceMemory;
      delete [] _elementNearFace;
      _elementNearFace = NULL;
    }

    unsigned jLocalElementNearFaceMemorySize = _localElementNearFaceMemorySize;
    MPI_Bcast( &jLocalElementNearFaceMemorySize, 1, MPI_UNSIGNED, jproc, MPI_COMM_WORLD );

    _elementNearFaceMemory = new int [jLocalElementNearFaceMemorySize];

    if( _iproc == jproc ) {
      for( unsigned i = 0; i < jLocalElementNearFaceMemorySize; i++ ) {
        _elementNearFaceMemory[i] = _localElementNearFaceMemory[i];
      }
    }
    MPI_Bcast( _elementNearFaceMemory, jLocalElementNearFaceMemorySize, MPI_INT, jproc, MPI_COMM_WORLD );

    unsigned bufferSize = _elementOffset[ jproc + 1 ] - _elementOffset[ jproc ];
    unsigned* buffer = new unsigned [ bufferSize ];

    if( _iproc == jproc ) {
      for( unsigned i = 0; i < _elementOwned; i++ ) {
        buffer[i] = _localElementNearFace[i + 1] - _localElementNearFace[i];
      }
    }

    MPI_Bcast( buffer , bufferSize, MPI_UNSIGNED, jproc, MPI_COMM_WORLD );

    _elementNearFace = new int * [bufferSize + 1];
    int* pt_int = _elementNearFaceMemory;
    for( unsigned iel = 0; iel < bufferSize; iel++ ) {
      _elementNearFace[iel] = pt_int;
      pt_int += buffer[iel];
    }
    _elementNearFace[bufferSize] = pt_int;

    delete [] buffer;

    _elementNearFaceIsScattered = false;
    _elementNearFaceIsLocalizedFromJproc = true;
    _jprocElementNearFaceIsLocalizedFrom = jproc;
  }


  void elem::FreeLocalizedElementNearFace() {

    _elementNearFace1.clearLocalized();
    
    _elementNearFaceIsScattered = true;
    _elementNearFaceIsLocalizedFromJproc = false;

    if( _elementNearFace != NULL ) {
      delete [] _elementNearFaceMemory;
      delete [] _elementNearFace;
      _elementNearFace = NULL;
    }
  }

} //end namespace femus


