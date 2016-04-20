/*=========================================================================

 Program: FEMuS
 Module: Marker
 Authors: Eugenio Aulisa and Giacomo Capodaglio

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_ism_Marker_hpp__
#define __femus_ism_Marker_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MarkerTypeEnum.hpp"
#include "ParallelObject.hpp"
#include "Mesh.hpp"

#include "vector"
#include "map"
#include "Mesh.hpp"


namespace femus {

  class Marker : public ParallelObject {
  public:
    Marker( std::vector < double > x, const MarkerType &markerType, Mesh *mesh){
      _x = x;
      _markerType = markerType;
      _mesh = mesh;
      GetElement();
    };
  
  std::vector < double > GetMarkerCoordinates(){ 
    return _x; 
  };  
  
  void GetElement();
  
  MarkerType GetMarkerType(){ 
    return _markerType; 
  };  
    
  private:
    
    int GetNextElement2D(const unsigned &dim, const int &iel,  const int &kel);
    int GetNextElement3D(const unsigned &dim, const int &iel,  const int &kel);
    
    std::vector < double > _x;
    MarkerType _markerType;
    const Mesh * _mesh;
    unsigned _elem;
  };
} //end namespace femus



#endif
