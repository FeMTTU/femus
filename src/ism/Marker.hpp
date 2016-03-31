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
#include "vector"
#include "map"
#include "Mesh.hpp"


namespace femus {

  class Marker : public ParallelObject {
  public:
    Marker( std::vector < double > x, const MarkerType &markerType){
      _x = x;
      _markerType = markerType;
    };
  
  std::vector < double > GetMarkerCoordinates(){ 
    return _x; 
  };  
  
  MarkerType GetMarkerType(){ 
    return _markerType; 
  };  
    
  private:
    std::vector < double > _x;
    MarkerType _markerType;
    
  };
} //end namespace femus



#endif
