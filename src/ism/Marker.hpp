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
    Marker( std::vector < double > x, const MarkerType &markerType, Mesh &mesh){
      _x = x;
      _markerType = markerType;
      _mesh = mesh;
    };
  
  std::vector < double > GetMarkerCoordinates(){ 
    return _x; 
  };  
  
  MarkerType GetMarkerType(){ 
    return _markerType; 
  };  
  
  // Gets 3 points in the plane and returns their orientation (the points must have only 2 components)
  int orientation2D(std::vector < double > P1, std::vector < double > P2, std::vector < double > P3){
    int val = (P2[1] - P1[1]) * (P3[0] - P2[0]) - (P2[0] - P1[0]) * (P3[1] - P2[1]);
    
    if (val == 0) return 0;  // colinear
    
    else return{(val > 0)? 1 : 2}; // clock or counterclock wise
  };
  
  // Given three colinear points P1, P2, P3 on the plane, the function checks if
  // point P2 lies on line segment P1P3
  bool onSegment2D(std::vector < double > P1, std::vector < double > P2, std::vector < double > P3) {
    if (P2[0] <= std::max(P1[0], P3[0]) && P2[0] >= std::min(P1[0], P3[0]) &&
        P2[1] <= std::max(P1[1], P3[1]) && P2[1] >= std::min(P1[1], P3[1]))
    return true;
 
    else return false;
  };
  
  // Tells if the line segments P1Q1 and P2Q2 intersect
  bool intersection2D( std::vector < double> P1, std::vector < double> Q1,
                       std::vector < double> P2, std::vector < double> Q2 ) {
  // Finds the four orientations needed for general and
  // special cases
   int o1 = orientation2D(P1, Q1, P2);
   int o2 = orientation2D(P1, Q1, Q2);
   int o3 = orientation2D(P2, Q2, P1);
   int o4 = orientation2D(P2, Q2, Q1);
 
    // General case
   if (o1 != o2 && o3 != o4)
   return true;
 
   // Special Cases
   // P1, Q1 and P2 are colinear and P2 lies on segment P1Q1
   else if (o1 == 0 && onSegment2D(P1, P2, Q1)) return true;
 
   // P1, Q1 and P2 are colinear and Q2 lies on segment P1Q1
   else if (o2 == 0 && onSegment2D(P1, Q2, Q1)) return true;
 
   // P2, Q2 and P1 are colinear and P1 lies on segment P2Q2
   else if (o3 == 0 && onSegment2D(P2, P1, Q2)) return true;
 
   // P2, Q2 and Q1 are colinear and Q1 lies on segment P2Q2
   else if (o4 == 0 && onSegment2D(P2, Q1, Q2)) return true;
 
   else return false; 
  };
  
  const unsigned GetMarkerElement(Marker &marker){
    std::vector < double > x = marker.GetMarkerCoordinates();
    std::vector < double > infty = x;
    infty[0] = 10000; 
    int count = 0;
    int markerElement = 0;
    while (count%2 == 0){ //we stop when count is odd so when count%2 == 1
    /*fare un loop sui nodi di ogni elemento del mesh, i nodi presi a due a due, chiamali inode e inode+1*/
    //if (intersection2D(inode, inode+1, x, infty)){
    // If the point x is colinear with line segment inode-(inode+1),
   // then check if it lies on segment. If it lies, return true,
   // otherwise false
    //if (orientation2D(inode, x, inode+1) == 0)
    // return onSegment2D(inode, x, inode+1);
    //count++;
  //}
    //markerElement = ielem
    }
    return markerElement;
  };
    
  private:
    std::vector < double > _x;
    MarkerType _markerType;
    Mesh _mesh;
    
  };
} //end namespace femus



#endif
