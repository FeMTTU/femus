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
    Marker( std::vector < double > x, const MarkerType &markerType, Mesh *mesh,const bool &debug = false){
      _x = x;
      _markerType = markerType;
      _mesh = mesh;
      GetElement(debug);
    };
  
  std::vector < double > GetMarkerCoordinates(){ 
    return _x; 
  };  
  

  void GetElement(const bool & debug = false);
  
  
  MarkerType GetMarkerType(){ 
    return _markerType; 
  };  
  
  void InverseMappingTEST(std::vector< double > &x);
  
  std::vector<double> GetPosition(std::vector<double> (*f)(std::vector<double>), int n, double T);
    
  private:
    
    
    std::vector< double > InverseMapping(const unsigned &currentElem, const unsigned &solutionType, const std::vector< double > &x);
    void InverseMapping(const unsigned &iel, const unsigned &solType, 
			      const std::vector< double > &x, std::vector< double > &xi);
    std::vector< double > InverseMappingTri(const unsigned &currentElem, const unsigned &solutionType,  std::vector< double > &x);
    std::vector< double > InverseMappingHex(const unsigned &currentElem, const unsigned &solutionType,  std::vector< double > &x);
    std::vector< double > InverseMappingTet(const unsigned &currentElem, const unsigned &solutionType,  std::vector< double > &x);
    std::vector< double > InverseMappingWedge(const unsigned &currentElem, const unsigned &solutionType, std::vector< double > &x);
    
    unsigned GetNextElement2D(const unsigned &iel,  const unsigned &kel);
    unsigned GetNextElement3D(const unsigned &iel,  const unsigned &kel);
    unsigned GetNextElement2DOLD(const unsigned &dim, const unsigned &iel,  const unsigned &kel);
    unsigned GetNextElement3DOLD(const unsigned &dim, const unsigned &iel,  const unsigned &kel);
    
    bool SPDCheck2D(std::vector< std::vector <double> > A);
    bool SPDCheck3D(std::vector< std::vector <double> > A);
     
     
    std::vector < double > _x;
    MarkerType _markerType;
    const Mesh * _mesh;
    unsigned _elem;

  };
} //end namespace femus



#endif
