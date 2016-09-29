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
    Marker( std::vector < double > x, const MarkerType &markerType, Mesh *mesh, const unsigned & solType, const bool &debug = false){
      _x = x;
      _markerType = markerType;
      _mesh = mesh;
      _solType = solType;
      
      GetElement(debug);
      
//       _xi.resize( _mesh->GetDimension() );
//       short unsigned elemType = _mesh->GetElementType(_elem);
//       for(int k = 0; k < _mesh->GetDimension(); k++) {
//         _xi[k] = _initialGuess[elemType][k];
//       }
//       for(int itype = 0; itype <= _solType; itype++) {
//         InverseMapping(_elem, itype,_x,_xi);
//       }
//       for(int k = 0; k < _mesh->GetDimension(); k++) {
// 	std::cout<< _xi[k] <<" ";
//       }
//       std::cout << std::endl;
      
     
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
    
    unsigned GetNextElement2D(const unsigned &iel);
    unsigned GetNextElement3D(const unsigned &iel);
    int FastForward(const unsigned &currentElem);
     
    std::vector < double > _x;
    std::vector < double > _xi;
    unsigned _solType;
    MarkerType _markerType;
    const Mesh * _mesh;
    unsigned _elem;
    
    static const double _initialGuess[6][3];

  };
} //end namespace femus



#endif
