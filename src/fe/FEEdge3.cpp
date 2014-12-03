#include "FEEdge3.hpp"


namespace femus {

FEEdge3::FEEdge3(std::vector<GeomEl> geomel_in) : FEElemBase(geomel_in) {
	    
	      _name[VV]="Edge_3";
	      _name[BB]="Edge_3_REMOVE"; 
             _pname[VV]="Polyline";
	     _pname[BB]="Polyline_REMOVE"; 

	     _ndof[VV]=3;    
	     _ndof[BB]=3; //REMOVE
	  }
	  
    FEEdge3::~FEEdge3() { }

} //end namespace femus


 
 
 
 