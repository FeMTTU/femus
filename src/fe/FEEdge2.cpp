#include "FEEdge2.hpp"


namespace femus {

FEEdge2::FEEdge2(std::vector<GeomEl> geomel_in) : FEElemBase(geomel_in) {
	    
	      _name[VV]="Edge_2";
             _pname[VV]="Polyline";
              _ndof[VV]=2;    

	      _name[BB]="Edge_2_REMOVE"; 
	     _pname[BB]="Polyline_REMOVE"; 
	      _ndof[BB]=2; //REMOVE
	  }
	  
    FEEdge2::~FEEdge2() { }

} //end namespace femus


 
 
 
 