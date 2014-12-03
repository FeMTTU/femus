#include "FEEdge1.hpp"


namespace femus {

FEEdge1::FEEdge1(std::vector<GeomEl> geomel_in) : FEElemBase(geomel_in) {
	    
	      _name[VV]="Edge_1";
             _pname[VV]="Polyline";
              _ndof[VV]=1;    

	      _name[BB]="Edge_1_REMOVE"; 
	     _pname[BB]="Polyline_REMOVE"; 
	      _ndof[BB]=1; //REMOVE
	  }
	  
    FEEdge1::~FEEdge1() { }

} //end namespace femus


 
 
 
 