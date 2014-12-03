#include "FEEdge2.hpp"


namespace femus {

FEEdge2::FEEdge2(std::vector<GeomEl> geomel_in) : FEElemBase(geomel_in) {
	    
              _ndof[VV]=2;    
	      _ndof[BB]=2; //REMOVE
	  }
	  
    FEEdge2::~FEEdge2() { }

} //end namespace femus


 
 
 
 