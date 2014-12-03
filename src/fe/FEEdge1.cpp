#include "FEEdge1.hpp"


namespace femus {

FEEdge1::FEEdge1(std::vector<GeomEl> geomel_in) : FEElemBase(geomel_in) {
	    
              _ndof[VV]=1;    
	      _ndof[BB]=1; //REMOVE
	  }
	  
    FEEdge1::~FEEdge1() { }

} //end namespace femus


 
 
 
 