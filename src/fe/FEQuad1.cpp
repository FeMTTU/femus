#include "FEQuad1.hpp"


namespace femus {




// =======================
FEQuad1::FEQuad1(std::vector<GeomEl> geomel_in) : FEElemBase(geomel_in) {
  
  _name[VV]="Quad_1";
  _name[BB]="Edge_1"; 
  _pname[VV]="Quadrilateral"; //TODO
  _pname[BB]="Polyline";      //TODO
    
   _ndof[VV]=1;    
   _ndof[BB]=1;
	  }
	  
// =======================
          FEQuad1::~FEQuad1() {    }

          
// =======================
      float FEQuad1::get_embedding_matrix(const uint a ,const uint b,const uint c )  {  return _embedding_matrix[a][b][c];  }
      

// =======================
//STATIC data member


 const float FEQuad1::_embedding_matrix[4][1][1] = {
   
  // embedding matrix for child 0
  {
    // 0  
    {1.0}, // 0
  },
  // embedding matrix for child 1
  {
    // 0  
    {1.0}, // 0
  },
  // embedding matrix for child 2
  {
    // 0  
    {1.0}, // 0
  },
  // embedding matrix for child 3
  {
    // 0  
    {1.0}, // 0
  }
};


} //end namespace femus


