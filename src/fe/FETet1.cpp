#include "FETet1.hpp"


namespace femus {




// =======================
FETet1::FETet1(GeomEl* geomel_in) : FEElemBase(geomel_in) {
	    
	      _name[VV]="Tet_1";
	      _name[BB]="Tri_1"; 
             _pname[VV]="Tetrahedron";
	     _pname[BB]="Triangle"; 

              _ndof[VV]=1;    
	      _ndof[BB]=1;
	  }
	  
// =======================
          FETet1::~FETet1() {    }

          
// =======================
      float FETet1::get_embedding_matrix(const uint a ,const uint b,const uint c )  {  return _embedding_matrix[a][b][c];  }
      

// =======================
//STATIC data member
	   const float FETet1::_embedding_matrix[8][1][1] =
 { // --------------------------------------------
    { // embedding matrix for child 0
      // 0    1    2    3  
      {1.0}, // 0
    },
    { // embedding matrix for child 1
      // 0    1    2    3  
      {1.0}, // 0
    },
    {// embedding matrix for child 2
      // 0    1    2    3  
      {1.0}, // 0
    }, 
    {// embedding matrix for child 3
      // 0    1    2    3  
      {1.0}, // 0
    },
    {  // embedding matrix for child 4
      // 0    1    2    3  
      {1.0}, // 0
    },
    {// embedding matrix for child 5
      // 0    1    2    3  
      {1.0}, // 0
    },
    { // embedding matrix for child 6
      // 0    1    2    3  
      {1.0}, // 0
    },
    // embedding matrix for child 7
    {
      // 0    1    2    3  
      {1.0}, // 0
    }
  };
 


} //end namespace femus


   