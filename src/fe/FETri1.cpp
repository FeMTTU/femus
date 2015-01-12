#include "FETri1.hpp"


namespace femus {




// =======================
FETri1::FETri1() : FEElemBase() {  }
	  
// =======================
          FETri1::~FETri1() {    }

          
// =======================
      float FETri1::get_embedding_matrix(const uint a ,const uint b,const uint c )  {  return _embedding_matrix[a][b][c];  }
      

// =======================
//STATIC data member
	   const float FETri1::_embedding_matrix[4][1][1] =
{ // -------------------------------------------
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


