#include "FETri3.hpp"


namespace femus {




// =======================
FETri3::FETri3() : FEElemBase() { }
	  
// =======================
          FETri3::~FETri3() {    }

          
// =======================
      float FETri3::get_embedding_matrix(const uint a ,const uint b,const uint c )  {  return _embedding_matrix[a][b][c];  }
      

// =======================
//STATIC data member
	   const float FETri3::_embedding_matrix[4][3][3] =
{ // -------------------------------------------
  {// embedding matrix for child 0
    // 0    1    2  
    {1.0, 0.0, 0.0}, // 0
    {0.5, 0.5, 0.0}, // 1
    {0.5, 0.0, 0.5}  // 2
  },
  {// embedding matrix for child 1
    // 0    1    2  
    {0.5, 0.5, 0.0}, // 0
    {0.0, 1.0, 0.0}, // 1
    {0.0, 0.5, 0.5}  // 2
  },
  { // embedding matrix for child 2
    // 0    1    2  
    {0.5, 0.0, 0.5}, // 0
    {0.0, 0.5, 0.5}, // 1
    {0.0, 0.0, 1.0}  // 2
  },
  { // embedding matrix for child 3
    // 0    1    2  
    {0.5, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.5}, // 1
    {0.5, 0.0, 0.5}  // 2
  }
};


  const double FETri3::_Prol[6*3/*NNDS*NNDSL*/] = { 
   1.,0.,0.,
   0.,1.,0.,
   0.,0.,1.,
   0.5,0.5,0.,
   0.,0.5,0.5,
   0.5,0.,0.5};
   
 

} //end namespace femus


  