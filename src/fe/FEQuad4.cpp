#include "FEQuad4.hpp"


namespace femus {




// =======================
FEQuad4::FEQuad4()  : FEElemBase() {  }
	  
// =======================
          FEQuad4::~FEQuad4() {    }

          
// =======================
      float FEQuad4::get_embedding_matrix(const uint a ,const uint b,const uint c )  {  return _embedding_matrix[a][b][c];  }
      

// =======================
//STATIC data member


 const float FEQuad4::_embedding_matrix[4][4][4] = {
   
  // embedding matrix for child 0
  {
    // 0    1    2    3
    {1.0, 0.0, 0.0, 0.0}, // 0
    {0.5, 0.5, 0.0, 0.0}, // 1
    {.25, .25, .25, .25}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  },
  // embedding matrix for child 1
  {
    // 0    1    2    3
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 1.0, 0.0, 0.0}, // 1
    {0.0, 0.5, 0.5, 0.0}, // 2
    {.25, .25, .25, .25}  // 3
  },
  // embedding matrix for child 2
  {
    // 0    1    2    3
    {0.5, 0.0, 0.0, 0.5}, // 0
    {.25, .25, .25, .25}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.0, 0.0, 1.0}  // 3
  },
  // embedding matrix for child 3
  {
    // 0    1    2    3
    {.25, .25, .25, .25}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.0, 0.0, 1.0, 0.0}, // 2
    {0.0, 0.0, 0.5, 0.5}  // 3
  }
};



  const double FEQuad4::_Prol[9*4/*NNDS*NNDSL*/] =  { 
   1.,0.,0.,0.,
   0.,1.,0.,0.,
   0.,0.,1.,0.,
   0.,0.,0.,1.,
   0.5,0.5,0.,0.,
   0.,0.5,0.5,0.,
   0.,0.,0.5,0.5,
   0.5,0.,0.,0.5,
   0.25,0.25,0.25,0.25
    };


} //end namespace femus


