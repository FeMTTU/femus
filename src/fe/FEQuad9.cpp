#include "FEQuad9.hpp"


namespace femus {




// =======================
FEQuad9::FEQuad9() : FEElemBase() {   }
	  
// =======================
          FEQuad9::~FEQuad9() {    }

          
// =======================
      float FEQuad9::get_embedding_matrix(const uint a ,const uint b,const uint c )  {  return _embedding_matrix[a][b][c];  }
      

// =======================
//STATIC data member


 const float FEQuad9::_embedding_matrix[4][9][9] =
{
  // embedding matrix for child 0
  {
    //         0           1           2           3           4           5           6           7           8
    {    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
    {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 3
    {   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
    {    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 5
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,   0.750000 }, // 6
    {   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 7
    {   0.140625, -0.0468750,  0.0156250, -0.0468750,   0.281250, -0.0937500, -0.0937500,   0.281250,   0.562500 }  // 8
  },
  // embedding matrix for child 1
  {
    //         0           1           2           3           4           5           6           7           8
    {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
    {    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 3
    {  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
    {    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 5
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,   0.750000 }, // 6
    {    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 7
    { -0.0468750,   0.140625, -0.0468750,  0.0156250,   0.281250,   0.281250, -0.0937500, -0.0937500,   0.562500 }  // 8
  },
  // embedding matrix for child 2
  {
    //         0           1           2           3           4           5           6           7           8
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 0
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 1
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,   0.750000 }, // 4
    {    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 5
    {    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 6
    {  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 7
    { -0.0468750,  0.0156250, -0.0468750,   0.140625, -0.0937500, -0.0937500,   0.281250,   0.281250,   0.562500 }  // 8
  },
  // embedding matrix for child 3
  {
    //         0           1           2           3           4           5           6           7           8
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 0
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 1
    {    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 3
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,   0.750000 }, // 4
    {    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 5
    {    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 6
    {    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 7
    {  0.0156250, -0.0468750,   0.140625, -0.0468750, -0.0937500,   0.281250,   0.281250, -0.0937500,   0.562500 }  // 8
  }
};



} //end namespace femus


