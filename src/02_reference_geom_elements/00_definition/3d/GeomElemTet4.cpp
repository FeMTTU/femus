#include "GeomElemTet4.hpp"


namespace femus {


          
// =======================
      float GeomElemTet4::get_embedding_matrix(const uint a ,const uint b,const uint c )  {  return _embedding_matrix[a][b][c];  }
      

// =======================
//STATIC data member
	   const float GeomElemTet4::_embedding_matrix[8][4][4] =
 { // --------------------------------------------
    { // embedding matrix for child 0
      // 0    1    2    3  
      {1.0, 0.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0, 0.0}, // 1
      {0.5, 0.0, 0.5, 0.0}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    },
    { // embedding matrix for child 1
      // 0    1    2    3  
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 1.0, 0.0, 0.0}, // 1
      {0.0, 0.5, 0.5, 0.0}, // 2
      {0.0, 0.5, 0.0, 0.5}  // 3
    },
    {// embedding matrix for child 2
      // 0    1    2    3  
      {0.5, 0.0, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.0, 0.0, 1.0, 0.0}, // 2
      {0.0, 0.0, 0.5, 0.5}  // 3
    }, 
    {// embedding matrix for child 3
      // 0    1    2    3  
      {0.5, 0.0, 0.0, 0.5}, // 0
      {0.0, 0.5, 0.0, 0.5}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.0, 0.0, 0.0, 1.0}  // 3
    },
    {  // embedding matrix for child 4
      // 0    1    2    3  
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 0.5, 0.0, 0.5}, // 1
      {0.5, 0.0, 0.5, 0.0}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    },
    {// embedding matrix for child 5
      // 0    1    2    3  
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.5, 0.0, 0.5, 0.0}, // 2
      {0.0, 0.5, 0.0, 0.5}  // 3
    },
    { // embedding matrix for child 6
      // 0    1    2    3  
      {0.5, 0.0, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.0, 0.5, 0.0, 0.5}  // 3
    },
  
    // embedding matrix for child 7
    {
      // 0    1    2    3  
      {0.5, 0.0, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.0, 0.5}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    }
  };	   
  
  
   const double GeomElemTet4::_Prol[10*4/*NNDS*NNDSL*/] = { 
   1.,0.,0.,0.,
   0.,1.,0.,0.,
   0.,0.,1.,0.,
   0.,0.,0.,1.,
   0.5,0.5,0.,0.,
      0.,0.5,0.5,0.,
   0.5,0.,0.5,0.,
   0.5,0.,0.,0.5,
    0.,0.5,0.,0.5,
        0.,0.,0.5,0.5
   };  
  


} //end namespace femus


   