#include "FEHex1.hpp"


namespace femus {




// =======================
FEHex1::FEHex1() : FEElemBase() {   }
	  
// =======================
          FEHex1::~FEHex1() {    }

          
// =======================
//this function is the same for all FE, but we had to do it identical for each one because of the fixed allocation
      float FEHex1::get_embedding_matrix(const uint a ,const uint b,const uint c )  {  return _embedding_matrix[a][b][c];  }
      

// =======================
//the embedding matrix picks every child of an element
// and, for every child node, gives the value of the FATHER SHAPE FNCS at the CHILD NODE
//STATIC data member
	   const float FEHex1::_embedding_matrix[8][1][1] =
      { // --------------------------------------
  // The 8 children of the Hex-type elements can be thought of as being
  // associated with the 8 vertices of the Hex.  Some of the children are
  // numbered the same as their corresponding vertex, while some are
  // not.  The children which are numbered differently have been marked
  // with ** in the comments below.
  {// embedding matrix for child 0 (child 0 is associated with vertex 0)
    //  0     1     2     3     4     5     6     7         //father----
    { 1.0} // 0  //child ----
  },
  { // embedding matrix for child 1 (child 1 is associated with vertex 1)
    //  0     1     2     3     4     5     6     7
    { 1.0}
  },
  { // embedding matrix for child 2 (child 2 is associated with vertex 3**)
    //  0      1    2     3     4     5     6     7
    { 1.0} // 0
  },
  {// embedding matrix for child 3 (child 3 is associated with vertex 2**)
    //  0      1    2     3     4     5     6     7
    { 1.0}
  },

  // embedding matrix for child 4 (child 4 is associated with vertex 4)
  {
    //  0      1    2     3     4     5     6     7
    { 1.0}
  },

  // embedding matrix for child 5 (child 5 is associated with vertex 5)
  {
    //  0      1    2     3     4     5     6     7
    { 1.0}
  },

  // embedding matrix for child 6 (child 6 is associated with vertex 7**)
  {
    //  0      1    2     3     4     5     6     7
    { 1.0}
  },

  // embedding matrix for child 7 (child 7 is associated with vertex 6**)
  {
    //  0      1    2     3     4     5     6     7
    {1.0}
  }
};



} //end namespace femus



