#ifndef __fehex27_h__
#define __fehex27_h__

#include <cstdlib>
#include <iostream>
//Class for the 3D Hex27


#include "FEElemBase.hpp"
#include "ElemType.hpp"

namespace femus {



class FEHex27 : public FEElemBase {

public:
  
     FEHex27();
     
    ~FEHex27();
  
      float get_embedding_matrix(const uint,const uint,const uint);

      static const float _embedding_matrix[8/*NCHILDS*/][27/*NNDS*/][27/*NNDS*/];   // (volume)

       double get_prol(const uint /*j*/) {std::cout << "Hex27: no prolongation needed\n"; abort();};


};


} //end namespace femus



#endif


//We could do the template over the number of dofs, as a non-type template parameter...
//but then if you want to templatize you still need the NUMBER of CHILDREN,
//so you should do two nontype templates


//clearly, since this is a SPECIFIC Finite element, some allocations must be EXPLICIT
//the fact of being STATIC is for other reasons, not for the specific fact i think...

//The number of CHILDS is more related to the  'dimension' and the GEOMETRY, it is 8 both for Hex8 and Hex27

//Prol is the prolongation of a FE over a GEOM, so Quad4 on Quad9, Quad8 on Quad9.
//so, is it to be associated to 