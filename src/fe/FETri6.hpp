#ifndef __fetri6_h__
#define __fetri6_h__

//Class for
#include <cstdlib>
#include <iostream>

#include "FEElemBase.hpp"


namespace femus {



class FETri6 : public FEElemBase  {

public:
  
     FETri6();
     
    ~FETri6();
  
      float get_embedding_matrix(const uint,const uint,const uint);

      static const float _embedding_matrix[4][6][6];   // (volume)

      double get_prol(const uint /*j*/) {std::cout << "Tri6: no prolongation needed\n"; abort();};

};


} //end namespace femus



#endif