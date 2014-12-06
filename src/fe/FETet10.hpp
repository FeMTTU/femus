#ifndef __fetet10_h__
#define __fetet10_h__

//Class for
#include <cstdlib>
#include <iostream>

#include "FEElemBase.hpp"


namespace femus {



class FETet10 : public FEElemBase  {

public:
  
     FETet10();
     
    ~FETet10();
  
      float get_embedding_matrix(const uint,const uint,const uint);

      static const float _embedding_matrix[8][10][10];   // (volume)
         double get_prol(const uint /*j*/) {std::cout << "Tet10: no prolongation needed\n"; abort();};

};


} //end namespace femus



#endif