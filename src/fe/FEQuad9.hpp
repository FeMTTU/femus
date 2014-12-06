#ifndef __fequad9_h__
#define __fequad9_h__

#include <cstdlib>
#include <iostream>


#include "FEElemBase.hpp"


namespace femus {



class FEQuad9 : public FEElemBase  {

public:
  
     FEQuad9();
     
    ~FEQuad9();
  
    
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[4][9][9];   // (volume)

         double get_prol(const uint /*j*/) {std::cout << "Quad9: no prolongation needed\n"; abort();};

};


} //end namespace femus



#endif