#ifndef __feedge2_h__
#define __feedge2_h__

//Class for
#include <cstdlib>
#include <iostream>

#include "FEElemBase.hpp"


namespace femus {



class FEEdge2 : public FEElemBase  {

public:
  
     FEEdge2();
     
    ~FEEdge2();
  
     float get_embedding_matrix(const uint,const uint,const uint){std::cout << "Edge2: embedding matrix not implemented\n"; abort();};
//       static const float _embedding_matrix[][][];   // (volume)

      double get_prol(const uint /*j*/) {std::cout << "Edge2: prolongation not implemented \n"; abort();};

};


} //end namespace femus



#endif