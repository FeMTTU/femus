#ifndef __feedge3_h__
#define __feedge3_h__

//Class for
#include <cstdlib>
#include <iostream>

#include "FEElemBase.hpp"


namespace femus {



class FEEdge3 : public FEElemBase  {

public:
  
     FEEdge3();
     
    ~FEEdge3();
  
     float get_embedding_matrix(const uint,const uint,const uint){std::cout << "Edge3: embedding matrix not implemented\n"; abort();};
//       static const float _embedding_matrix[][][];   // (volume)

      double get_prol(const uint /*j*/) {std::cout << "Edge3: no prolongation needed\n"; abort();};

};


} //end namespace femus



#endif