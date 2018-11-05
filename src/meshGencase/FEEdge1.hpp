#ifndef __femus_meshGencase_FEEdge1_hpp__
#define __femus_meshGencase_FEEdge1_hpp__

//Class for
#include <cstdlib>
#include <iostream>

#include "FEElemBase.hpp"


namespace femus {



class FEEdge1 : public FEElemBase  {

public:
  
     FEEdge1();
     
    ~FEEdge1();
  
    unsigned int  get_dimension() const { return 1; };
    unsigned int n_nodes()        const { return 1; };
    
     float get_embedding_matrix(const uint,const uint,const uint){std::cout << "Edge1: embedding matrix not implemented\n"; abort();};
//       static const float _embedding_matrix[][][];   // (volume)

      double get_prol(const uint /*j*/) {std::cout << "Edge1: prolongation not implemented \n"; abort();};

};


} //end namespace femus



#endif
