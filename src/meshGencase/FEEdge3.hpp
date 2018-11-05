#ifndef __femus_meshGencase_FEEdge3_hpp__
#define __femus_meshGencase_FEEdge3_hpp__

//Class for
#include <cstdlib>
#include <iostream>

#include "FEElemBase.hpp"


namespace femus {



class FEEdge3 : public FEElemBase  {

public:
  
     FEEdge3();
     
    ~FEEdge3();
  
    unsigned int  get_dimension() const { return 1; };
    unsigned int n_nodes()        const { return 3; };
    std::string   get_name_med()  const { return "SE3"; };
    std::string   get_name_xdmf() const { return "Edge_3"; };
    
     float get_embedding_matrix(const uint,const uint,const uint){std::cout << "Edge3: embedding matrix not implemented\n"; abort();};
//       static const float _embedding_matrix[][][];   // (volume)

      double get_prol(const uint /*j*/) {std::cout << "Edge3: no prolongation needed\n"; abort();};

};


} //end namespace femus



#endif
