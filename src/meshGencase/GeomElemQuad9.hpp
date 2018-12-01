#ifndef __femus_meshGencase_FEQuad9_hpp__
#define __femus_meshGencase_FEQuad9_hpp__

#include <cstdlib>
#include <iostream>


#include "GeomElemBase.hpp"


namespace femus {



class FEQuad9 : public GeomElemBase  {

public:
  
     FEQuad9();
     
    ~FEQuad9();
  
    unsigned int  get_dimension() const { return 2; };
    unsigned int n_nodes()        const { return 9; };
    std::string   get_name_med()  const { return "QU9"; };
    std::string   get_name_xdmf() const { return "Quadrilateral_9"; };
    
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[4][9][9];   // (volume)

         double get_prol(const uint /*j*/) {std::cout << "Quad9: no prolongation needed\n"; abort();};

};


} //end namespace femus



#endif
