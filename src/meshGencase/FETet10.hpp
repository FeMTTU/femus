#ifndef __femus_meshGencase_FETet10_hpp__
#define __femus_meshGencase_FETet10_hpp__

//Class for
#include <cstdlib>
#include <iostream>

#include "GeomElemBase.hpp"


namespace femus {



class FETet10 : public GeomElemBase  {

public:
  
     FETet10();
     
    ~FETet10();
  
    unsigned int  get_dimension() const { return 3; };
    unsigned int n_nodes()        const { return 10; };
    std::string   get_name_med()  const { return "T10"; };
    std::string   get_name_xdmf() const { return "Tetrahedron_10"; };
    
      float get_embedding_matrix(const uint,const uint,const uint);

      static const float _embedding_matrix[8][10][10];   // (volume)
         double get_prol(const uint /*j*/) {std::cout << "Tet10: no prolongation needed\n"; abort();};

};


} //end namespace femus



#endif
