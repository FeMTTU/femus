#ifndef __femus_meshGencase_GeomElemEdge2_hpp__
#define __femus_meshGencase_GeomElemEdge2_hpp__

//Class for
#include <cstdlib>
#include <iostream>

#include "GeomElemBase.hpp"


namespace femus {



class GeomElemEdge2 : public GeomElemBase  {

public:
  
     GeomElemEdge2();
     
    ~GeomElemEdge2();
  
    unsigned int  get_dimension() const { return 1; };
    unsigned int n_nodes()        const { return 2; };
    std::string   get_name_med()  const { return "SE2"; };
    std::string   get_name_xdmf() const { return "Polyline"; };
    
     float get_embedding_matrix(const uint,const uint,const uint){std::cout << "Edge2: embedding matrix not implemented\n"; abort();};

      double get_prol(const uint /*j*/) {std::cout << "Edge2: prolongation not implemented \n"; abort();};

private:
    
//       static const float _embedding_matrix[][][];   // (volume)

};


} //end namespace femus



#endif
