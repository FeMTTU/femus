#ifndef __femus_meshGencase_GeomElemEdge2_hpp__
#define __femus_meshGencase_GeomElemEdge2_hpp__


#include "GeomElemEdge.hpp"


namespace femus {



class GeomElemEdge2 : public GeomElemEdge  {

public:

      GeomElemEdge2() : GeomElemEdge() { };
      
    unsigned int n_nodes()        const { return 2; };
    
      std::vector<unsigned> get_nodes_of_face(const unsigned f) const { std::cout << "Not implemented FE" << __func__ << std::endl; abort(); };
    
// Refinement - BEGIN ===
public:

      float get_embedding_matrix(const uint,const uint,const uint){std::cout << "Edge2: embedding matrix not implemented\n"; abort();};

      double get_prol(const uint /*j*/) {std::cout << "Edge2: prolongation not implemented \n"; abort();};

private:
    
//       static const float _embedding_matrix[][][];   // (volume)
// Refinement - END ===

      
// File names - BEGIN ===
    std::string   get_name_med()  const { return "SE2"; };
    std::string   get_name_xdmf() const { return "Polyline"; };
// File names - END ===
    
      
};


} //end namespace femus



#endif
