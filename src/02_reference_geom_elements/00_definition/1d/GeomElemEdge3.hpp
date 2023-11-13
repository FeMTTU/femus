#ifndef __femus_meshGencase_GeomElemEdge3_hpp__
#define __femus_meshGencase_GeomElemEdge3_hpp__


#include "GeomElemEdge.hpp"


namespace femus {



class GeomElemEdge3 : public GeomElemEdge  {

public:
  
      GeomElemEdge3() : GeomElemEdge() { };
        
    unsigned int n_nodes()        const { return 3; };
    
      std::vector<unsigned> get_nodes_of_face(const unsigned f) const { std::cout << "Not implemented FE" << __func__ << std::endl; abort(); };
    
// Refinement - BEGIN ===
public:
     float get_embedding_matrix(const uint,const uint,const uint){std::cout << "Edge3: embedding matrix not implemented\n"; abort();};

      double get_prol(const uint /*j*/) {std::cout << "Edge3: no prolongation needed\n"; abort();};

private:
    
//       static const float _embedding_matrix[][][];   // (volume)
      
// Refinement - END ===

      
// File names - BEGIN ===
    std::string   get_name_med()  const { return "SE3"; };
    std::string   get_name_xdmf() const { return "Edge_3"; };
// File names - END ===
      
};


} //end namespace femus



#endif
