#ifndef __femus_meshGencase_GeomElemQuad9_hpp__
#define __femus_meshGencase_GeomElemQuad9_hpp__


#include "GeomElemQuad.hpp"


namespace femus {



class GeomElemQuad9 : public GeomElemQuad  {

public:
  
  
    
    unsigned int n_nodes()        const { return 9; };
    
    std::vector<unsigned> get_nodes_of_face(const unsigned f) const { std::vector<unsigned> my_faces(_faces[f],_faces[f] + 3);  return my_faces; }; 
    
    
private:
    
    static const unsigned _faces[4][3];
    
// Refinement - BEGIN ===
public:
  
      float get_embedding_matrix(const uint,const uint,const uint);

         double get_prol(const uint /*j*/) {std::cout << "Quad9: no prolongation needed\n"; abort();};

private:
     static const float _embedding_matrix[4][9][9];   // (volume)

// Refinement - END ===


// File names - BEGIN ===
    std::string   get_name_med()  const { return "QU9"; };
    std::string   get_name_xdmf() const { return "Quadrilateral_9"; };
// File names - END ===
     
};


} //end namespace femus



#endif
