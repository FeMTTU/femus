#ifndef __femus_meshGencase_GeomElemHex8_hpp__
#define __femus_meshGencase_GeomElemHex8_hpp__




#include "GeomElemHex.hpp"


namespace femus {



class  GeomElemHex8 : public GeomElemHex  {

public:
  
       GeomElemHex8() : GeomElemHex() { };

    
    unsigned int n_nodes()        const { return 8; };
    
    
      std::vector<unsigned> get_nodes_of_face(const unsigned f) const { std::cout << "Not implemented FE" << __func__ << std::endl; abort(); };
    
    
// Refinement - BEGIN ===
public:

      float get_embedding_matrix(const uint,const uint,const uint);

       double get_prol(const uint j) {return _Prol[j];};
       
private:
      
       static const double _Prol[/*NNDS*/27*8/*NNDSL*/];
    
      static const float _embedding_matrix[8/*NCHILDS*/][8/*NNDS*/][8/*NNDS*/];   // (volume)
// Refinement - END ===


// File names - BEGIN ===
    std::string   get_name_med()  const { return "HE8"; };
    std::string   get_name_xdmf() const { return "Hexahedron"; };
// File names - END ===
    
       
};


} //end namespace femus



#endif
