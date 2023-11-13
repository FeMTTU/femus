#ifndef __femus_meshGencase_GeomElemQuad4_hpp__
#define __femus_meshGencase_GeomElemQuad4_hpp__



#include "GeomElemQuad.hpp"


namespace femus {



class GeomElemQuad4 : public GeomElemQuad  {

public:
  
      GeomElemQuad4() : GeomElemQuad() { };
  
        
    unsigned int n_nodes()        const { return 4; };
    
      std::vector<unsigned> get_nodes_of_face(const unsigned f) const { std::cout << "Not implemented FE" << __func__ << std::endl; abort(); };

// Refinement - BEGIN ===
public:
      
    float get_embedding_matrix(const uint,const uint,const uint);

               double get_prol(const uint j) {return _Prol[j];};
     
private:
      
      static const double _Prol[/*NNDS*/9*4/*NNDSL*/];
   
      static const float _embedding_matrix[4][4][4];   // (volume)
// Refinement - END ===


// File names - BEGIN ===
    std::string   get_name_med()  const { return "QU4"; };
    std::string   get_name_xdmf() const { return "Quadrilateral"; };
// File names - END ===


};


} //end namespace femus



#endif
