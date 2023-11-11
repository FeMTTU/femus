#ifndef __femus_meshGencase_GeomElemTri3_hpp__
#define __femus_meshGencase_GeomElemTri3_hpp__



#include "GeomElemTri.hpp"


namespace femus {



class GeomElemTri3 : public GeomElemTri  {

public:
    
    
    unsigned int n_nodes()        const { return 3; };
    
    
      std::vector<unsigned> get_nodes_of_face(const unsigned f) const { std::cout << "Not implemented FE" << __func__ << std::endl; abort(); };
    
    
// Refinement - BEGIN ===
public:
  
       float get_embedding_matrix(const uint,const uint,const uint);

                double get_prol(const uint j) {return _Prol[j];};
      static const double _Prol[/*NNDS*/6*3/*NNDSL*/];

private:
    
      static const float _embedding_matrix[4][3][3];   // (volume)
// Refinement - END ===


// File names - BEGIN ===
    std::string   get_name_med()  const { return "TR3"; };
    std::string   get_name_xdmf() const { return "Triangle"; };
// File names - END ===

};


} //end namespace femus



#endif
