#ifndef __femus_meshGencase_GeomElemTet4_hpp__
#define __femus_meshGencase_GeomElemTet4_hpp__




#include "GeomElemTet.hpp"


namespace femus {



class GeomElemTet4 : public GeomElemTet  {

public:
    
    unsigned int n_nodes()        const { return 4; };
    
      std::vector<unsigned> get_nodes_of_face(const unsigned f) const { std::cout << "Not implemented FE" << __func__ << std::endl; abort(); };
 
// Refinement - BEGIN ===
    float get_embedding_matrix(const uint,const uint,const uint);

              double get_prol(const uint j) {return _Prol[j];};

private:
      
    static const double _Prol[/*NNDS*/10*4/*NNDSL*/];
    
      static const float _embedding_matrix[8][4][4];   // (volume)
// Refinement - END ===
   
// File names - BEGIN ===
    std::string   get_name_med()  const { return "TE4"; };
    std::string   get_name_xdmf() const { return "Tetrahedron"; };
// File names - END ===
    
};


} //end namespace femus



#endif
