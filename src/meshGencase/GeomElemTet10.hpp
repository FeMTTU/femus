#ifndef __femus_meshGencase_GeomElemTet10_hpp__
#define __femus_meshGencase_GeomElemTet10_hpp__

//Class for
#include <cstdlib>
#include <iostream>

#include "GeomElemBase.hpp"


namespace femus {



class GeomElemTet10 : public GeomElemBase  {

public:
  
     GeomElemTet10();
     
    ~GeomElemTet10();
  
    unsigned int  get_dimension() const { return 3; };
    unsigned int n_nodes()        const { return 10; };
    unsigned int n_nodes_linear() const { return 4; };
    std::string   get_name_med()  const { return "T10"; };
    std::string   get_name_xdmf() const { return "Tetrahedron_10"; };
    std::vector<unsigned> get_face (const unsigned f) const { std::vector<unsigned> my_faces(_faces[f],_faces[f] + 6);  return my_faces; }; 
    
      float get_embedding_matrix(const uint,const uint,const uint);

         double get_prol(const uint /*j*/) {std::cout << "Tet10: no prolongation needed\n"; abort();};

private:
    
    static const unsigned _faces[4][6];
    
    
      static const float _embedding_matrix[8][10][10];   // (volume)

};


} //end namespace femus



#endif
