#ifndef __femus_meshGencase_GeomElemTri7_hpp__
#define __femus_meshGencase_GeomElemTri7_hpp__


#include "GeomElemTri.hpp"



namespace femus {



class GeomElemTri7 : public GeomElemTri  {

public:
    
       GeomElemTri7() : GeomElemTri() { };
    
    unsigned int n_nodes()        const { return 7; };
    
    
    std::vector<unsigned> get_nodes_of_face(const unsigned f) const { std::vector<unsigned> my_faces(_faces[f],_faces[f] + 3);  return my_faces; }; 

    
private:
    
    static const unsigned _faces[3][3];
    

// Refinement - BEGIN ===
public:

  float get_embedding_matrix(const uint,const uint,const uint);


      double get_prol(const uint /*j*/) {std::cout << "Tri6: no prolongation needed\n"; abort();};

private:
    
    static const float _embedding_matrix[4][6][6];   // (volume)
// Refinement - END ===


// File names - BEGIN ===
    std::string   get_name_med()  const { return "TR7"; };
    std::string   get_name_xdmf() const { abort(); /*return "Triangle_6";*/ };
// File names - END ===


};


} //end namespace femus



#endif
