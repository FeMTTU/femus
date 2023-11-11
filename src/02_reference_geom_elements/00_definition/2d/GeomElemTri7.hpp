#ifndef __femus_meshGencase_GeomElemTri7_hpp__
#define __femus_meshGencase_GeomElemTri7_hpp__


#include "GeomElemBase.hpp"



namespace femus {



class GeomElemTri7 : public GeomElemBase  {

public:
  
     GeomElemTri7();
     
    ~GeomElemTri7();
  
    unsigned int  get_dimension() const { return 2; };
    unsigned int n_nodes_linear() const { return 3; };
    
    unsigned int n_nodes()        const { return 7; };
    
    std::string   get_name_med()  const { return "TR7"; };
    std::string   get_name_xdmf() const { abort(); /*return "Triangle_6";*/ };
    
    std::vector<unsigned> get_face (const unsigned f) const { std::vector<unsigned> my_faces(_faces[f],_faces[f] + 3);  return my_faces; }; 

    
private:
    
    static const unsigned _faces[3][3];
    

// Refinement - BEGIN ===
public:

  float get_embedding_matrix(const uint,const uint,const uint);


      double get_prol(const uint /*j*/) {std::cout << "Tri6: no prolongation needed\n"; abort();};

private:
    
    static const float _embedding_matrix[4][6][6];   // (volume)
// Refinement - END ===
      
};


} //end namespace femus



#endif
