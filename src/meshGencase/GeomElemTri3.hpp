#ifndef __femus_meshGencase_GeomElemTri3_hpp__
#define __femus_meshGencase_GeomElemTri3_hpp__

//Class for


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemTri3 : public GeomElemBase  {

public:
  
     GeomElemTri3();
     
    ~GeomElemTri3();
  
    unsigned int  get_dimension() const { return 2; };
    unsigned int n_nodes()        const { return 3; };
    std::string   get_name_med()  const { return "TR3"; };
    std::string   get_name_xdmf() const { return "Triangle"; };
    
      float get_embedding_matrix(const uint,const uint,const uint);

                double get_prol(const uint j) {return _Prol[j];};
      static const double _Prol[/*NNDS*/6*3/*NNDSL*/];

private:
    
      static const float _embedding_matrix[4][3][3];   // (volume)

};


} //end namespace femus



#endif
