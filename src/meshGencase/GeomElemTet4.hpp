#ifndef __femus_meshGencase_GeomElemTet4_hpp__
#define __femus_meshGencase_GeomElemTet4_hpp__

//Class for


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemTet4 : public GeomElemBase  {

public:
  
     GeomElemTet4();
     
    ~GeomElemTet4();
  
    unsigned int  get_dimension() const { return 3; };
    unsigned int n_nodes()        const { return 4; };
    std::string   get_name_med()  const { return "TE4"; };
    std::string   get_name_xdmf() const { return "Tetrahedron"; };

    float get_embedding_matrix(const uint,const uint,const uint);

              double get_prol(const uint j) {return _Prol[j];};
    static const double _Prol[/*NNDS*/10*4/*NNDSL*/];

private:
    
      static const float _embedding_matrix[8][4][4];   // (volume)
    
};


} //end namespace femus



#endif
