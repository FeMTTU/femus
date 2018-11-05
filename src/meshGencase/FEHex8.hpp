#ifndef __femus_meshGencase_FEHex8_hpp__
#define __femus_meshGencase_FEHex8_hpp__




#include "FEElemBase.hpp"


namespace femus {



class  FEHex8 : public FEElemBase  {

public:
  
     FEHex8();
     
    ~FEHex8();
  
    unsigned int  get_dimension() const { return 3; };
    unsigned int n_nodes()        const { return 8; };
    std::string   get_name_med()  const { return "HE8"; };
    std::string   get_name_xdmf() const { return "Hexahedron"; };
    
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[8/*NCHILDS*/][8/*NNDS*/][8/*NNDS*/];   // (volume)

       double get_prol(const uint j) {return _Prol[j];};
       static const double _Prol[/*NNDS*/27*8/*NNDSL*/];
};


} //end namespace femus



#endif
