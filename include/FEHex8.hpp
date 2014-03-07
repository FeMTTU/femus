#ifndef __fehex8_h__
#define __fehex8_h__




#include "FEElemBase.hpp"

class  FEHex8 : public FEElemBase  {

public:
  
     FEHex8(GeomEl* geomel_in);
     
    ~FEHex8();
  
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[8/*NCHILDS*/][8/*NNDS*/][8/*NNDS*/];   // (volume)

       double get_prol(const uint j) {return _Prol[j];};
       static const double _Prol[/*NNDS*/27*8/*NNDSL*/];
};

#endif