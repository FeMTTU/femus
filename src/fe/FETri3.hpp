#ifndef __fetri3_h__
#define __fetri3_h__

//Class for


#include "FEElemBase.hpp"


namespace femus {



class FETri3 : public FEElemBase  {

public:
  
     FETri3();
     
    ~FETri3();
  
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[4][3][3];   // (volume)

                double get_prol(const uint j) {return _Prol[j];};
      static const double _Prol[/*NNDS*/6*3/*NNDSL*/];

};


} //end namespace femus



#endif