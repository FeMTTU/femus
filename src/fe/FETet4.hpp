#ifndef __fetet4_h__
#define __fetet4_h__

//Class for


#include "FEElemBase.hpp"


namespace femus {



class FETet4 : public FEElemBase  {

public:
  
     FETet4();
     
    ~FETet4();
  
      float get_embedding_matrix(const uint,const uint,const uint);

      static const float _embedding_matrix[8][4][4];   // (volume)

              double get_prol(const uint j) {return _Prol[j];};
    static const double _Prol[/*NNDS*/10*4/*NNDSL*/];

  
};


} //end namespace femus



#endif