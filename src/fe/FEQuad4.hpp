#ifndef __fequad4_h__
#define __fequad4_h__




#include "FEElemBase.hpp"


namespace femus {



class FEQuad4 : public FEElemBase  {

public:
  
     FEQuad4();
     
    ~FEQuad4();
  
    
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[4][4][4];   // (volume)

               double get_prol(const uint j) {return _Prol[j];};
     static const double _Prol[/*NNDS*/9*4/*NNDSL*/];
     
     //Shapes at quadrature points
     //The quadrature is known at runtime... so, we should SWITCH in each element with all the possible
     //quadrature rules... beh, intanto possiamo fare la porcata e fare una routine di lettura con tutti gli switch nella classe padre
     
};


} //end namespace femus



#endif