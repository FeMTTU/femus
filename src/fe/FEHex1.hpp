#ifndef __fehex1_h__
#define __fehex1_h__


#include <cstdlib>
#include <iostream>


#include "FEElemBase.hpp"


namespace femus {



class FEHex1 : public FEElemBase  {

public:
  
     FEHex1();
     
    ~FEHex1();
  
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[8/*NCHILDS*/][1/*NNDS*/][1/*NNDS*/];   // (volume)

      double get_prol(const uint /*j*/) {/*return _Prol[j];*/std::cout << "FEHex1: no prolongation needed\n"; abort();};
//        static const double _Prol[/*NNDS*/27*8/*NNDSL*/];

};


} //end namespace femus



#endif
