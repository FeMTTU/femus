#ifndef __fetet1_h__
#define __fetet1_h__

#include <cstdlib>
#include <iostream>

#include"FEElemBase.hpp"


namespace femus {



class FETet1 : public FEElemBase  {

public:
  
     FETet1();
     
    ~FETet1();
  
      float get_embedding_matrix(const uint,const uint,const uint);

      static const float _embedding_matrix[8][1][1];   // (volume)

                 double get_prol(const uint /*j*/) {/*return _Prol[j];*/std::cout << "FETet1: no prolongation needed\n"; abort();};
//     static const double _Prol[/*NNDS*/10*4/*NNDSL*/];

  
};


} //end namespace femus



#endif