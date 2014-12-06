#ifndef __fetri1_h__
#define __fetri1_h__

#include <cstdlib>
#include <iostream>


#include "FEElemBase.hpp"


namespace femus {



class FETri1 : public FEElemBase  {

public:
  
     FETri1();
     
    ~FETri1();
  
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[4][1][1];   // (volume)

                double get_prol(const uint /*j*/) {/*return _Prol[j];*/std::cout << "FETri1: no prolongation needed\n"; abort();};
//       static const double _Prol[/*NNDS*/6*3/*NNDSL*/];

};


} //end namespace femus



#endif