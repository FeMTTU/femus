#ifndef __femus_meshGencase_FETri1_hpp__
#define __femus_meshGencase_FETri1_hpp__

#include <cstdlib>
#include <iostream>


#include "GeomElemBase.hpp"


namespace femus {



class FETri1 : public GeomElemBase  {

public:
  
     FETri1();
     
    ~FETri1();
  
    unsigned int  get_dimension() const { return 2; };
    unsigned int n_nodes()        const { return 1; };
    
      float get_embedding_matrix(const uint,const uint,const uint);
      static const float _embedding_matrix[4][1][1];   // (volume)

                double get_prol(const uint /*j*/) {/*return _Prol[j];*/std::cout << "FETri1: no prolongation needed\n"; abort();};
//       static const double _Prol[/*NNDS*/6*3/*NNDSL*/];

};


} //end namespace femus



#endif
