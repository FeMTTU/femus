#ifndef __femus_quadrature_Hexahedron_hpp__
#define __femus_quadrature_Hexahedron_hpp__


#include "quadrature_interface.hpp"



namespace femus {
    

  class hex_gauss {
      
  public:
      
    static const unsigned GaussPoints[ femus::Gauss::_NUM_QUAD_RULE_HEX_QUAD_LINE ];
    static const double *Gauss[ femus::Gauss::_NUM_QUAD_RULE_HEX_QUAD_LINE ];

  private:
    
    static const double Gauss0[4][1];
    static const double Gauss1[4][8];
    static const double Gauss2[4][27];
    static const double Gauss3[4][64];
    static const double Gauss4[4][125];
    static const double Gauss5[4][216];
  };
  
  
}


#endif
