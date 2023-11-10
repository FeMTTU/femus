#ifndef __femus_quadrature_Quadrangle_hpp__
#define __femus_quadrature_Quadrangle_hpp__

#include "quadrature_interface.hpp"

 namespace femus {
     
   
  class quad_gauss {
      
  public:
      
    static const unsigned GaussPoints[ femus::Gauss::_NUM_QUAD_RULE_HEX_QUAD_LINE];
    static const double *Gauss[ femus::Gauss::_NUM_QUAD_RULE_HEX_QUAD_LINE ];
    
  private:
    
    static const double Gauss0[3][1];
    static const double Gauss1[3][4];
    static const double Gauss2[3][9];
    static const double Gauss3[3][16];
    static const double Gauss4[3][25];
    static const double Gauss5[3][36];
  };
  
  
     
     
 }


#endif
