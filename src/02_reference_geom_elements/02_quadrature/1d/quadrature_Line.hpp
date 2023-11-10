#ifndef __femus_quadrature_Line_hpp__
#define __femus_quadrature_Line_hpp__



#include "quadrature_interface.hpp"


namespace femus {

     
  class line_gauss {
      
  public:
      
    static const unsigned GaussPoints[ femus::Gauss::_NUM_QUAD_RULE_HEX_QUAD_LINE ];
    static const double *Gauss[ femus::Gauss::_NUM_QUAD_RULE_HEX_QUAD_LINE ];
    
  private:
    
    static const double Gauss0[2][1];
    static const double Gauss1[2][2];
    static const double Gauss2[2][3];
    static const double Gauss3[2][4];
    static const double Gauss4[2][5];
    static const double Gauss5[2][6];
  };  
   
}


#endif
