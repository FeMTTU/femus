#ifndef __femus_quadrature_Point_hpp__
#define __femus_quadrature_Point_hpp__

namespace femus {


  class point_gauss {
      
  public:
      
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];
    
  private:
    
    static const double Gauss0[2][1];
    static const double Gauss1[2][1];
    static const double Gauss2[2][1];
    static const double Gauss3[2][1];
    static const double Gauss4[2][1];
  };  
  
}


#endif
