#ifndef __femus_quadrature_Wedge_hpp__
#define __femus_quadrature_Wedge_hpp__


namespace femus {
  

  class wedge_gauss {
      
  public:
      
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];
    
  private:
    
    static const double Gauss0[4][1];
    static const double Gauss1[4][8];
    static const double Gauss2[4][21];
    static const double Gauss3[4][52];
    static const double Gauss4[4][95];
  };  

  
}


#endif
