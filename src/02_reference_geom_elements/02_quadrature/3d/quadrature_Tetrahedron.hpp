#ifndef __femus_quadrature_Tetrahedron_hpp__
#define __femus_quadrature_Tetrahedron_hpp__



namespace femus {


    
  class tet_gauss {
      
  public:
      
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];
    
  private:
    
    static const double Gauss0[4][1];
    static const double Gauss1[4][5];
    static const double Gauss2[4][15];
    static const double Gauss3[4][31];
    static const double Gauss4[4][45];
  };

  
}



#endif
