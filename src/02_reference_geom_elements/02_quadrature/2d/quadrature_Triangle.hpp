#ifndef __femus_quadrature_Triangle_hpp__
#define __femus_quadrature_Triangle_hpp__


 namespace femus {

     
     class tri_gauss {
         
  public:
      
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];
    
  private:
    
    static const double Gauss0[3][1];
    static const double Gauss1[3][4];
    static const double Gauss2[3][7];
    static const double Gauss3[3][13];
    static const double Gauss4[3][19];
  };
 
  
 }

 
#endif
