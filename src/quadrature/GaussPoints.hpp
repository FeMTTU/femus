/*=========================================================================

 Program: FEMUS
 Module: FemusInit
 Authors: Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_quadrature_GaussPoints_hpp__
#define __femus_quadrature_GaussPoints_hpp__

#include <vector>
#include <string>

namespace femus {

  class hex_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[4][1];
    static const double Gauss1[4][8];
    static const double Gauss2[4][27];
    static const double Gauss3[4][64];
    static const double Gauss4[4][125];
  };
  
  
  class wedge_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[4][1];
    static const double Gauss1[4][8];
    static const double Gauss2[4][21];
    static const double Gauss3[4][52];
    static const double Gauss4[4][95];
  };  
  
  
  class tet_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[4][1];
    static const double Gauss1[4][5];
    static const double Gauss2[4][15];
    static const double Gauss3[4][31];
    static const double Gauss4[4][45];
  };

  class quad_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[3][1];
    static const double Gauss1[3][4];
    static const double Gauss2[3][9];
    static const double Gauss3[3][16];
    static const double Gauss4[3][25];
  };
  

  class tri_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[3][1];
    static const double Gauss1[3][4];
    static const double Gauss2[3][7];
    static const double Gauss3[3][13];
    static const double Gauss4[3][19];
  };
  
  
  class line_gauss {
  public:
    static const unsigned GaussPoints[6];
    static const double *Gauss[6];  
    static const double Gauss0[2][1];
    static const double Gauss1[2][2];
    static const double Gauss2[2][3];
    static const double Gauss3[2][4];
    static const double Gauss4[2][5];
    static const double Gauss5[2][6];
  };  
  
  class point_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[2][1];
    static const double Gauss1[2][1];
    static const double Gauss2[2][1];
    static const double Gauss3[2][1];
    static const double Gauss4[2][1];
  };  
  
  
  
  class Gauss {
     
  public:

    Gauss(const char *geom_elem, const char *order_gauss);
    
  inline const double *  GetGaussWeightsPointer() const {
    return GaussWeight;
  };
  
  inline const double *  GetGaussCoordinatePointer (const unsigned &k) const {
    return GaussWeight + (k+1) * GaussPoints;
  };
  
  
  inline const double  GetGaussWeight(const unsigned ig) const {
    return GaussWeight[ig];
  };
  
  inline const unsigned GetGaussPointsNumber() const {
      return GaussPoints;
  };     

  inline const std::string  GetGaussOrderString() const {
    return _order;
  };
  
  inline int  GetGaussOrderIdx() const {
    return gauss_order;
  };
  
  protected:
    
    int gauss_order;
    std::string _order;
    unsigned GaussPoints;  
    const double *GaussWeight;
   
  };
     

     
} //end namespace femus     


#endif
