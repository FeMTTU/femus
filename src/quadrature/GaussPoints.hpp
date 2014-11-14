/*=========================================================================

 Program: FEMUS
 Module: FemTTUInit
 Authors: Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __gauss_points_hpp__
#define __gauss_points_hpp__

#include <vector>

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
    static const double Gauss0[4][1];
    static const double Gauss1[4][4];
    static const double Gauss2[4][9];
    static const double Gauss3[4][16];
    static const double Gauss4[4][25];
  };
  

  class tri_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[4][1];
    static const double Gauss1[4][4];
    static const double Gauss2[4][7];
    static const double Gauss3[4][13];
    static const double Gauss4[4][19];
  };
  
  
  class line_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[4][1];
    static const double Gauss1[4][2];
    static const double Gauss2[4][3];
    static const double Gauss3[4][4];
    static const double Gauss4[4][5];
  };  
  
  
  
  //   class Gauss {
//     
//   public: 
//     
//    Gauss(unsigned int order, unsigned int geom_el);
//    unsigned get_qpoints() { return n_qpoints;}; 
//    unsigned get_order()   { return order;}; 
// 
//   private:      
//     unsigned order;
//     unsigned n_qpoints;
//     std::vector<double>  weights;
//     std::vector< std::vector<double> >  qpoints;
// 
//   };



#endif