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
  
  //   class Quadrature {
//     
//   public: 
//     
//    Quadrature(unsigned int order, unsigned int geom_el);
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


// ************** LINE ***************
//first row-weights, second row: x-coordinates
const double GaussLine0[4][1]= {{2},
  {0}
};

const double GaussLine1[4][2]= {{1,1},
  {-0.57735026918963,0.57735026918963}
};

const double GaussLine2[4][3]= {{0.55555555555556,0.88888888888889,0.55555555555556},
  {-0.77459666924148,0,0.77459666924148}
};

const double GaussLine3[4][4]= {{0.34785484513745,0.65214515486255,0.65214515486255,0.34785484513745},
  {-0.86113631159405,-0.33998104358486,0.33998104358486,0.86113631159405}
};

const double GaussLine4[4][5]= {{0.23692688505619,0.47862867049937,0.56888888888889,0.47862867049937,0.23692688505619},
  {-0.90617984593866,-0.53846931010568,0,0.53846931010568,0.90617984593866}
};

//number of gauss point for the previous  LINE vector
const unsigned GaussPointsLine[5]= {1,2,3,4,5};





#endif