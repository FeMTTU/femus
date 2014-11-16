/*=========================================================================

  Program: FEMUS
  Module: Triangle
  Authors: Eugenio Aulisa
 
  Copyright (c) FEMTTU
  All rights reserved. 

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Basis.hpp"


namespace femus {

  // triangle const vectors
  const double tri_lag::X[15][2]= {
    {0, 0},      {1, 0},      {0, 1},
    {0.5, 0},    {0.5, 0.5},  {0, 0.5},
    {0.25,0},    {0.25,0.25}, {0,0.25},
    {0.75,0.25}, {0.5,0.25},  {0.75,0},
    {0,0.75},    {0.25,0.5},  {0.25,0.75}
  };

  const int tri_lag::IND[6][2]= {
    {0, 0},{2, 0},{0, 2},
    {1, 0},{1, 1},{0, 1}
  };

  const int tri_lag::KVERT_IND[15][2]= {
    {0,0},{1,0},{2,0},
    {0,1},{1,1},{2,1},
    {0,3},{0,4},{0,5},
    {1,3},{1,4},{1,5},
    {2,3},{2,4},{2,5}
  };
  
  const double tri_const::X[12][2]={ 
    {0.166666666667,0.166666666667},{0.666666666667,0.166666666667},{0.166666666667,0.666666666667},{0.333333333333,0.333333333333},
    {0.166666666667,0.166666666667},{0.666666666667,0.166666666667},{0.166666666667,0.666666666667},{0.333333333333,0.333333333333},
    {0.166666666667,0.166666666667},{0.666666666667,0.166666666667},{0.166666666667,0.666666666667},{0.333333333333,0.333333333333},
  };
  
  const int tri_const::IND[3][2]= {{1, 0},{0,1},{0,0}};

  const int tri_const::KVERT_IND[12][2]={ 
    {0,0},{1,0},{2,0},{3,0},
    {0,1},{1,1},{2,1},{3,1},
    {0,2},{1,2},{2,2},{3,2}
  };
   
  
  double tripwl::eval_phi(const int *I,const double* x) const {
    return 1.*I[0]+x[0]*(1.-I[0])*I[1]+x[1]*(1.-I[0])*(1.-I[1]);
  }

  double tripwl::eval_dphidx(const int *I,const double* x) const {
    return (1.-I[0])*I[1];
  }

  double tripwl::eval_dphidy(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1]);
  }
  
  
  
  double tri1::eval_phi(const int *I,const double* x) const {
    return triangle1(x[0],x[1],I[0],I[1]);
  }

  double tri1::eval_dphidx(const int *I,const double* x) const {
    return dtriangle1dx(x[0],x[1],I[0],I[1]);
  }

  double tri1::eval_dphidy(const int *I,const double* x) const {
    return dtriangle1dy(x[0],x[1],I[0],I[1]);
  }

  //************************************************************

  double tri2::eval_phi(const int *I,const double* x) const {
    return triangle2(x[0],x[1],I[0],I[1]);
  }

  double tri2::eval_dphidx(const int *I,const double* x) const {
    return dtriangle2dx(x[0],x[1],I[0],I[1]);
  }

  double tri2::eval_dphidy(const int *I,const double* x) const {
    return dtriangle2dy(x[0],x[1],I[0],I[1]);
  }

  double tri2::eval_d2phidx2(const int *I,const double* x) const {
    return d2triangle2dx2(x[0],x[1],I[0],I[1]);
  }

  double tri2::eval_d2phidy2(const int *I,const double* x) const {
    return d2triangle2dy2(x[0],x[1],I[0],I[1]);
  }

  double tri2::eval_d2phidxdy(const int *I,const double* x) const {
    return d2triangle2dxdy(x[0],x[1],I[0],I[1]);
  }
  
} //end namespace femus

