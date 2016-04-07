/*=========================================================================

  Program: FEMUS
  Module: Quadrilateral
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
#include "cmath"


namespace femus {
  
  //************************************************************
  
  const double quad_lag::X[25][2]= { 
    {-1,-1},{1,-1},{1, 1},{-1, 1},
    { 0,-1},{1, 0},{0, 1},{-1, 0},{0, 0},
    {-0.5,-1},{0,-0.5},{-0.5,0},{-1,-0.5},
    { 0.5,-1},{1,-0.5},{ 0.5,0},
    { 1, 0.5},{0.5, 1},{0, 0.5},
    {-0.5, 1},{-1,0.5},
    {-0.5,-0.5},{0.5, -0.5},{0.5, 0.5},{-0.5, 0.5}
  };


  const int quad_lag::IND[9][2]= { 
    {0, 0},{2, 0},{2, 2},{0, 2},
    {1, 0},{2, 1},{1, 2},{0, 1},
    {1, 1}
  };

  const int quad_lag::KVERT_IND[25][2]= {
    {0,0},{1,1},{2,2},{3,3},
    {0,1},{1,2},{2,3},{3,0},{0,2},
    {0,4},{0,5},{0,6},{0,7},
    {1,4},{1,5},{1,6},
    {2,5},{2,6},{2,7},
    {3,6},{3,7},
    {0,8},{1,8},{2,8},{3,8}
  };
  
  //************************************************************
  
  const double quad_const::X[12][2]={ 
    {-0.5,-0.5},{0.5, -0.5},{0.5, 0.5},{-0.5, 0.5},
    {-0.5,-0.5},{0.5, -0.5},{0.5, 0.5},{-0.5, 0.5},
    {-0.5,-0.5},{0.5, -0.5},{0.5, 0.5},{-0.5, 0.5}
  };
  
  const int quad_const::IND[3][2]= {{0, 0},{1, 0},{0, 1}};

  const int quad_const::KVERT_IND[12][2]={ 
    {0,0},{1,0},{2,0},{3,0},
    {0,1},{1,1},{2,1},{3,1},
    {0,2},{1,2},{2,2},{3,2}
  };
  
  //************************************************************
  
  double quadpwl::eval_phi(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1]) + 
	    x[0]*eval_dphidx(I,x) + 
	    x[1]*eval_dphidy(I,x);
  }

  double quadpwl::eval_dphidx(const int *I,const double* x) const {
    return I[0];
  }

  double quadpwl::eval_dphidy(const int *I,const double* x) const {
    return I[1];
  }

  //************************************************************

  double quadLinear::eval_phi(const int *I,const double* x) const {
    return lagLinear(x[0],I[0])*lagLinear(x[1],I[1]);
  }

  double quadLinear::eval_dphidx(const int *I,const double* x) const {
    return dlagLinear(x[0],I[0])*lagLinear(x[1],I[1]);
  }

  double quadLinear::eval_dphidy(const int *I,const double* x) const {
    return lagLinear(x[0],I[0])*dlagLinear(x[1],I[1]);
  }

  double quadLinear::eval_d2phidxdy(const int *I,const double* x) const {
    return dlagLinear(x[0],I[0])*dlagLinear(x[1],I[1]);
  }

  //************************************************************

  double quadBiquadratic::eval_phi(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1]);
  }

  double quadBiquadratic::eval_dphidx(const int *I,const double* x) const {
    return dlagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1]);
  }

  double quadBiquadratic::eval_dphidy(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*dlagBiquadratic(x[1],I[1]);
  }

  double quadBiquadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2lagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1]);
  }

  double quadBiquadratic::eval_d2phidy2(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*d2lagBiquadratic(x[1],I[1]);
  }

  double quadBiquadratic::eval_d2phidxdy(const int *I,const double* x) const {
    return dlagBiquadratic(x[0],I[0])*dlagBiquadratic(x[1],I[1]);
  }

  //************************************************************

  double quadQuadratic::eval_phi(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)? 
      th2(x[0],I[0])*th2(x[1],I[1])
      : 
      (-1.+ix*x[0]+jx*x[1])*th2(x[0],I[0])*th2(x[1],I[1]);
  }

  double quadQuadratic::eval_dphidx(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      dth2(x[0],I[0])*th2(x[1],I[1])
      :
      th2(x[1],I[1])*(ix*th2(x[0],I[0]) + (-1.+ix*x[0]+jx*x[1])*dth2(x[0],I[0])) ;
  }

  double quadQuadratic::eval_dphidy(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      th2(x[0],I[0])*dth2(x[1],I[1])
      :
      th2(x[0],I[0])*(jx*th2(x[1],I[1]) + (-1.+ix*x[0]+jx*x[1])*dth2(x[1],I[1]));
  }

  double quadQuadratic::eval_d2phidx2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      d2th2(x[0],I[0])*th2(x[1],I[1])
      :
      th2(x[1],I[1])*( 2.*ix*dth2(x[0],I[0]) + (-1.+ix*x[0]+jx*x[1])*d2th2(x[0],I[0]));
  }

  double quadQuadratic::eval_d2phidy2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      th2(x[0],I[0])*d2th2(x[1],I[1])
      :
      th2(x[0],I[0])*( 2.*jx*dth2(x[1],I[1]) + (-1.+ix*x[0]+jx*x[1])*d2th2(x[1],I[1]));
  }

  double quadQuadratic::eval_d2phidxdy(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      dth2(x[0],I[0])*dth2(x[1],I[1]) 
      :
      ix*th2(x[0],I[0])*dth2(x[1],I[1])+ 
      jx*th2(x[1],I[1])*dth2(x[0],I[0])+
      (-1.+ix*x[0]+jx*x[1])*dth2(x[0],I[0])*dth2(x[1],I[1]);
  }

} //end namespace femus


