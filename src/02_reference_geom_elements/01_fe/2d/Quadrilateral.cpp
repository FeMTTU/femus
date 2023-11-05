/*=========================================================================

  Program: FEMUS
  Module: Quadrilateral
  Authors: Eugenio Aulisa and Sara Calandrini
 
  Copyright (c) FEMTTU
  All rights reserved. 

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/



#include "Basis.hpp"
#include <cmath>


namespace femus {
  
  //************************************************************
  
  const double quad_lag::Xc[9][2]= { 
    {-1,-1},{1,-1},{1, 1},{-1, 1},
    { 0,-1},{1, 0},{0, 1},{-1, 0},{0, 0}
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
  
  const unsigned quad_lag::fine2CoarseVertexMapping[4][4]= { // coarse mesh dof = f2CVM[element type][fine element][fine vertex]
    {0,4,8,7},
    {4,1,5,8},
    {8,5,2,6},
    {7,8,6,3}
  };
  
  const unsigned quad_lag::faceDofs[4][3]={
    {0, 1, 4},
    {1, 2, 5},
    {2, 3, 6},
    {3, 0, 7}
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
  
  double quadpwLinear::eval_phi(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1]) + 
	    x[0]*eval_dphidx(I,x) + 
	    x[1]*eval_dphidy(I,x);
  }

  double quadpwLinear::eval_dphidx(const int *I,const double* x) const {
    return I[0];
  }

  double quadpwLinear::eval_dphidy(const int *I,const double* x) const {
    return I[1];
  }

  //************************************************************

  double QuadLinear::eval_phi(const int *I,const double* x) const {
    return lagLinear(x[0],I[0])*lagLinear(x[1],I[1]);
  }

  double QuadLinear::eval_dphidx(const int *I,const double* x) const {
    return dlagLinear(x[0],I[0])*lagLinear(x[1],I[1]);
  }

  double QuadLinear::eval_dphidy(const int *I,const double* x) const {
    return lagLinear(x[0],I[0])*dlagLinear(x[1],I[1]);
  }

  double QuadLinear::eval_d2phidxdy(const int *I,const double* x) const {
    return dlagLinear(x[0],I[0])*dlagLinear(x[1],I[1]);
  }

  //************************************************************

  double QuadBiquadratic::eval_phi(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1]);
  }

  double QuadBiquadratic::eval_dphidx(const int *I,const double* x) const {
    return dlagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1]);
  }

  double QuadBiquadratic::eval_dphidy(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*dlagBiquadratic(x[1],I[1]);
  }

  double QuadBiquadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2lagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1]);
  }

  double QuadBiquadratic::eval_d2phidy2(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*d2lagBiquadratic(x[1],I[1]);
  }

  double QuadBiquadratic::eval_d2phidxdy(const int *I,const double* x) const {
    return dlagBiquadratic(x[0],I[0])*dlagBiquadratic(x[1],I[1]);
  }

  //************************************************************

  double QuadQuadratic::eval_phi(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)? 
      lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])
      : 
      (-1.+ix*x[0]+jx*x[1])*lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1]);
  }

  double QuadQuadratic::eval_dphidx(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      dlagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])
      :
      lagQuadratic(x[1],I[1])*(ix*lagQuadratic(x[0],I[0]) + (-1.+ix*x[0]+jx*x[1])*dlagQuadratic(x[0],I[0])) ;
  }

  double QuadQuadratic::eval_dphidy(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      lagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1])
      :
      lagQuadratic(x[0],I[0])*(jx*lagQuadratic(x[1],I[1]) + (-1.+ix*x[0]+jx*x[1])*dlagQuadratic(x[1],I[1]));
  }

  double QuadQuadratic::eval_d2phidx2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      d2lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])
      :
      lagQuadratic(x[1],I[1])*( 2.*ix*dlagQuadratic(x[0],I[0]) + (-1.+ix*x[0]+jx*x[1])*d2lagQuadratic(x[0],I[0]));
  }

  double QuadQuadratic::eval_d2phidy2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      lagQuadratic(x[0],I[0])*d2lagQuadratic(x[1],I[1])
      :
      lagQuadratic(x[0],I[0])*( 2.*jx*dlagQuadratic(x[1],I[1]) + (-1.+ix*x[0]+jx*x[1])*d2lagQuadratic(x[1],I[1]));
  }

  double QuadQuadratic::eval_d2phidxdy(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.;
    return (fabs(ix*jx)==0)?
      dlagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1]) 
      :
      ix*lagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1])+ 
      jx*lagQuadratic(x[1],I[1])*dlagQuadratic(x[0],I[0])+
      (-1.+ix*x[0]+jx*x[1])*dlagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1]);
  }

} //end namespace femus


