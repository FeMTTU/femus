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
  /*const double tri_lag::X[15][2]= {
    {0, 0},      {1, 0},      {0, 1},
    {0.5, 0},    {0.5, 0.5},  {0, 0.5},
    {0.25,0},    {0.25,0.25}, {0,0.25},
    {0.75,0},    {0.75,0.25}, {0.5,0.25},  
    {0.25,0.5},  {0.25,0.75}, {0,0.75} 
  };*/
  
   const double tri_lag::X[19][2]= {
    {0, 0},         {1, 0},         {0, 1},
    {0.5, 0},       {0.5, 0.5},     {0, 0.5},   
    {0.25,0},       {0.25,0.25},    {0,0.25},   
    {0.75,0},       {0.75,0.25},    {0.5,0.25}, 
    {0.25,0.5},     {0.25,0.75},    {0,0.75},   
    {1./3., 1./3.}, {1./6., 1./6.}, {2./3., 1./6.}, {1./6, 2./3.}
  };

  /*const int tri_lag::IND[6][2]= {
    {0, 0},{2, 0},{0, 2},
    {1, 0},{1, 1},{0, 1}
  };*/
  
  const int tri_lag::IND[7][2]= {
    {0, 0},{2, 0},{0, 2},
    {1, 0},{1, 1},{0, 1},
    {7,7}
  };

  /*const int tri_lag::KVERT_IND[15][2]= {
    {0,0},{1,1},{2,2},
    {0,1},{1,2},{2,0},
    {0,3},{0,4},{0,5},
    {1,3},{1,4},{1,5},
    {2,3},{2,4},{2,5}
  };*/
  
  const int tri_lag::KVERT_IND[19][2]= {
    {0,0},{1,1},{2,2},
    {0,1},{1,2},{2,0},
    {0,3},{0,4},{0,5},
    {1,3},{1,4},{1,5},
    {2,3},{2,4},{2,5},
    {3,6},{0,6},{1,6},{2,6},
  };
  
  const double tri_const::X[12][2]={ 
    {0.,0.},{1.,0.},{0.,1.},{0.5,0.5},
    {0.,0.},{1.,0.},{0.,1.},{0.5,0.5},
    {0.,0.},{1.,0.},{0.,1.},{0.5,0.5}
  };
  
  const int tri_const::IND[3][2]= {{0, 0},{1,0},{0,1}};

  const int tri_const::KVERT_IND[12][2]={ 
    {0,0},{1,0},{2,0},{3,0},
    {0,1},{1,1},{2,1},{3,1},
    {0,2},{1,2},{2,2},{3,2}
  };
   
    
  double tripwl::eval_phi(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1]) + 
	    x[0]*eval_dphidx(I,x) + 
	    x[1]*eval_dphidy(I,x);
  }

  double tripwl::eval_dphidx(const int *I,const double* x) const {
    return I[0];
  }

  double tripwl::eval_dphidy(const int *I,const double* x) const {
    return I[1];
  }
  
   
  
  double triLinear::eval_phi(const int *I,const double* x) const {
    return triangleLinear(x[0],x[1],I[0],I[1]);
  }

  double triLinear::eval_dphidx(const int *I,const double* x) const {
    return dtriangleLineardx(x[0],x[1],I[0],I[1]);
  }

  double triLinear::eval_dphidy(const int *I,const double* x) const {
    return dtriangleLineardy(x[0],x[1],I[0],I[1]);
  }

  //************************************************************

  double triQuadratic::eval_phi(const int *I,const double* x) const {
    return triangleQuadratic(x[0],x[1],I[0],I[1]);
  }

  double triQuadratic::eval_dphidx(const int *I,const double* x) const {
    return dtriangleQuadraticdx(x[0],x[1],I[0],I[1]);
  }

  double triQuadratic::eval_dphidy(const int *I,const double* x) const {
    return dtriangleQuadraticdy(x[0],x[1],I[0],I[1]);
  }

  double triQuadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2triangleQuadraticdx2(x[0],x[1],I[0],I[1]);
  }

  double triQuadratic::eval_d2phidy2(const int *I,const double* x) const {
    return d2triangleQuadraticdy2(x[0],x[1],I[0],I[1]);
  }

  double triQuadratic::eval_d2phidxdy(const int *I,const double* x) const {
    return d2triangleQuadraticdxdy(x[0],x[1],I[0],I[1]);
  }
  
  //****************************************************************
  
  double triBiquadratic::eval_phi(const int *I,const double* x) const {
    return triangleBiquadratic(x[0],x[1],I[0],I[1]);
  }

  double triBiquadratic::eval_dphidx(const int *I,const double* x) const {
    return dtriangleBiquadraticdx(x[0],x[1],I[0],I[1]);
  }

  double triBiquadratic::eval_dphidy(const int *I,const double* x) const {
    return dtriangleBiquadraticdy(x[0],x[1],I[0],I[1]);
  }

  double triBiquadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2triangleBiquadraticdx2(x[0],x[1],I[0],I[1]);
  }

  double triBiquadratic::eval_d2phidy2(const int *I,const double* x) const {
    return d2triangleBiquadraticdy2(x[0],x[1],I[0],I[1]);
  }

  double triBiquadratic::eval_d2phidxdy(const int *I,const double* x) const {
    return d2triangleBiquadraticdxdy(x[0],x[1],I[0],I[1]);
  }
  
} //end namespace femus

