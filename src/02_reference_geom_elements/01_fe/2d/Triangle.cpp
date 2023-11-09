/*=========================================================================

  Program: FEMUS
  Module: Triangle
  Authors: Eugenio Aulisa and Sara Calandrini
 
  Copyright (c) FEMTTU
  All rights reserved. 

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/



#include "Triangle.hpp"


namespace femus {

   const double tri_lag::Xc[7][2]= {
    {0, 0},         {1, 0},         {0, 1},
    {0.5, 0},       {0.5, 0.5},     {0, 0.5},  
    {1./3., 1./3.}
   };
  
  const int tri_lag::IND[7][2]= {
    {0, 0},{2, 0},{0, 2},
    {1, 0},{1, 1},{0, 1},
    {7,7}
  };
  
  const int tri_lag::KVERT_IND[19][2]= {
    {0,0},{1,1},{2,2},
    {0,1},{1,2},{2,0},
    {0,3},{0,4},{0,5},
    {1,3},{1,4},{1,5},
    {2,3},{2,4},{2,5},
    {3,6},{0,6},{1,6},{2,6},
  };
  
  const unsigned tri_lag::fine2CoarseVertexMapping[4][3]= { 
    {0,3,5},
    {3,1,4},
    {5,4,2},
    {4,5,3} 
  };
    
  const unsigned tri_lag::faceDofs[3][3] = { 
    {0, 1, 3},
    {1, 2, 4},
    {2, 0, 5}
  };
  
  //******************************************************************
  
  const double tri_const::X[12][2]={ 
    {1./6.,1./6.},{2./3.,1./6.},{1./6.,2./3.},{1./3.,1./3.},
    {1./6.,1./6.},{2./3.,1./6.},{1./6.,2./3.},{1./3.,1./3.},
    {1./6.,1./6.},{2./3.,1./6.},{1./6.,2./3.},{1./3.,1./3.}
  };
  
  const int tri_const::IND[3][2]= {{0, 0},{1,0},{0,1}};

  const int tri_const::KVERT_IND[12][2]={ 
    {0,0},{1,0},{2,0},{3,0},
    {0,1},{1,1},{2,1},{3,1},
    {0,2},{1,2},{2,2},{3,2}
  };
   
  //*****************************************************************
    
  double tripwLinear::eval_phi(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1]) + 
	    (x[0]-1./3.)*eval_dphidx(I,x) + 
	    (x[1]-1./3.)*eval_dphidy(I,x);
  }

  double tripwLinear::eval_dphidx(const int *I,const double* x) const {
    return I[0];
  }

  double tripwLinear::eval_dphidy(const int *I,const double* x) const {
    return I[1];
  }
  
   
  
  double TriLinear::eval_phi(const int *I,const double* x) const {
    return triangleLinear(x[0],x[1],I[0],I[1]);
  }

  double TriLinear::eval_dphidx(const int *I,const double* x) const {
    return dtriangleLineardx(x[0],x[1],I[0],I[1]);
  }

  double TriLinear::eval_dphidy(const int *I,const double* x) const {
    return dtriangleLineardy(x[0],x[1],I[0],I[1]);
  }

  //************************************************************

  double TriQuadratic::eval_phi(const int *I,const double* x) const {
    return triangleQuadratic(x[0],x[1],I[0],I[1]);
  }

  double TriQuadratic::eval_dphidx(const int *I,const double* x) const {
    return dtriangleQuadraticdx(x[0],x[1],I[0],I[1]);
  }

  double TriQuadratic::eval_dphidy(const int *I,const double* x) const {
    return dtriangleQuadraticdy(x[0],x[1],I[0],I[1]);
  }

  double TriQuadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2triangleQuadraticdx2(x[0],x[1],I[0],I[1]);
  }

  double TriQuadratic::eval_d2phidy2(const int *I,const double* x) const {
    return d2triangleQuadraticdy2(x[0],x[1],I[0],I[1]);
  }

  double TriQuadratic::eval_d2phidxdy(const int *I,const double* x) const {
    return d2triangleQuadraticdxdy(x[0],x[1],I[0],I[1]);
  }
  
  //****************************************************************
  
  double TriBiquadratic::eval_phi(const int *I,const double* x) const {
    return triangleBiquadratic(x[0],x[1],I[0],I[1]);
  }

  double TriBiquadratic::eval_dphidx(const int *I,const double* x) const {
    return dtriangleBiquadraticdx(x[0],x[1],I[0],I[1]);
  }

  double TriBiquadratic::eval_dphidy(const int *I,const double* x) const {
    return dtriangleBiquadraticdy(x[0],x[1],I[0],I[1]);
  }

  double TriBiquadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2triangleBiquadraticdx2(x[0],x[1],I[0],I[1]);
  }

  double TriBiquadratic::eval_d2phidy2(const int *I,const double* x) const {
    return d2triangleBiquadraticdy2(x[0],x[1],I[0],I[1]);
  }

  double TriBiquadratic::eval_d2phidxdy(const int *I,const double* x) const {
    return d2triangleBiquadraticdxdy(x[0],x[1],I[0],I[1]);
  }
  
} //end namespace femus

