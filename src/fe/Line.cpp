/*=========================================================================

  Program: FEMUS
  Module: Line
  Authors: Eugenio Aulisa and Sara Calandrini
 
  Copyright (c) FEMTTU
  All rights reserved. 

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/



#include "Basis.hpp"


namespace femus {
  
  // line const vectors
  const double line_lag::Xc[3][1]= {{-1},{1},{0}/*,{-0.5},{0.5}*/};

  const int line_lag::IND[3][1]= {{0},{2},{1}};

  const int line_lag::KVERT_IND[5][2]= {{0,0},{1,1},{0,1},{0,2},{1,2}};
  
  const unsigned line_lag::fine2CoarseVertexMapping[2][2]={
  {0,2},
  {2,1} };
  
  //************************************************************
  
  const double line_const::X[4][1]={ 
    {-0.5},{0.5},
    {-0.5},{0.5}
  };
  
  const int line_const::IND[2][1]= {{0},{1}};

  const int line_const::KVERT_IND[4][2]={ 
    {0,0},{1,0},
    {0,1},{1,1}
  };
  
  //************************************************************
  
  double linepwLinear::eval_phi(const int *I,const double* x) const {
    return (1.-I[0]) + x[0]*eval_dphidx(I,x);
  }

  double linepwLinear::eval_dphidx(const int *I,const double* x) const {
    return I[0];
  }
 
  //************************************************************

  double LineLinear::eval_phi(const int *I,const double* x) const {
    return lagLinear(x[0],I[0]);
  }

  double LineLinear::eval_dphidx(const int *I,const double* x) const {
    return dlagLinear(x[0],I[0]);
  }

  //************************************************************

  double LineBiquadratic::eval_phi(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0]);
  }

  double LineBiquadratic::eval_dphidx(const int *I,const double* x) const {
    return dlagBiquadratic(x[0],I[0]);
  }

  double LineBiquadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2lagBiquadratic(x[0],I[0]);
  }

} //end namespace femus

