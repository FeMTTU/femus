/*=========================================================================

  Program: FEMUS
  Module: Line
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
  
  // line const vectors
  const double line_lag::X[5][1]= {{-1},{1},{0},{-0.5},{0.5}};


  const int line_lag::IND[3][1]= {{0},{2},{1}};


  const int line_lag::KVERT_IND[5][2]= {{0,0},{1,1},{0,1},{0,2},{1,2}};
  
  
  
  const double line_const::X[4][1]={ 
    {-0.5},{0.5},
    {-0.5},{0.5}
  };
  
  const int line_const::IND[2][1]= {{1},{0}};

  const int line_const::KVERT_IND[4][2]={ 
    {0,0},{1,0},
    {0,1},{1,1}
  };
  
  
  double linepwl::eval_phi(const int *I,const double* x) const {
    return  I[0]+x[0]*(1.-I[0]);
  }

  double linepwl::eval_dphidx(const int *I,const double* x) const {
    return (1.-I[0]);
  }

  
  
  

  double line1::eval_phi(const int *I,const double* x) const {
    return lag1(x[0],I[0]);
  }

  double line1::eval_dphidx(const int *I,const double* x) const {
    return dlag1(x[0],I[0]);
  }

  //************************************************************

  double line2::eval_phi(const int *I,const double* x) const {
    return lag2(x[0],I[0]);
  }

  double line2::eval_dphidx(const int *I,const double* x) const {
    return dlag2(x[0],I[0]);
  }

  double line2::eval_d2phidx2(const int *I,const double* x) const {
    return d2lag2(x[0],I[0]);
  }

} //end namespace femus

