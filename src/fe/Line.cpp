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

