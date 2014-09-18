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


//************************************************************

double line0::eval_phi(const int *I,const double* x) const {
  return lag0(x[0],I[0]);
}

double line0::eval_dphidx(const int *I,const double* x) const {
  return dlag0(x[0],I[0]);
}

double line0::eval_dphidy(const int *I,const double* x) const {
  return 0;
}

double line0::eval_dphidz(const int *I,const double* x) const {
  return 0.;
}

double line0::lag0(const double& x, const int& i) const {
  return 1.;
}
double line0::dlag0(const double& x, const int& i) const {
  return 0.;
}


//************************************************************

double line1::eval_phi(const int *I,const double* x) const {
  return lag1(x[0],I[0]);
}

double line1::eval_dphidx(const int *I,const double* x) const {
  return dlag1(x[0],I[0]);
}

double line1::eval_dphidy(const int *I,const double* x) const {
  return 0;
}

double line1::eval_dphidz(const int *I,const double* x) const {
  return 0.;
}

double line1::lag1(const double& x, const int& i) const {
  return (!i)*0.5*(1.-x)+!(i-2)*0.5*(1.+x);
}
double line1::dlag1(const double& x, const int& i) const {
  return (!i)*(-0.5)+!(i-2)*0.5;
}

//************************************************************

double line2::eval_phi(const int *I,const double* x) const {
  return lag2(x[0],I[0]);
}

double line2::eval_dphidx(const int *I,const double* x) const {
  return dlag2(x[0],I[0]);
}

double line2::eval_dphidy(const int *I,const double* x) const {
  return 0;
}

double line2::eval_dphidz(const int *I,const double* x) const {
  return 0.;
}

double line2::lag2(const double& x, const int& i) const {
  return !i*0.5*x*(x-1.) + !(i-1)*(1.-x)*(1.+x) + !(i-2)*0.5*x*(1.+x);
}

double line2::dlag2(const double& x, const int& i) const {
  return !i*(x-0.5) + !(i-1)*(-2.*x) + !(i-2)*(x+0.5);
}


} //end namespace femus

