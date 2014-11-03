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

double quad0::eval_phi(const int *I,const double* x) const {
  return 1;
}



//************************************************************
double quadpwl::eval_phi(const int *I,const double* x) const {
  return 1.*I[0]+x[0]*(1.-I[0])*I[1]+x[1]*(1.-I[0])*(1.-I[1]);
}

double quadpwl::eval_dphidx(const int *I,const double* x) const {
  return (1.-I[0])*I[1];
}

double quadpwl::eval_dphidy(const int *I,const double* x) const {
  return (1.-I[0])*(1.-I[1]);
}

//************************************************************


//************************************************************

double quad1::eval_phi(const int *I,const double* x) const {
  return lag1(x[0],I[0])*lag1(x[1],I[1]);
}

double quad1::eval_dphidx(const int *I,const double* x) const {
  return dlag1(x[0],I[0])*lag1(x[1],I[1]);
}

double quad1::eval_dphidy(const int *I,const double* x) const {
  return lag1(x[0],I[0])*dlag1(x[1],I[1]);
}

double quad1::eval_d2phidxdy(const int *I,const double* x) const {
  return dlag1(x[0],I[0])*dlag1(x[1],I[1]);
}

double quad1::lag1(const double& x, const int& i) const {
  return (!i)*0.5*(1.-x)+!(i-2)*0.5*(1.+x);
}
double quad1::dlag1(const double& x, const int& i) const {
  return (!i)*(-0.5)+!(i-2)*0.5;
}

//************************************************************

double quad2::eval_phi(const int *I,const double* x) const {
  return lag2(x[0],I[0])*lag2(x[1],I[1]);
}

double quad2::eval_dphidx(const int *I,const double* x) const {
  return dlag2(x[0],I[0])*lag2(x[1],I[1]);
}

double quad2::eval_dphidy(const int *I,const double* x) const {
  return lag2(x[0],I[0])*dlag2(x[1],I[1]);
}

double quad2::eval_d2phidx2(const int *I,const double* x) const {
  return d2lag2(x[0],I[0])*lag2(x[1],I[1]);
}

double quad2::eval_d2phidy2(const int *I,const double* x) const {
  return lag2(x[0],I[0])*d2lag2(x[1],I[1]);
}

double quad2::eval_d2phidxdy(const int *I,const double* x) const {
  return dlag2(x[0],I[0])*dlag2(x[1],I[1]);
}

double quad2::lag2(const double& x, const int& i) const {
  return !i*0.5*x*(x-1.) + !(i-1)*(1.-x)*(1.+x) + !(i-2)*0.5*x*(1.+x);
}

double quad2::dlag2(const double& x, const int& i) const {
  return !i*(x-0.5) + !(i-1)*(-2.*x) + !(i-2)*(x+0.5);
}

double quad2::d2lag2(const double& x, const int& i) const {
  return !i + !(i-1)*(-2.) + !(i-2);
}


//************************************************************

double quadth::eval_phi(const int *I,const double* x) const {
  const double ix=I[0]-1.,jx=I[1]-1.;
  return (fabs(ix*jx)==0)? 
	  th2(x[0],I[0])*th2(x[1],I[1])
	  : 
	  (-1.+ix*x[0]+jx*x[1])*th2(x[0],I[0])*th2(x[1],I[1]);
}

double quadth::eval_dphidx(const int *I,const double* x) const {
  const double ix=I[0]-1.,jx=I[1]-1.;
  return (fabs(ix*jx)==0)?
          dth2(x[0],I[0])*th2(x[1],I[1])
	  :
          th2(x[1],I[1])*(ix*th2(x[0],I[0]) + (-1.+ix*x[0]+jx*x[1])*dth2(x[0],I[0])) ;
}

double quadth::eval_dphidy(const int *I,const double* x) const {
  const double ix=I[0]-1.,jx=I[1]-1.;
  return (fabs(ix*jx)==0)?
          th2(x[0],I[0])*dth2(x[1],I[1])
	  :
          th2(x[0],I[0])*(jx*th2(x[1],I[1]) + (-1.+ix*x[0]+jx*x[1])*dth2(x[1],I[1]));
}

double quadth::eval_d2phidx2(const int *I,const double* x) const {
  const double ix=I[0]-1.,jx=I[1]-1.;
  return (fabs(ix*jx)==0)?
          d2th2(x[0],I[0])*th2(x[1],I[1])
	  :
          th2(x[1],I[1])*( 2.*ix*dth2(x[0],I[0]) + (-1.+ix*x[0]+jx*x[1])*d2th2(x[0],I[0]));
}

double quadth::eval_d2phidy2(const int *I,const double* x) const {
  const double ix=I[0]-1.,jx=I[1]-1.;
  return (fabs(ix*jx)==0)?
          th2(x[0],I[0])*d2th2(x[1],I[1])
	  :
          th2(x[0],I[0])*( 2.*jx*dth2(x[1],I[1]) + (-1.+ix*x[0]+jx*x[1])*d2th2(x[1],I[1]));
}

double quadth::eval_d2phidxdy(const int *I,const double* x) const {
  const double ix=I[0]-1.,jx=I[1]-1.;
  return (fabs(ix*jx)==0)?
          dth2(x[0],I[0])*dth2(x[1],I[1]) 
	  :
          ix*th2(x[0],I[0])*dth2(x[1],I[1])+ 
	  jx*th2(x[1],I[1])*dth2(x[0],I[0])+
	  (-1.+ix*x[0]+jx*x[1])*dth2(x[0],I[0])*dth2(x[1],I[1]);
}

double quadth::th2(const double& x, const int& i) const {
  return !i*(0.5)*(1.-x) + !(i-1)*(1.-x)*(1.+x) + !(i-2)*(0.5)*(1.+x);
}

double quadth::dth2(const double& x, const int& i) const {
  return (!i)*(-0.5) + !(i-1)*(-2.*x) + !(i-2)*(0.5);
}

double quadth::d2th2(const double& x, const int& i) const {
  return  !(i-1)*(-2.);
}

} //end namespace femus


