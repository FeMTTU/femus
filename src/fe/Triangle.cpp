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


double tri0::eval_phi(const int *I,const double* x) const {
  return 1.;
}

//************************************************************

double tri1::eval_phi(const int *I,const double* x) const {
  return tri1a(x[0],x[1],I[0],I[1]);
}

double tri1::eval_dphidx(const int *I,const double* x) const {
  return dtri1dx(x[0],x[1],I[0],I[1]);
}

double tri1::eval_dphidy(const int *I,const double* x) const {
  return dtri1dy(x[0],x[1],I[0],I[1]);
}

double tri1::tri1a(const double& x, const double& y, const int& i,const int& j) const {
  return (!i*!j)*(1.-x-y)+ !(i-2)*x + !(j-2)*y;
}

double tri1::dtri1dx(const double& x, const double& y, const int& i,const int& j) const {
  return -(!i*!j) + !(i-2);
}

double tri1::dtri1dy(const double& x, const double& y, const int& i,const int& j) const {
  return -(!i*!j) + !(j-2);
}

//************************************************************

double tri2::eval_phi(const int *I,const double* x) const {
  return tri2a(x[0],x[1],I[0],I[1]);
}

double tri2::eval_dphidx(const int *I,const double* x) const {
  return dtri2dx(x[0],x[1],I[0],I[1]);
}

double tri2::eval_dphidy(const int *I,const double* x) const {
  return dtri2dy(x[0],x[1],I[0],I[1]);
}

double tri2::eval_d2phidx2(const int *I,const double* x) const {
  return d2tri2dx2(x[0],x[1],I[0],I[1]);
}

double tri2::eval_d2phidy2(const int *I,const double* x) const {
  return d2tri2dy2(x[0],x[1],I[0],I[1]);
}

double tri2::eval_d2phidxdy(const int *I,const double* x) const {
  return d2tri2dxdy(x[0],x[1],I[0],I[1]);
}

double tri2::tri2a(const double& x, const double& y, const int& i,const int& j) const {
  
  return !i     * (!j* (1.-x-y)*(1.-2.*x-2.*y) + !(j-1)* 4.*y*(1.-x-y) + !(j-2)*(-y+2.*y*y)) +
	 !(i-1) * (!j* 4.*x*(1.-x-y) + !(j-1) * 4.*x*y) +
	 !(i-2) * (!j*(-x+2.*x*x));   
//   double f=0.;
//   switch (i) {
//   case(0):
//     switch (j) {
//     case(0):
//       f=(1.-x-y)*(1.-2.*x-2.*y);
//       break;
//     case(1):
//       f=4.*y*(1.-x-y);
//       break;
//     case(2):
//       f=(-y+2.*y*y);
//       break;
//     }
//     break;
//   case(1):
//     switch (j) {
//     case(0):
//       f=4.*x*(1.-x-y);
//       break;
//     case(1):
//       f=4.*x*y;
//       break;
//     }
//     break;
//   case(2):
//     f=(-x+2.*x*x);
//     break;
//   }
//   return f;
}

double tri2::dtri2dx(const double& x, const double& y, const int& i,const int& j) const {
  
  return !i     * (!j* (-3.+4.*x+4.*y) + !(j-1)*y*(-4.) ) +
	 !(i-1) * (!j* 4.*(1.-2.*x-y)  + !(j-1)*y*(4.)) +
	 !(i-2) * (!j*(-1 + 4.*x));
    
//   double f=0.;
//   switch (i) {
//   case(0):
//     switch (j) {
//     case(0):
//       f=(-3.+4.*x+4.*y);
//       break;
//     case(1):
//       f=-4.*y;
//       break;
//     case(2):
//       f=0;
//       break;
//     }
//     break;
//   case(1):
//     switch (j) {
//     case(0):
//       f=4.*(1-2.*x-y);
//       break;
//     case(1):
//       f=4.*y;
//       break;
//     }
//     break;
//   case(2):
//     f=-1.+4.*x;
//     break;
//   }
//   return f;
}


double tri2::dtri2dy(const double& x, const double& y, const int& i,const int& j) const {
  
  return !j     * (!i* (-3.+4.*y+4.*x) + !(i-1)*x*(-4.) ) +
	 !(j-1) * (!i* 4.*(1.-2.*y-x)  + !(i-1)*x*(4.)) +
	 !(j-2) * (!i*(-1 + 4.*y));	 
	  
//   double f=0.;
//   switch (i) {
//   case(0):
//     switch (j) {
//     case(0):
//       f=(-3.+4.*x+4.*y);
//       break;
//     case(1):
//       f=4.*(1-x-2.*y);
//       break;
//     case(2):
//       f=-1.+4.*y;
//       break;
//     }
//     break;
//   case(1):
//     switch (j) {
//     case(0):
//       f=-4.*x;
//       break;
//     case(1):
//       f=4.*x;
//       break;
//     }
//     break;
//   case(2):
//     f=0.;
//     break;
//   }
//   return f;
}

double tri2::d2tri2dx2(const double& x, const double& y, const int& i,const int& j) const {
  return !j*( (!i)*4. +!(i-1)*(-8.) + !(i-2)*4. );
}

double tri2::d2tri2dy2(const double& x, const double& y, const int& i,const int& j) const {
  return !i*( (!j)*4. +!(j-1)*(-8.) + !(j-2)*4. );
}

double tri2::d2tri2dxdy(const double& x, const double& y, const int& i,const int& j) const {
  return (!i*!j + !(i-1)*!(j-1))*4. + ((!i-1)*!j + !i*!(j-1))*(-4.);
}

} //end namespace femus

