/*=========================================================================

 Program: FEMUS
 Module: Wedge
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




double wedge1::eval_phi(const int *I,const double* x) const {
  return tri1(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
}

double wedge1::eval_dphidx(const int *I,const double* x) const {
  return dtri1dx(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
}

double wedge1::eval_dphidy(const int *I,const double* x) const {
  return dtri1dy(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
}

double wedge1::eval_dphidz(const int *I,const double* x) const {
  return tri1(x[0],x[1],I[0],I[1])*dlag1(x[2],I[2]);
}

double wedge1::tri1(const double& x, const double& y, const int& i,const int& j) const {
  return (!i*!j)*(1.-x-y)+ !(i-2)*x + !(j-2)*y;
}

double wedge1::dtri1dx(const double& x, const double& y, const int& i,const int& j) const {
  return -(!i*!j) + !(i-2);
}

double wedge1::dtri1dy(const double& x, const double& y, const int& i,const int& j) const {
  return -(!i*!j) + !(j-2);
}

double wedge1::lag1(const double& x, const int& i) const {
  return (!i)*0.5*(1.-x)+!(i-2)*0.5*(1.+x);
}

double wedge1::dlag1(const double& x, const int& i) const {
  return (!i)*(-0.5)+!(i-2)*0.5;
}

//************************************************************

double wedge2::eval_phi(const int *I,const double* x) const {
  return tri2(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_dphidx(const int *I,const double* x) const {
  return dtri2dx(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_dphidy(const int *I,const double* x) const {
  return dtri2dy(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_dphidz(const int *I,const double* x) const {
  return tri2(x[0],x[1],I[0],I[1])*dlag2(x[2],I[2]);
}

double wedge2::eval_d2phidx2(const int *I,const double* x) const {
  return d2tri2dx2(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_d2phidy2(const int *I,const double* x) const {
  return d2tri2dy2(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_d2phidz2(const int *I,const double* x) const {
  return tri2(x[0],x[1],I[0],I[1])*d2lag2(x[2],I[2]);
}

double wedge2::eval_d2phidxdy(const int *I,const double* x) const {
  return d2tri2dxdy(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_d2phidydz(const int *I,const double* x) const {
  return dtri2dy(x[0],x[1],I[0],I[1])*dlag2(x[2],I[2]);
}

double wedge2::eval_d2phidzdx(const int *I,const double* x) const {
  return dtri2dx(x[0],x[1],I[0],I[1])*dlag2(x[2],I[2]);
}

double wedge2::lag2(const double& x, const int& i) const {
  return !i*0.5*x*(x-1.) + !(i-1)*(1.-x)*(1.+x) + !(i-2)*0.5*x*(1.+x);
}

double wedge2::dlag2(const double& x, const int& i) const {
  return !i*(x-0.5) + !(i-1)*(-2.*x) + !(i-2)*(x+0.5);
}

double wedge2::d2lag2(const double& x, const int& i) const {
  return !i + !(i-1)*(-2.) + !(i-2);
}

double wedge2::tri2(const double& x, const double& y, const int& i,const int& j) const {
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

double wedge2::dtri2dx(const double& x, const double& y, const int& i,const int& j) const {
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


double wedge2::dtri2dy(const double& x, const double& y, const int& i,const int& j) const {
  
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

double wedge2::d2tri2dx2(const double& x, const double& y, const int& i,const int& j) const {
  return !j*( (!i)*4. +!(i-1)*(-8.) + !(i-2)*4. );
}

double wedge2::d2tri2dy2(const double& x, const double& y, const int& i,const int& j) const {
  return !i*( (!j)*4. +!(j-1)*(-8.) + !(j-2)*4. );
}

double wedge2::d2tri2dxdy(const double& x, const double& y, const int& i,const int& j) const {
  return ( (!i)*(!j) + !(i-1)*!(j-1) )*4. + ( !(i-1)*(!j) + (!i)*!(j-1) )*(-4.);
}


//************************************************************
double  wedgeth::eval_phi(const int *I,const double* X) const {
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  
  double t=1.-(x+y);
  return
    !i     *( !j     *( !k *t*(-2.+2.*t-z)*(1.-z)*0.5   + !(k-1) *t*(1.-z*z) 	+ !(k-2) *t*(-2.+2.*t+z)*(1.+z)*0.5 	) +
	      !(j-1) *( !k *2.*y*t*(1.-z) 				        + !(k-2) *2.*y*t*(1.+z)          	) +
	      !(j-2) *( !k *y*(-2.+2.*y-z)*(1.-z)*0.5   + !(k-1) *y*(1.-z*z)      + !(k-2) *y*(-2.+2.*y+z)*(1.+z)*0.5   ) ) +
    !(i-1) *( !j     *( !k *2.*x*t*(1.-z)						+ !(k-2) *2.*x*t*(1.+z)	        )+
              !(j-1) *( !k *2.*x*y*(1.-z)						+ !(k-2) *2.*x*y*(1.+z)		) ) +
    !(i-2) *(         ( !k *x*(-2.+2.*x-z)*(1.-z)*0.5   + !(k-1) *x*(1.-z*z)	+ !(k-2) *x*(-2.+2.*x+z)*(1.+z)*0.5    	) );
  
  ////////////////////////////
//  
}

double  wedgeth::eval_dphidx(const int *I,const double* X) const {
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  
  double t=1.-(x+y);
  return
    !i     *( !j   *( !k *(1.-2.*t+0.5*z)*(1.-z)      + !(k-1) *(z*z-1.) 	+ !(k-2) *(1.-2.*t-0.5*z)*(1.+z)      )+
	    !(j-1) *( !k *(-2.)*y*(1.-z) 				        + !(k-2) *(-2.)*y*(1.+z)              ) ) +
    !(i-1) *( !j   *( !k *2.*(1.-z)*(1.-2.*x-y)					+ !(k-2) *2.*(1.+z)*(1.-2.*x-y)       )+
	    !(j-1) *( !k *2.*y*(1.-z)						+ !(k-2) *2.*y*(1.+z)		      ) ) +
    !(i-2) *(	    ( !k *(-1.+2.*x-0.5*z)*(1.-z)     + !(k-1) *(1.-z*z)	+ !(k-2) *(-1.+2.*x+0.5*z)*(1.+z)     ) );
  
  ////////////////////
//   double f=0.;
//   t=x+y;
//   
//   switch (i) {
//   case 0:
//     switch (j) {
//     case 0:
//       switch (k) {
//       case 0:
//         f=(-1.+2.*t+0.5*z)*(1.-z);
//         break;
//       case 1:
//         f=(z*z-1.);
//         break;
//       case 2:
//         f=(-1.+2.*t-0.5*z)*(1.+z);
//         break;
//       }
//       break;
//     case 1:
//       switch (k) {
//       case 0:
//         f=-2.*y*(1-z);
//         break;
//       case 2:
//         f=-2.*y*(1+z);
//         break;
//       }
//       break;
//     case 2:
//       switch (k) {
//       case 0:
//         f=0.;
//         break;
//       case 1:
//         f=0.;
//         break;
//       case 2:
//         f=0.;
//         break;
//       }
//       break;
//     }
//     break;
//   case 1:
//     switch (j) {
//     case 0:
//       switch (k) {
//       case 0:
//         f=2*(1.-z)*(1.-2.*x-y);
//         break;
//       case 2:
//         f=2*(1.+z)*(1.-2.*x-y);
//         break;
//       }
//       break;
//     case 1:
//       switch (k) {
//       case 0:
//         f=2.*y*(1.-z);
//         break;
//       case 2:
//         f=2.*y*(1.+z);
//         break;
//       }
//       break;
//     }
//     break;
//   case 2:
//     switch (k) {
//     case 0:
//       f=(-2.+4.*x-z)*(1.-z)*0.5;
//       break;
//     case 1:
//       f=(1.-z*z);
//       break;
//     case 2:
//       f=(-2.+4.*x+z)*(1+z)*0.5;
//       break;
//     }
//     break;
//   }
//   return f;
  
  
}

double  wedgeth::eval_dphidy(const int *I,const double* X) const {
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  
   double t=1.-(x+y);
  return
    !i     *( !j     *( !k *(1.-2.*t+0.5*z)*(1.-z)      + !(k-1) *(z*z-1.)	+ !(k-2) *(1.-2.*t-0.5*z)*(1.+z)     ) +
	      !(j-1) *( !k *2.*(1-z)*(1.-x-2.*y) 				+ !(k-2) *2.*(1+z)*(1.-x-2.*y)        ) +
	      !(j-2) *( !k *(-1.+2.*y-0.5*z)*(1.-z)     + !(k-1) *(1.-z*z)      + !(k-2) *(-1.+2.*y+0.5*z)*(1.+z)     ) ) +
    !(i-1) *( !j     *( !k *(-2.)*x*(1.-z)					+ !(k-2) *(-2.)*x*(1.+z)	      )+
	      !(j-1) *( !k *2.*x*(1.-z)						+ !(k-2) *2.*x*(1.+z)		      ) );
  /////////////////////////
//   double f=0.;
//   t=x+y;
//   
//   switch (i) {
//   case 0:
//     switch (j) {
//     case 0:
//       switch (k) {
//       case 0:
//         f=(-1.+2.*t+0.5*z)*(1.-z);
//         break;
//       case 1:
//         f=(z*z-1.);
//         break;
//       case 2:
//         f=(-1.+2.*t-0.5*z)*(1.+z);
//         break;
//       }
//       break;
//     case 1:
//       switch (k) {
//       case 0:
//         f=2.*(1-z)*(1.-x-2.*y);
//         break;
//       case 2:
//         f=2.*(1+z)*(1.-x-2.*y);
//         break;
//       }
//       break;
//     case 2:
//       switch (k) {
//       case 0:
//         f=(-2.+4.*y-z)*(1.-z)*0.5;
//         break;
//       case 1:
//         f=(1.-z*z);
//         break;
//       case 2:
//         f=(-2.+4.*y+z)*(1.+z)*0.5;
//         break;
//       }
//       break;
//     }
//     break;
//   case 1:
//     switch (j) {
//     case 0:
//       switch (k) {
//       case 0:
//         f=-2.*x*(1.-z);
//         break;
//       case 2:
//         f=-2.*x*(1.+z);
//         break;
//       }
//       break;
//     case 1:
//       switch (k) {
//       case 0:
//         f=2.*x*(1.-z);
//         break;
//       case 2:
//         f=2.*x*(1.+z);
//         break;
//       }
//       break;
//     }
//     break;
//   case 2:
//     switch (k) {
//     case 0:
//       f=0.;
//       break;
//     case 1:
//       f=0.;
//       break;
//     case 2:
//       f=0.;
//       break;
//     }
//     break;
//   }
//   return f;
}

double  wedgeth::eval_dphidz(const int *I,const double* X) const {
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  double t=1.-(x+y);
  return
    !i     *( !j   *( !k *t*(0.5-t+z)	   + !(k-1) *(-2.)*t*z  	+ !(k-2) *t*(-0.5+t+z) 	  ) +
	    !(j-1) *( !k *(-2.)*y*t				        + !(k-2) *2.*y*t	  ) +
	    !(j-2) *( !k *y*(0.5-y+z)	   + !(k-1) *(-2.)*y*z		+ !(k-2) *y*(-0.5+y+z)	  ) ) +
    !(i-1) *( !j   *( !k *(-2.)*x*t					+ !(k-2) *2.*x*t    	  ) +
	    !(j-1) *( !k *(-2.)*x*y					+ !(k-2) *2.*x*y	  ) ) +
    !(i-2) *(	    ( !k *x*(0.5-x+z)	   + !(k-1) *(-2.)*x*z		+ !(k-2) *x*(-0.5+x+z)    ) );
 
  
  /////////////////
  
//   double f=0.; 
//   switch (i) {
//   case 0:
//     switch (j) {
//     case 0:
//       switch (k) {
//       case 0:
//         f=t*(0.5-t+z);
//         break;
//       case 1:
//         f=-2.*t*z;
//         break;
//       case 2:
//         f=t*(-0.5+t+z);
//         break;
//       }
//       break;
//     case 1:
//       switch (k) {
//       case 0:
//         f=-2.*y*t;
//         break;
//       case 2:
//         f=2.*y*t;
//         break;
//       }
//       break;
//     case 2:
//       switch (k) {
//       case 0:
//         f=y*(0.5-y+z);
//         break;
//       case 1:
//         f=-2.*y*z;
//         break;
//       case 2:
//         f=y*(-0.5+y+z);
//         break;
//       }
//       break;
//     }
//     break;
//   case 1:
//     switch (j) {
//     case 0:
//       switch (k) {
//       case 0:
//         f=-2.*x*t;
//         break;
//       case 2:
//         f=2.*x*t;
//         break;
//       }
//       break;
//     case 1:
//       switch (k) {
//       case 0:
//         f=-2.*x*y;
//         break;
//       case 2:
//         f=2.*x*y;
//         break;
//       }
//       break;
//     }
//     break;
//   case 2:
//     switch (k) {
//     case 0:
//       f=x*(0.5-x+z);
//       break;
//     case 1:
//       f=-2.*x*z;
//       break;
//     case 2:
//       f=x*(-0.5+x+z);
//       break;
//     }
//     break;
//   }
//   return f;
  
}

double wedgeth::eval_d2phidx2(const int *I,const double* X) const {
  
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  
  double t=1.-(x+y);
  return
    !i     *( !j  *( !k *(2.)*(1.-z)  + !(k-2) *(2.)*(1.+z)   ) ) +
    !(i-1) *( !j  *( !k *(-4.)*(1.-z) + !(k-2) *(-4.)*(1.+z)  ) ) +
    !(i-2) *(	   ( !k *(2.)*(1.-z)  + !(k-2) *(2.)*(1.+z)   ) );
}

double wedgeth::eval_d2phidy2(const int *I,const double* X) const {
  
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  
  double t=1.-(x+y);
  return
    !i *( !j     *( !k *(2.) *(1.-z)	+ !(k-2) *(2.) *(1.+z)   ) +
	  !(j-1) *( !k *(-4.)*(1.-z)	+ !(k-2) *(-4.)*(1.+z)   ) +
	  !(j-2) *( !k *(2.) *(1.-z)	+ !(k-2) *(2.) *(1.+z)   ) );
  
}

double wedgeth::eval_d2phidz2(const int *I,const double* X) const {
  
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  
  double t=1.-(x+y);
  return
    !i     *( !j     *( !k *t  + !(k-1) *(-2.)*t  + !(k-2) *t  ) +
	      !(j-2) *( !k *y  + !(k-1) *(-2.)*y  + !(k-2) *y  ) ) +
    !(i-2) *(	      ( !k *x  + !(k-1) *(-2.)*x  + !(k-2) *x  ) );
}

double  wedgeth::eval_d2phidxdy(const int *I,const double* X) const {
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  
  double t=1.-(x+y);
  return
    !i     *( !j   *( !k *(2.)*(1.-z)      + !(k-2) *(2.)*(1.+z)      )+
	    !(j-1) *( !k *(-2.)*(1.-z) 	   + !(k-2) *(-2.)*(1.+z)     ) ) +
    !(i-1) *( !j   *( !k *(-2.)*(1.-z)	   + !(k-2) *(-2.)*(1.+z)     )+
	    !(j-1) *( !k *(2.)*(1.-z)	   + !(k-2) *(2.)*(1.+z)      ) );
  
}

double  wedgeth::eval_d2phidydz(const int *I,const double* X) const {
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  
   double t=1.-(x+y);
  return
    !i     *( !j     *( !k *(-0.5+2.*t-z)	+ !(k-1) *(2.)*z	+ !(k-2) *(0.5-2.*t-z) 	 ) +
	      !(j-1) *( !k *(-2.)*(1.-x-2.*y) 			        + !(k-2) *2.*(1.-x-2.*y) ) +
	      !(j-2) *( !k *(0.5-2.*y+z)	+ !(k-1) *(-2.)*z       + !(k-2) *(-0.5+2.*y+z)   ) ) +
    !(i-1) *( !j     *( !k *(2.)*x					+ !(k-2) *(-2.)*x	 )+
	      !(j-1) *( !k *(-2.)*x					+ !(k-2) *2.*x		 ) );
}

double  wedgeth::eval_d2phidzdx(const int *I,const double* X) const {
  const double x=X[0];   const double y=X[1];   const double z=X[2];
  const int i=I[0];      const int j=I[1];      const int k=I[2];
  double t=1.-(x+y);
  return
    !i     *( !j   *( !k *(-0.5+2*t-z)	   + !(k-1) *(2.)*z  	 + !(k-2) *(0.5-2.*t-z)   ) +
	    !(j-1) *( !k *(2.)*y				 + !(k-2) *(-2.)*y	  ) ) +
    !(i-1) *( !j   *( !k *(-2.)*(1-2.*x-y)			 + !(k-2) *2.*(1-2.*x-y)  ) +
	    !(j-1) *( !k *(-2.)*y				 + !(k-2) *2.*y	  	  ) ) +
    !(i-2) *(	    ( !k *(0.5-2.*x+z)	   + !(k-1) *(-2.)*z	 + !(k-2) *(-0.5+2.*x+z)  ) );
}


} //end namespace femus


