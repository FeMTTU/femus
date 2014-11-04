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
    return triangle1(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
  }

  double wedge1::eval_dphidx(const int *I,const double* x) const {
    return dtriangle1dx(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
  }

  double wedge1::eval_dphidy(const int *I,const double* x) const {
    return dtriangle1dy(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
  }

  double wedge1::eval_dphidz(const int *I,const double* x) const {
    return triangle1(x[0],x[1],I[0],I[1])*dlag1(x[2],I[2]);
  }

  //************************************************************

  double wedge2::eval_phi(const int *I,const double* x) const {
    return triangle2(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
  }

  double wedge2::eval_dphidx(const int *I,const double* x) const {
    return dtriangle2dx(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
  }

  double wedge2::eval_dphidy(const int *I,const double* x) const {
    return dtriangle2dy(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
  }

  double wedge2::eval_dphidz(const int *I,const double* x) const {
    return triangle2(x[0],x[1],I[0],I[1])*dlag2(x[2],I[2]);
  }

  double wedge2::eval_d2phidx2(const int *I,const double* x) const {
    return d2triangle2dx2(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
  }

  double wedge2::eval_d2phidy2(const int *I,const double* x) const {
    return d2triangle2dy2(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
  }

  double wedge2::eval_d2phidz2(const int *I,const double* x) const {
    return triangle2(x[0],x[1],I[0],I[1])*d2lag2(x[2],I[2]);
  }

  double wedge2::eval_d2phidxdy(const int *I,const double* x) const {
    return d2triangle2dxdy(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
  }

  double wedge2::eval_d2phidydz(const int *I,const double* x) const {
    return dtriangle2dy(x[0],x[1],I[0],I[1])*dlag2(x[2],I[2]);
  }

  double wedge2::eval_d2phidzdx(const int *I,const double* x) const {
    return dtriangle2dx(x[0],x[1],I[0],I[1])*dlag2(x[2],I[2]);
  }

  //************************************************************
  
  double  wedgeth::eval_phi(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j     *( !k *t*(-2.+2.*t-z)*(1.-z)*0.5   + !(k-1) *t*(1.-z*z) 	+ !(k-2) *t*(-2.+2.*t+z)*(1.+z)*0.5 	) +
		!(j-1) *( !k *2.*y*t*(1.-z) 				        + !(k-2) *2.*y*t*(1.+z)          	) +
		!(j-2) *( !k *y*(-2.+2.*y-z)*(1.-z)*0.5   + !(k-1) *y*(1.-z*z)  + !(k-2) *y*(-2.+2.*y+z)*(1.+z)*0.5   ) ) +
      !(i-1) *( !j     *( !k *2.*x*t*(1.-z)					+ !(k-2) *2.*x*t*(1.+z)	        )+
		!(j-1) *( !k *2.*x*y*(1.-z)					+ !(k-2) *2.*x*y*(1.+z)		) ) +
      !(i-2) *(         ( !k *x*(-2.+2.*x-z)*(1.-z)*0.5   + !(k-1) *x*(1.-z*z)	+ !(k-2) *x*(-2.+2.*x+z)*(1.+z)*0.5    	) );
    
  }

  double  wedgeth::eval_dphidx(const int *I,const double* X) const {
  
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j   *( !k *(1.-2.*t+0.5*z)*(1.-z)      + !(k-1) *(z*z-1.) 	+ !(k-2) *(1.-2.*t-0.5*z)*(1.+z)      )+
		!(j-1) *( !k *(-2.)*y*(1.-z) 				        + !(k-2) *(-2.)*y*(1.+z)              ) ) +
      !(i-1) *( !j   *( !k *2.*(1.-z)*(1.-2.*x-y)				+ !(k-2) *2.*(1.+z)*(1.-2.*x-y)       )+
		!(j-1) *( !k *2.*y*(1.-z)					+ !(k-2) *2.*y*(1.+z)		      ) ) +
      !(i-2) *(	    ( !k *(-1.+2.*x-0.5*z)*(1.-z)     + !(k-1) *(1.-z*z)	+ !(k-2) *(-1.+2.*x+0.5*z)*(1.+z)     ) );
  }

  double  wedgeth::eval_dphidy(const int *I,const double* X) const {
  
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j     *( !k *(1.-2.*t+0.5*z)*(1.-z)      + !(k-1) *(z*z-1.)	+ !(k-2) *(1.-2.*t-0.5*z)*(1.+z)      ) +
		!(j-1) *( !k *2.*(1-z)*(1.-x-2.*y) 				+ !(k-2) *2.*(1+z)*(1.-x-2.*y)        ) +
		!(j-2) *( !k *(-1.+2.*y-0.5*z)*(1.-z)     + !(k-1) *(1.-z*z)    + !(k-2) *(-1.+2.*y+0.5*z)*(1.+z)     ) ) +
      !(i-1) *( !j     *( !k *(-2.)*x*(1.-z)					+ !(k-2) *(-2.)*x*(1.+z)	      )+
		!(j-1) *( !k *2.*x*(1.-z)					+ !(k-2) *2.*x*(1.+z)		      ) );
  
  }

  double  wedgeth::eval_dphidz(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j     *( !k *t*(0.5-t+z) + !(k-1) *(-2.)*t*z  + !(k-2) *t*(-0.5+t+z) 	) +
		!(j-1) *( !k *(-2.)*y*t				+ !(k-2) *2.*y*t	) +
		!(j-2) *( !k *y*(0.5-y+z) + !(k-1) *(-2.)*y*z	+ !(k-2) *y*(-0.5+y+z)	) ) +
      !(i-1) *( !j     *( !k *(-2.)*x*t				+ !(k-2) *2.*x*t    	) +
		!(j-1) *( !k *(-2.)*x*y				+ !(k-2) *2.*x*y	) ) +
      !(i-2) *(	        ( !k *x*(0.5-x+z) + !(k-1) *(-2.)*x*z	+ !(k-2) *x*(-0.5+x+z)  ) );
 
  }

  double wedgeth::eval_d2phidx2(const int *I,const double* X) const {
  
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j  *( !k *(2.)*(1.-z)  + !(k-2) *(2.)*(1.+z)   ) ) +
      !(i-1) *( !j  *( !k *(-4.)*(1.-z) + !(k-2) *(-4.)*(1.+z)  ) ) +
      !(i-2) *(	     ( !k *(2.)*(1.-z)  + !(k-2) *(2.)*(1.+z)   ) );
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
      !(i-2) *(	        ( !k *x  + !(k-1) *(-2.)*x  + !(k-2) *x  ) );
  }

  double  wedgeth::eval_d2phidxdy(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j     *( !k *(2.)*(1.-z)	+ !(k-2) *(2.)*(1.+z)      )+
		!(j-1) *( !k *(-2.)*(1.-z)	+ !(k-2) *(-2.)*(1.+z)     ) ) +
      !(i-1) *( !j     *( !k *(-2.)*(1.-z)	+ !(k-2) *(-2.)*(1.+z)     )+
		!(j-1) *( !k *(2.)*(1.-z)	+ !(k-2) *(2.)*(1.+z)      ) );
  
  }

  double  wedgeth::eval_d2phidydz(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
  
    return
      !i     *( !j     *( !k *(-0.5+2.*t-z)	+ !(k-1) *(2.)*z	+ !(k-2) *(0.5-2.*t-z) 	 ) +
		!(j-1) *( !k *(-2.)*(1.-x-2.*y) 		  	+ !(k-2) *2.*(1.-x-2.*y) ) +
		!(j-2) *( !k *(0.5-2.*y+z)	+ !(k-1) *(-2.)*z       + !(k-2) *(-0.5+2.*y+z)   ) ) +
      !(i-1) *( !j     *( !k *(2.)*x					+ !(k-2) *(-2.)*x	 )+
		!(j-1) *( !k *(-2.)*x					+ !(k-2) *2.*x		 ) );
  }

  double  wedgeth::eval_d2phidzdx(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j     *( !k *(-0.5+2*t-z)   + !(k-1) *(2.)*z  	 + !(k-2) *(0.5-2.*t-z)   ) +
		!(j-1) *( !k *(2.)*y				 + !(k-2) *(-2.)*y	  ) ) +
      !(i-1) *( !j     *( !k *(-2.)*(1-2.*x-y)			 + !(k-2) *2.*(1-2.*x-y)  ) +
		!(j-1) *( !k *(-2.)*y				 + !(k-2) *2.*y	  	  ) ) +
      !(i-2) *(	        ( !k *(0.5-2.*x+z)   + !(k-1) *(-2.)*z	 + !(k-2) *(-0.5+2.*x+z)  ) );
  }

} //end namespace femus


