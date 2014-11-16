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

  //************************************************************

  const double wedge_lag::X[75][3]= {
    {0, 0, -1},      {1, 0, -1},     {0, 1, -1},
    {0, 0, 1},       {1, 0, 1},      {0, 1, 1},
    {0.5, 0, -1},    {0.5, 0.5, -1}, {0, 0.5, -1},
    {0.5, 0, 1},     {0.5, 0.5, 1},  {0, 0.5, 1},
    {0, 0, 0},       {1, 0, 0} ,     {0, 1, 0},
    {0.5, 0, 0},     {0.5, 0.5,0},   {0, 0.5, 0},  //0->17
    {0.25,0,-1},     {0.25,0.25,-1}, {0,0.25,-1},
    {0.25,0,0},      {0.25,0.25,0},  {0,0.25,0},
    {0.,0,-0.5},     {0.5,0.,-0.5},  {0,0.5,-0.5},  //18->26
    {0.75,0.25,-1},  {0.5,0.25,-1},  {0.75,0,-1},
    {0.75,0.25,0},   {0.5,0.25,0},   {0.75,0,0},
    {1,0.,-0.5},     {0.5,0.5,-0.5}, //27->34
    {0,0.75,-1},     {0.25,0.5,-1},  {0.25,0.75,-1},
    {0,0.75,0},      {0.25,0.5,0},   {0.25,0.75,0},
    {0,1.,-0.5},                    //34->41
    {0.25,0,1},      {0.25,0.25,1},  {0,0.25,1},
    {0.,0,0.5},      {0.5,0.,0.5},   {0,0.5,0.5},  //42->47
    {0.75,0.25,1},   {0.5,0.25,1},   {0.75,0,1},
    {1.0,0.,0.5},    {0.5,0.5,0.5},  //48->52
    {0,0.75,1},      {0.25,0.5,1}, {0.25,0.75,1},
    {0,1.,0.5},                      //53->56
    {0.25,0,-0.5},   {0.25,0.25,-0.5},{0,0.25,-0.5},
    {0.75,0.25,-0.5},{0.5,0.25,-0.5}, {0.75,0.,-0.5},
    {0,0.75,-0.5},   {0.25,0.5,-0.5}, {0.25,0.75,-0.5},// 57->65
    {0.25,0,0.5},    {0.25,0.25,0.5}, {0,0.25,0.5},
    {0.75,0.25,0.5}, {0.5,0.25,0.5},  {0.75,0.,0.5},
    {0,0.75,0.5},    {0.25,0.5,0.5},  {0.25,0.75,0.5}// 66->74
  };

  const int wedge_lag::IND[18][3]= {
    {0, 0, 0},{2, 0, 0},{0, 2, 0},
    {0, 0, 2},{2, 0, 2},{0, 2, 2},
    {1, 0, 0},{1, 1, 0},{0, 1, 0},
    {1, 0, 2},{1, 1, 2},{0, 1, 2},
    {0, 0, 1},{2, 0, 1},{0, 2, 1},
    {1, 0, 1},{1, 1, 1},{0, 1, 1}
  };

  const int wedge_lag::KVERT_IND[75][2]= {
    {0,0},{1,0},{2,0},{4,3},{5,3},{6,3},
    {0,1},{1,1},{2,1},{4,4},{5,4},{6,4},
    {0,3},{1,3},{2,3},
    {0,4},{1,4},{2,4},
    {0,6},{0,7},{0,8},{0,9},{0,10},{0,11},{0,12},{0,13},{0,14},
    {1,6},{1,7},{1,8},{1,9},{1,10},{1,11},{1,12},{1,13},
    {2,6},{2,7},{2,8},{2,9},{2,10},{2,11},{2,12},
    {4,9},{4,10},{4,11},{4,12},{4,13},{4,14},
    {5,9},{5,10},{5,11},{5,12},{5,13},
    {6,9},{6,10},{6,11},{6,12},
    {0,15},{0,16},{0,17},
    {1,15},{1,16},{1,17},
    {2,15},{2,16},{2,17},
    {4,15},{4,16},{4,17},
    {5,15},{5,16},{5,17},
    {6,15},{6,16},{6,17}
  };

  //************************************************************

  const double wedge_const::X[32][3]= {
    {0.166666666667,0.166666666667,-0.5},{0.666666666667,0.166666666667,-0.5},{0.166666666667,0.666666666667,-0.5},{0.333333333333,0.333333333333,-0.5},
    {0.166666666667,0.166666666667, 0.5},{0.666666666667,0.166666666667, 0.5},{0.166666666667,0.666666666667, 0.5},{0.333333333333,0.333333333333, 0.5},
    {0.166666666667,0.166666666667,-0.5},{0.666666666667,0.166666666667,-0.5},{0.166666666667,0.666666666667,-0.5},{0.333333333333,0.333333333333,-0.5},
    {0.166666666667,0.166666666667, 0.5},{0.666666666667,0.166666666667, 0.5},{0.166666666667,0.666666666667, 0.5},{0.333333333333,0.333333333333, 0.5},
    {0.166666666667,0.166666666667,-0.5},{0.666666666667,0.166666666667,-0.5},{0.166666666667,0.666666666667,-0.5},{0.333333333333,0.333333333333,-0.5},
    {0.166666666667,0.166666666667, 0.5},{0.666666666667,0.166666666667, 0.5},{0.166666666667,0.666666666667, 0.5},{0.333333333333,0.333333333333, 0.5},
    {0.166666666667,0.166666666667,-0.5},{0.666666666667,0.166666666667,-0.5},{0.166666666667,0.666666666667,-0.5},{0.333333333333,0.333333333333,-0.5},
    {0.166666666667,0.166666666667, 0.5},{0.666666666667,0.166666666667, 0.5},{0.166666666667,0.666666666667, 0.5},{0.333333333333,0.333333333333, 0.5}
  };

  const int wedge_const::IND[4][3]= { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} };

  const int wedge_const::KVERT_IND[32][2]= {
    {0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0},
    {0,1},{1,1},{2,1},{3,1},{4,1},{5,1},{6,1},{7,1},
    {0,2},{1,2},{2,2},{3,2},{4,2},{5,2},{6,2},{7,2},
    {0,3},{1,3},{2,3},{3,3},{4,3},{5,3},{6,3},{7,3}
  };
  
  //************************************************************
  
  double wedgepwl::eval_phi(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1])*(1.-I[2]) + 
	    x[0]*eval_dphidx(I,x) + 
	    x[1]*eval_dphidy(I,x) + 
	    x[2]*eval_dphidz(I,x);
  }

  double wedgepwl::eval_dphidx(const int *I,const double* x) const {
    return I[0];
  }

  double wedgepwl::eval_dphidy(const int *I,const double* x) const {
    return I[1];
  }

  double wedgepwl::eval_dphidz(const int *I,const double* x) const {
    return I[2];
  }

  //************************************************************
  
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


