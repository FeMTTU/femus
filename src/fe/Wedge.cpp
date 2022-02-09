/*=========================================================================

  Program: FEMUS
  Module: Wedge
  Authors: Eugenio Aulisa and Sara Calandrini
 
  Copyright (c) FEMTTU
  All rights reserved. 

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/



#include "Basis.hpp"


namespace femus {

  //************************************************************

  const double wedge_lag::Xc[21][3]= {
    {0, 0, -1},      {1, 0, -1},      {0, 1, -1},                     //vertici triangoli
    {0, 0, 1},       {1, 0, 1},       {0, 1, 1}, 
    {0.5, 0, -1},    {0.5, 0.5, -1},  {0, 0.5, -1},                   //midpoints triangoli
    {0.5, 0, 1},     {0.5, 0.5, 1},   {0, 0.5, 1},
    {0, 0, 0},       {1, 0, 0} ,      {0, 1, 0},                      //midpoints quadrati
    {0.5, 0, 0},     {0.5, 0.5,0},    {0, 0.5, 0},  //0->17           //facce quadrati
    {1./3.,1./3.,-1},{1./3.,1./3.,1}, {1./3.,1./3.,0} 
  };
    
//     {0.25,0,-1},     {0.25,0.25,-1}, {0,0.25,-1},                    //midpoints triangoli 
//     {0.25,0,0},      {0.25,0.25,0},  {0,0.25,0},                     
//     {0.,0,-0.5},     {0.5,0.,-0.5},  {0,0.5,-0.5},  //18->26         //midpoints quadrati
//     {0.75,0.25,-1},  {0.5,0.25,-1},  {0.75,0,-1},                    //midpoints triangoli
//     {0.75,0.25,0},   {0.5,0.25,0},   {0.75,0,0},
//     {1,0.,-0.5},     {0.5,0.5,-0.5}, //27->34                        //midpoint quadrato + faccia
//     {0,0.75,-1},     {0.25,0.5,-1},  {0.25,0.75,-1},
//     {0,0.75,0},      {0.25,0.5,0},   {0.25,0.75,0},
//     {0,1.,-0.5},                    //34->41
//     {0.25,0,1},      {0.25,0.25,1},  {0,0.25,1},                     //midpoints triangolo alto
//     {0.,0,0.5},      {0.5,0.,0.5},   {0,0.5,0.5},  //42->47          //midpoint quadrato + faccia
//     {0.75,0.25,1},   {0.5,0.25,1},   {0.75,0,1},
//     {1.0,0.,0.5},    {0.5,0.5,0.5},  //48->52
//     {0,0.75,1},      {0.25,0.5,1}, {0.25,0.75,1},
//     {0,1.,0.5},                      //53->56
//     {0.25,0,-0.5},   {0.25,0.25,-0.5},{0,0.25,-0.5},
//     {0.75,0.25,-0.5},{0.5,0.25,-0.5}, {0.75,0.,-0.5},
//     {0,0.75,-0.5},   {0.25,0.5,-0.5}, {0.25,0.75,-0.5}, //57->65
//     {0.25,0,0.5},    {0.25,0.25,0.5}, {0,0.25,0.5},
//     {0.75,0.25,0.5}, {0.5,0.25,0.5},  {0.75,0.,0.5},
//     {0,0.75,0.5},    {0.25,0.5,0.5},  {0.25,0.75,0.5}, //66->74 
//     {1./3.,1./3.,-1},  {1./6.,1./6.,-1},  {2./3.,1./6.,-1},   {1./6.,2./3.,-1},      //nodo centrale triangoli
//     {1./3.,1./3.,0},   {1./6.,1./6.,0},   {2./3.,1./6.,0},    {1./6.,2./3.,0},
//     {1./3.,1./3.,1},   {1./6.,1./6.,1},   {2./3.,1./6.,1},    {1./6.,2./3.,1}, //75->86
//     {1./3.,1./3.,-0.5},{1./6.,1./6.,-0.5},{2./3.,1./6.,-0.5}, {1./6.,2./3.,-0.5},    //nodo centrale wedges
//     {1./3.,1./3.,0.5}, {1./6.,1./6.,0.5}, {2./3.,1./6.,0.5},  {1./6.,2./3.,0.5}, //87->94
//   };

  const int wedge_lag::IND[21][3]= {
    {0, 0, 0},{2, 0, 0},{0, 2, 0},
    {0, 0, 2},{2, 0, 2},{0, 2, 2},
    {1, 0, 0},{1, 1, 0},{0, 1, 0},
    {1, 0, 2},{1, 1, 2},{0, 1, 2},
    {0, 0, 1},{2, 0, 1},{0, 2, 1},
    {1, 0, 1},{1, 1, 1},{0, 1, 1},
    {7, 7, 0},{7, 7, 2},{7, 7, 1},
  };
  
  //           5
//          /|\
//         / | \
//        /  |  \
//      11  19   10
//      /   14    \
//     /     |     \
//    /      |      \
//   3-------9-------4
//   |  17  20  16   |
//   |       2       |
//   |      / \      |
//   |     /   \     |
//  12    / 15  \   13
//   |   8       7   |
//   |  /   18    \  |
//   | /           \ |
//   |/             \|
//   0-------6-------1
  
  
    const int wedge_lag::KVERT_IND[95][2]= { //nuova numerazione
    {0,0},{1,1},{2,2},{4,3},{5,4},{6,5}, 
    {0,1},{1,2},{2,0},{4,4},{5,5},{6,3}, 
    {0,3},{1,4},{2,5}, 
    {0,4},{1,5},{2,3}, //0->17 
    {0,6},{0,7},{0,8},{0,9},{0,10},{0,11},{0,12},{0,13},{0,14}, //18->26 
    {1,7},{1,8},{1,6},{1,10},{1,11},{1,9},{1,13},{1,14}, //27->34 
    {2,8},{2,6},{2,7},{2,11},{2,9},{2,10},{2,14}, //34->41 
    {4,9},{4,10},{4,11},{4,12},{4,13},{4,14}, //42->47 
    {5,10},{5,11},{5,9},{5,13},{5,14}, //48->52
    {6,11},{6,9},{6,10},{6,14}, //53->56 
    {0,15},{0,16},{0,17},
    {1,16},{1,17},{1,15},
    {2,17},{2,15},{2,16}, //57->65 
    {4,15},{4,16},{4,17},
    {5,16},{5,17},{5,15},
    {6,17},{6,15},{6,16}, //66->74 
    {3,18},{0,18},{1,18},{2,18},
    {7,18},{4,18},{5,18},{6,18},
    {7,19},{4,19},{5,19},{6,19}, //75->86 
    {3,20},{0,20},{1,20},{2,20},{7,20},{4,20},{5,20},{6,20}, //87->94
  };
  
  const unsigned wedge_lag::fine2CoarseVertexMapping[8][8]= { //nuova numerazione
    {0,6,8,12,15,17},
    {6,1,7,15,13,16},
    {8,7,2,17,16,14},
    {7,8,6,16,17,15},
    {12,15,17,3,9,11},
    {15,13,16,9,4,10},
    {17,16,14,11,10,5},
    {16,17,15,10,11,9}
  };
  
    const unsigned wedge_lag::faceDofs[5][9] = { 
      {0, 1, 4, 3, 6, 13, 9, 12, 15},
      {1, 2, 5, 4, 7, 14, 10, 13, 16},
      {2, 0, 3, 5, 8, 12, 11, 14, 17},
      {0, 2, 1, 8, 7, 6, 18},
      {3, 4, 5, 9, 10, 11, 19}
  };

  //************************************************************
  
  const double wedge_const::X[32][3]= {
    {1./6.,1./6.,-0.5},{2./3.,1./6.,-0.5},{1./6.,2./3.,-0.5},{1./3.,1./3.,-0.5},{1./6.,1./6.,0.5},{2./3.,1./6.,0.5},{1./6.,2./3.,0.5},{1./3.,1./3.,0.5},
    {1./6.,1./6.,-0.5},{2./3.,1./6.,-0.5},{1./6.,2./3.,-0.5},{1./3.,1./3.,-0.5},{1./6.,1./6.,0.5},{2./3.,1./6.,0.5},{1./6.,2./3.,0.5},{1./3.,1./3.,0.5},
    {1./6.,1./6.,-0.5},{2./3.,1./6.,-0.5},{1./6.,2./3.,-0.5},{1./3.,1./3.,-0.5},{1./6.,1./6.,0.5},{2./3.,1./6.,0.5},{1./6.,2./3.,0.5},{1./3.,1./3.,0.5},
    {1./6.,1./6.,-0.5},{2./3.,1./6.,-0.5},{1./6.,2./3.,-0.5},{1./3.,1./3.,-0.5},{1./6.,1./6.,0.5},{2./3.,1./6.,0.5},{1./6.,2./3.,0.5},{1./3.,1./3.,0.5},
  };
   
  const int wedge_const::IND[4][3]= { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} };

  const int wedge_const::KVERT_IND[32][2]= {
    {0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0},
    {0,1},{1,1},{2,1},{3,1},{4,1},{5,1},{6,1},{7,1},
    {0,2},{1,2},{2,2},{3,2},{4,2},{5,2},{6,2},{7,2},
    {0,3},{1,3},{2,3},{3,3},{4,3},{5,3},{6,3},{7,3}
  };
  
  //************************************************************
  
  double wedgepwLinear::eval_phi(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1])*(1.-I[2]) + 
	   (x[0]-1./3.)*eval_dphidx(I,x) + 
	   (x[1]-1./3.)*eval_dphidy(I,x) + 
	          x[2] *eval_dphidz(I,x);
  }

  double wedgepwLinear::eval_dphidx(const int *I,const double* x) const {
    return I[0];
  }

  double wedgepwLinear::eval_dphidy(const int *I,const double* x) const {
    return I[1];
  }

  double wedgepwLinear::eval_dphidz(const int *I,const double* x) const {
    return I[2];
  }

  //************************************************************
  
  double WedgeLinear::eval_phi(const int *I,const double* x) const {
    return triangleLinear(x[0],x[1],I[0],I[1])*lagLinear(x[2],I[2]);
  }

  double WedgeLinear::eval_dphidx(const int *I,const double* x) const {
    return dtriangleLineardx(x[0],x[1],I[0],I[1])*lagLinear(x[2],I[2]);
  }

  double WedgeLinear::eval_dphidy(const int *I,const double* x) const {
    return dtriangleLineardy(x[0],x[1],I[0],I[1])*lagLinear(x[2],I[2]);
  }

  double WedgeLinear::eval_dphidz(const int *I,const double* x) const {
    return triangleLinear(x[0],x[1],I[0],I[1])*dlagLinear(x[2],I[2]);
  }

  //************************************************************

  double WedgeBiquadratic::eval_phi(const int *I,const double* x) const {
    return triangleBiquadratic(x[0],x[1],I[0],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_dphidx(const int *I,const double* x) const {
    return dtriangleBiquadraticdx(x[0],x[1],I[0],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_dphidy(const int *I,const double* x) const {
    return dtriangleBiquadraticdy(x[0],x[1],I[0],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_dphidz(const int *I,const double* x) const {
    return triangleBiquadratic(x[0],x[1],I[0],I[1])*dlagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2triangleBiquadraticdx2(x[0],x[1],I[0],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_d2phidy2(const int *I,const double* x) const {
    return d2triangleBiquadraticdy2(x[0],x[1],I[0],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_d2phidz2(const int *I,const double* x) const {
    return triangleBiquadratic(x[0],x[1],I[0],I[1])*d2lagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_d2phidxdy(const int *I,const double* x) const {
    return d2triangleBiquadraticdxdy(x[0],x[1],I[0],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_d2phidydz(const int *I,const double* x) const {
    return dtriangleBiquadraticdy(x[0],x[1],I[0],I[1])*dlagBiquadratic(x[2],I[2]);
  }

  double WedgeBiquadratic::eval_d2phidzdx(const int *I,const double* x) const {
    return dtriangleBiquadraticdx(x[0],x[1],I[0],I[1])*dlagBiquadratic(x[2],I[2]);
  }

  //************************************************************
  
  double  WedgeQuadratic::eval_phi(const int *I,const double* X) const {
    
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

  double  WedgeQuadratic::eval_dphidx(const int *I,const double* X) const {
  
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

  double  WedgeQuadratic::eval_dphidy(const int *I,const double* X) const {
  
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

  double  WedgeQuadratic::eval_dphidz(const int *I,const double* X) const {
    
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

  double WedgeQuadratic::eval_d2phidx2(const int *I,const double* X) const {
  
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j  *( !k *(2.)*(1.-z)  + !(k-2) *(2.)*(1.+z)   ) ) +
      !(i-1) *( !j  *( !k *(-4.)*(1.-z) + !(k-2) *(-4.)*(1.+z)  ) ) +
      !(i-2) *(	     ( !k *(2.)*(1.-z)  + !(k-2) *(2.)*(1.+z)   ) );
  }

  double WedgeQuadratic::eval_d2phidy2(const int *I,const double* X) const {
  
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i *( !j     *( !k *(2.) *(1.-z)	+ !(k-2) *(2.) *(1.+z)   ) +
	    !(j-1) *( !k *(-4.)*(1.-z)	+ !(k-2) *(-4.)*(1.+z)   ) +
	    !(j-2) *( !k *(2.) *(1.-z)	+ !(k-2) *(2.) *(1.+z)   ) );
  
  }

  double WedgeQuadratic::eval_d2phidz2(const int *I,const double* X) const {
  
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j     *( !k *t  + !(k-1) *(-2.)*t  + !(k-2) *t  ) +
		!(j-2) *( !k *y  + !(k-1) *(-2.)*y  + !(k-2) *y  ) ) +
      !(i-2) *(	        ( !k *x  + !(k-1) *(-2.)*x  + !(k-2) *x  ) );
  }

  double  WedgeQuadratic::eval_d2phidxdy(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y);
    
    return
      !i     *( !j     *( !k *(2.)*(1.-z)	+ !(k-2) *(2.)*(1.+z)      )+
		!(j-1) *( !k *(-2.)*(1.-z)	+ !(k-2) *(-2.)*(1.+z)     ) ) +
      !(i-1) *( !j     *( !k *(-2.)*(1.-z)	+ !(k-2) *(-2.)*(1.+z)     )+
		!(j-1) *( !k *(2.)*(1.-z)	+ !(k-2) *(2.)*(1.+z)      ) );
  
  }

  double  WedgeQuadratic::eval_d2phidydz(const int *I,const double* X) const {
    
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

  double  WedgeQuadratic::eval_d2phidzdx(const int *I,const double* X) const {
    
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


