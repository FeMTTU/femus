/*=========================================================================

  Program: FEMUS
  Module: Tetrahedral
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

  double tet1::eval_phi(const int *I,const double* X) const {
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    return (!i*!j*!k)*(1.-x-y-z)+ !(i-2)*x + !(j-2)*y + !(k-2)*z;
  }

  double tet1::eval_dphidx(const int *I,const double* X) const {
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    return -(!i*!j*!k) + !(i-2);
  }

  double tet1::eval_dphidy(const int *I,const double* X) const {
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    return -(!i*!j*!k) + !(j-2);
  }

  double tet1::eval_dphidz(const int *I,const double* X) const {
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    return -(!i*!j*!k) + !(k-2);
  }

  //************************************************************

  double  tet2::eval_phi(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y+z);
 
    return
      !i     *( !j     *( !k*t*(2.*t-1.) + !(k-1)*4.*z*t + !(k-2)*(-z+2.*z*z) )+
		!(j-1) *( !k*4.*y*t      + !(k-1)*4.*y*z		      )+
		!(j-2) *( !k*(-y+2.*y*y) 				      ) )+
      !(i-1) *( !j     *( !k*4.*x*t      + !(k-1)*4.*x*z		      )+
		!(j-1) *( !k*4.*x*y					      ) )+
      !(i-2) *( !j     *( !k*(-x+2.*x*x)				      ) );
  }

  double  tet2::eval_dphidx(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y+z);
  
    return
      !i    *( !j     *( !k*(-4.*t+1.) + !(k-1)*(-4.)*z )+
	       !(j-1) *( !k*(-4.)*y 		      	) )+
      !(i-1)*( !j     *( !k*4.*(t-x)   + !(k-1)*4.*z    )+
	       !(j-1) *( !k*4.*y			) )+
      !(i-2)*( !j     *( !k*(-1.+4.*x)		      	) );
  }

  double  tet2::eval_dphidy(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y+z);
  
    return
      !i     *( !j     *( !k*(-4.*t+1.) + !(k-1)*(-4.)*z )+
	        !(j-1) *( !k*4.*(t-y)   + !(k-1)*4.*z    )+
	        !(j-2) *( !k*(-1.+4.*y)                  ) )+
      !(i-1) *( !j     *( !k*(-4.)*x	                 )+
	        !(j-1) *( !k*4.*x		         ) );
  }


  double  tet2::eval_dphidz(const int *I,const double* X) const {
    
    const double x=X[0];   const double y=X[1];   const double z=X[2];
    const int i=I[0];      const int j=I[1];      const int k=I[2];
    double t=1.-(x+y+z);
  
    return
      !i     *( !j     *( !k*(-4.*t+1.) + !(k-1)*4.*(t-z) + !(k-2)*(-1+4.*z) )+
		!(j-1) *( !k*(-4.)*y    + !(k-1)*4.*y			     ) )+
      !(i-1) *( !j     *( !k*(-4.)*x    + !(k-1)*4.*x			     ) );
    }

  double  tet2::eval_d2phidx2(const int *I,const double* X) const {
  
    const int i=I[0];  const int j=I[1];  const int k=I[2];
    return
      !i    *( !j *( !k*(4.)   ) )+
      !(i-1)*( !j *( !k*(-8.)  ) )+
      !(i-2)*( !j *( !k*(4.)   ) );   
  }

  double  tet2::eval_d2phidy2(const int *I,const double* X) const {
  
    const int i=I[0];  const int j=I[1];  const int k=I[2];
    return
      !i  *( !j     *( !k*(4.)  )+
	     !(j-1) *( !k*(-8.) )+
	     !(j-2) *( !k*(4.)  ) );
  }


  double  tet2::eval_d2phidz2(const int *I,const double* X) const {
  
    const int i=I[0];  const int j=I[1];  const int k=I[2];
    return
      !i *( !j *( !k*(4.) + !(k-1)*(-8.) + !(k-2)*(4.) ));
  }

  double  tet2::eval_d2phidxdy(const int *I,const double* X) const {
  
    const int i=I[0]; const int j=I[1];  const int k=I[2];
    return
      !i    *( !j     *( !k*(4.)  )+
	       !(j-1) *( !k*(-4.) ) )+
      !(i-1)*( !j     *( !k*(-4.) )+
	       !(j-1) *( !k*(4.)  ) );
  }

  double  tet2::eval_d2phidydz(const int *I,const double* X) const {
 
    const int i=I[0];  const int j=I[1];  const int k=I[2];
    return
      !i   *( !j     *( !k*(4.)   + !(k-1)*(-4.) )+
	      !(j-1) *( !k*(-4.)  + !(k-1)*(4.)	) );
  }

  double  tet2::eval_d2phidzdx(const int *I,const double* X) const {
 
    const int i=I[0];  const int j=I[1];  const int k=I[2];
    return
      !i    *( !j *( !k*(4.)  + !(k-1)*(-4.) ) )+
      !(i-1)*( !j *( !k*(-4.) + !(k-1)*4.    ) );
  
  }

} //end namespace femus


