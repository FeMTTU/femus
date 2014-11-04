/*=========================================================================

  Program: FEMUS
  Module: Hexaedron
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
#include <cmath>


namespace femus {

  double hexpwl::eval_phi(const int *I,const double* x) const {
    return I[0]+ x[0]*I[1] + x[1]*I[2] + x[2]*(1.-I[0])*(1.-I[1])*(1.-I[2]);
  }

  double hexpwl::eval_dphidx(const int *I,const double* x) const {
    return I[1];
  }

  double hexpwl::eval_dphidy(const int *I,const double* x) const {
    return I[2];
  }

  double hexpwl::eval_dphidz(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1])*(1.-I[2]);
  }

  // ************************************************************

  double hex1::eval_phi(const int *I,const double* x) const {
    return lag1(x[0],I[0])*lag1(x[1],I[1])*lag1(x[2],I[2]);
  }

  double hex1::eval_dphidx(const int *I,const double* x) const {
    return dlag1(x[0],I[0])*lag1(x[1],I[1])*lag1(x[2],I[2]);
  }

  double hex1::eval_dphidy(const int *I,const double* x) const {
    return lag1(x[0],I[0])*dlag1(x[1],I[1])*lag1(x[2],I[2]);
  }

  double hex1::eval_dphidz(const int *I,const double* x) const {
    return lag1(x[0],I[0])*lag1(x[1],I[1])*dlag1(x[2],I[2]);
  }

  double hex1::eval_d2phidxdy(const int *I,const double* x) const {
    return dlag1(x[0],I[0])*dlag1(x[1],I[1])*lag1(x[2],I[2]);
  }

  double hex1::eval_d2phidydz(const int *I,const double* x) const {
    return lag1(x[0],I[0])*dlag1(x[1],I[1])*dlag1(x[2],I[2]);
  }

  double hex1::eval_d2phidzdx(const int *I,const double* x) const {
    return dlag1(x[0],I[0])*lag1(x[1],I[1])*dlag1(x[2],I[2]);
  }

  // ************************************************************

  double hex2::eval_phi(const int *I,const double* x) const {
    return lag2(x[0],I[0])*lag2(x[1],I[1])*lag2(x[2],I[2]);
  }

  double hex2::eval_dphidx(const int *I,const double* x) const {
    return dlag2(x[0],I[0])*lag2(x[1],I[1])*lag2(x[2],I[2]);
  }

  double hex2::eval_dphidy(const int *I,const double* x) const {
    return lag2(x[0],I[0])*dlag2(x[1],I[1])*lag2(x[2],I[2]);
  }

  double hex2::eval_dphidz(const int *I,const double* x) const {
    return lag2(x[0],I[0])*lag2(x[1],I[1])*dlag2(x[2],I[2]);
  }

  double hex2::eval_d2phidx2(const int *I,const double* x) const {
    return d2lag2(x[0],I[0])*lag2(x[1],I[1])*lag2(x[2],I[2]);
  }

  double hex2::eval_d2phidy2(const int *I,const double* x) const {
    return lag2(x[0],I[0])*d2lag2(x[1],I[1])*lag2(x[2],I[2]);
  }

  double hex2::eval_d2phidz2(const int *I,const double* x) const {
    return lag2(x[0],I[0])*lag2(x[1],I[1])*d2lag2(x[2],I[2]);
  }

  double hex2::eval_d2phidxdy(const int *I,const double* x) const {
    return dlag2(x[0],I[0])*dlag2(x[1],I[1])*lag2(x[2],I[2]);
  }

  double hex2::eval_d2phidydz(const int *I,const double* x) const {
    return lag2(x[0],I[0])*dlag2(x[1],I[1])*dlag2(x[2],I[2]);
  }

  double hex2::eval_d2phidzdx(const int *I,const double* x) const {
    return dlag2(x[0],I[0])*lag2(x[1],I[1])*dlag2(x[2],I[2]);
  }

  // ************************************************************

  double hexth::eval_phi(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)? 
      th2(x[0],I[0])*th2(x[1],I[1])*th2(x[2],I[2])
      : 
      (-2.+ix*x[0]+jx*x[1]+kx*x[2])*th2(x[0],I[0])*th2(x[1],I[1])*th2(x[2],I[2]);
  }

  double hexth::eval_dphidx(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      dth2(x[0],I[0])*th2(x[1],I[1])*th2(x[2],I[2]) 
      :
      th2(x[1],I[1])*th2(x[2],I[2])*(ix*th2(x[0],I[0]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dth2(x[0],I[0]));
  }

  double hexth::eval_dphidy(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      th2(x[0],I[0])*dth2(x[1],I[1])*th2(x[2],I[2]) 
      :
      th2(x[0],I[0])*th2(x[2],I[2])*(jx*th2(x[1],I[1]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dth2(x[1],I[1]));
  }

  double hexth::eval_dphidz(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      th2(x[0],I[0])*th2(x[1],I[1])*dth2(x[2],I[2]) 
      :
      th2(x[0],I[0])*th2(x[1],I[1])*(kx*th2(x[2],I[2]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dth2(x[2],I[2]));
  }

  double hexth::eval_d2phidx2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      d2th2(x[0],I[0])*th2(x[1],I[1])*th2(x[2],I[2]) 
      :
      th2(x[1],I[1])*th2(x[2],I[2])*( 2.*ix*dth2(x[0],I[0]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*d2th2(x[0],I[0]) );
  }

  double hexth::eval_d2phidy2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      th2(x[0],I[0])*d2th2(x[1],I[1])*th2(x[2],I[2]) 
      :
      th2(x[2],I[2])*th2(x[0],I[0])*( 2.*jx*dth2(x[1],I[1]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*d2th2(x[1],I[1]) );
  }

  double hexth::eval_d2phidz2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      th2(x[0],I[0])*th2(x[1],I[1])*d2th2(x[2],I[2]) 
      :
      th2(x[0],I[0])*th2(x[1],I[1])*( 2.*kx*dth2(x[2],I[2]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*d2th2(x[2],I[2]) );
  }

  double hexth::eval_d2phidxdy(const int *I,const double* x) const {
    const double ix=I[0]-1., jx=I[1]-1., kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      dth2(x[0],I[0])*dth2(x[1],I[1])*th2(x[2],I[2]) 
      :
      th2(x[2],I[2])*( ix*th2(x[0],I[0])*dth2(x[1],I[1]) + 
		       jx*th2(x[1],I[1])*dth2(x[0],I[0]) +
		       (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dth2(x[0],I[0])*dth2(x[1],I[1]));
	      
  }

  double hexth::eval_d2phidydz(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      th2(x[0],I[0])*dth2(x[1],I[1])*dth2(x[2],I[2]) 
      :
      th2(x[0],I[0])*( jx*th2(x[1],I[1])*dth2(x[2],I[2]) + 
		       kx*th2(x[2],I[2])*dth2(x[1],I[1]) +
		       (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dth2(x[1],I[1])*dth2(x[2],I[2]));
	    
  }

  double hexth::eval_d2phidzdx(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      dth2(x[0],I[0])*th2(x[1],I[1])*dth2(x[2],I[2]) 
      :
      th2(x[1],I[1])*( kx*th2(x[2],I[2])*dth2(x[0],I[0]) + 
		       ix*th2(x[0],I[0])*dth2(x[2],I[2]) +
		       (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dth2(x[2],I[2])*dth2(x[0],I[0]));
  }

} //end namespace femus


