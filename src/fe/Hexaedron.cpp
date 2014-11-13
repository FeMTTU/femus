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
  
  //hex const vectors
  const double hex_lag::X[125][3]= {{-1, -1, -1},{1, -1, -1},{1, 1, -1},{-1, 1, -1},{-1, -1, 1},{1, -1, 1},{1, 1, 1},{-1, 1, 1},{0, -1, -1},
					{1, 0, -1},{0, 1, -1},{-1, 0, -1},{0, -1, 1},{1, 0, 1},{0, 1, 1},{-1, 0, 1},{-1, -1, 0},
					{1, -1, 0},{1, 1, 0},{-1, 1, 0},{0, -1, 0},{1, 0, 0},{0, 1, 0},{-1, 0, 0},{0, 0, -1},
					{0, 0, 1},{0, 0, 0},{-0.5, -1, -1},{0, -0.5, -1},{-0.5, 0, -1},{-1, -0.5, -1},{-0.5, -1, 0},{0, -0.5, 0},
					{-0.5, 0, 0},{-1, -0.5, 0},{-1, -1, -0.5},{0, -1, -0.5},{0, 0, -0.5},{-1, 0, -0.5},{0.5, -1, -1},{1, -0.5, -1},
					{0.5, 0, -1},{0.5, -1, 0},{1, -0.5, 0},{0.5, 0, 0},{1, -1, -0.5},{1, 0, -0.5},{1, 0.5, -1},{0.5, 1, -1},
					{0, 0.5, -1},{1, 0.5, 0},{0.5, 1, 0},{0, 0.5, 0},{1, 1, -0.5},{0, 1, -0.5},{-0.5, 1, -1},{-1, 0.5, -1},
					{-0.5, 1, 0},{-1, 0.5, 0},{-1, 1, -0.5},{-0.5, -1, 1},{0, -0.5, 1},{-0.5, 0, 1},{-1, -0.5, 1},{-1, -1, 0.5},
					{0, -1, 0.5},{0, 0, 0.5},{-1, 0, 0.5},{0.5, -1, 1},{1, -0.5, 1},{0.5, 0, 1},{1, -1, 0.5},{1, 0, 0.5},
					{1, 0.5, 1},{0.5, 1, 1},{0, 0.5, 1},{1, 1, 0.5},{0, 1, 0.5},{-0.5, 1, 1},{-1, 0.5, 1},{-1, 1, 0.5},
					{-0.5, -1, -0.5},{0, -0.5, -0.5},{-0.5, 0, -0.5},{-1, -0.5, -0.5},{-0.5, -0.5, -1},{-0.5, -0.5, 0},{0.5, -1, -0.5},{1, -0.5, -0.5},
					{0.5, 0, -0.5},{0.5, -0.5, -1},{0.5, -0.5, 0},{1, 0.5, -0.5},{0.5, 1, -0.5},{0, 0.5, -0.5},{0.5, 0.5, -1},{0.5, 0.5, 0},
					{-0.5, 1, -0.5},{-1, 0.5, -0.5},{-0.5, 0.5, -1},{-0.5, 0.5, 0},{-0.5, -1, 0.5},{0, -0.5, 0.5},{-0.5, 0, 0.5},{-1, -0.5, 0.5},
					{-0.5, -0.5, 1},{0.5, -1, 0.5},{1, -0.5, 0.5},{0.5, 0, 0.5},{0.5, -0.5, 1},{1, 0.5, 0.5},{0.5, 1, 0.5},{0, 0.5, 0.5},
					{0.5, 0.5, 1},{-0.5, 1, 0.5},{-1, 0.5, 0.5},{-0.5, 0.5, 1},{-0.5, -0.5, -0.5},{0.5, -0.5, -0.5},{0.5, 0.5, -0.5},{-0.5, 0.5, -0.5},
					{-0.5, -0.5, 0.5},{0.5, -0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5}
					};

  const int hex_lag::IND[27][3]= {{0, 0, 0},{2, 0, 0},{2, 2, 0},{0, 2, 0},
				      {0, 0, 2},{2, 0, 2},{2, 2, 2},{0, 2, 2},
				      {1, 0, 0},{2, 1, 0},{1, 2, 0},{0, 1, 0},
				      {1, 0, 2},{2, 1, 2},{1, 2, 2},{0, 1, 2},
				      {0, 0, 1},{2, 0, 1},{2, 2, 1},{0, 2, 1},
				      {1, 0, 1},{2, 1, 1},{1, 2, 1},{0, 1, 1},{1, 1, 0},{1, 1, 2},{1, 1, 1}
				      };
  
  const int hex_lag::KVERT_IND[125][2]= {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7}
					    ,{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4}
					    ,{0,4},{1,5},{2,6},{3,7}
					    ,{0,5},{1,6},{2,7},{3,4},{0,2},{4,6},{0,6}
					    ,{0,8},{0,9},{0,10},{0,11},{0,12},{0,13},{0,14},{0,15},{0,16},{0,17},{0,18},{0,19} //27->38
					    ,{1,8},{1,9},{1,10},       {1,12},{1,13},{1,14},              {1,17},{1,18}        //39->46
					    ,      {2,9},{2,10},{2,11},       {2,13},{2,14},{2,15},              {2,18},{2,19} //47->54
					    ,            {3,10},{3,11},              {3,14},{3,15},                     {3,19} //55->59
					    ,                          {4,12},{4,13},{4,14},{4,15},{4,16},{4,17},{4,18},{4,19} //60->67
					    ,                          {5,12},{5,13},{5,14},              {5,17},{5,18}        //68->72
					    ,                                 {6,13},{6,14},{6,15},              {6,18},{6,19} //73->77
					    ,                                        {7,14},{7,15},                     {7,19} //78->80
					    ,{0,20},{0,21},{0,22},{0,23},{0,24},{0,25} //81->86
					    ,{1,20},{1,21},{1,22},       {1,24},{1,25} //87->91
					    ,       {2,21},{2,22},{2,23},{2,24},{2,25} //92->96
					    ,              {3,22},{3,23},{3,24},{3,25} //97->100
					    ,{4,20},{4,21},{4,22},{4,23}       ,{4,25} //101->105
					    ,{5,20},{5,21},{5,22}              ,{5,25} //106->109
					    ,       {6,21},{6,22},{6,23}       ,{6,25} //110->113
					    ,              {7,22},{7,23}       ,{7,25} //114->116
					    ,{0,26},{1,26},{2,26},{3,26},{4,26},{5,26},{6,26},{7,26} //117->124
					    }; 
  
  
  
  const double hex_const::X[32][3]= {{-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
					 {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
					 {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
					 {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
					 {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
					 {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
					 {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
					 {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5}
					 };

  const int hex_const::IND[4][3]= {{1,0,0},{0, 1,0},{0, 0,1},{0,0,0}};

  const int hex_const::KVERT_IND[32][2]= {{0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0},
					      {0,1},{1,1},{2,1},{3,1},{4,1},{5,1},{6,1},{7,1},
					      {0,2},{1,2},{2,2},{3,2},{4,2},{5,2},{6,2},{7,2},
					      {0,3},{1,3},{2,3},{3,3},{4,3},{5,3},{6,3},{7,3}
					      };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

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


