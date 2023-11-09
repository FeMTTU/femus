/*=========================================================================

  Program: FEMUS
  Module: Hexaedron
  Authors: Eugenio Aulisa and Sara Calandrini
 
  Copyright (c) FEMTTU
  All rights reserved. 

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/



#include "Hexahedron.hpp"

#include <cmath>


namespace femus {
  
  
  //hex const vectors
  const double hex_lag::Xc[27][3]= {
    {-1, -1, -1},{1, -1, -1},{1, 1, -1},{-1, 1, -1},{-1, -1, 1},{1, -1, 1},{1, 1, 1},{-1, 1, 1},{0, -1, -1},
    {1, 0, -1},{0, 1, -1},{-1, 0, -1},{0, -1, 1},{1, 0, 1},{0, 1, 1},{-1, 0, 1},{-1, -1, 0},
    {1, -1, 0},{1, 1, 0},{-1, 1, 0},{0, -1, 0},{1, 0, 0},{0, 1, 0},{-1, 0, 0},{0, 0, -1},
    {0, 0, 1},{0, 0, 0}
  }; 
    

  const int hex_lag::IND[27][3]= {
    {0, 0, 0},{2, 0, 0},{2, 2, 0},{0, 2, 0},
    {0, 0, 2},{2, 0, 2},{2, 2, 2},{0, 2, 2},
    {1, 0, 0},{2, 1, 0},{1, 2, 0},{0, 1, 0},
    {1, 0, 2},{2, 1, 2},{1, 2, 2},{0, 1, 2},
    {0, 0, 1},{2, 0, 1},{2, 2, 1},{0, 2, 1},
    {1, 0, 1},{2, 1, 1},{1, 2, 1},{0, 1, 1},{1, 1, 0},{1, 1, 2},{1, 1, 1}
  };
  
  const int hex_lag::KVERT_IND[125][2]= {
    {0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},
    {0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},
    {0,4},{1,5},{2,6},{3,7},
    {0,5},{1,6},{2,7},{3,4},{0,2},{4,6},{0,6},
    {0,8},{0,9},{0,10},{0,11},{0,12},{0,13},{0,14},{0,15},{0,16},{0,17},{0,18},{0,19}, //27->38
    {1,8},{1,9},{1,10},       {1,12},{1,13},{1,14},              {1,17},{1,18},        //39->46
	  {2,9},{2,10},{2,11},       {2,13},{2,14},{2,15},              {2,18},{2,19}, //47->54
		{3,10},{3,11},              {3,14},{3,15},                     {3,19}, //55->59
			      {4,12},{4,13},{4,14},{4,15},{4,16},{4,17},{4,18},{4,19}, //60->67
			      {5,12},{5,13},{5,14},              {5,17},{5,18},        //68->72
				      {6,13},{6,14},{6,15},              {6,18},{6,19}, //73->77
					     {7,14},{7,15},                     {7,19}, //78->80
					     {0,20},{0,21},{0,22},{0,23},{0,24},{0,25}, //81->86
					     {1,20},{1,21},{1,22},       {1,24},{1,25}, //87->91
					            {2,21},{2,22},{2,23},{2,24},{2,25}, //92->96
					                   {3,22},{3,23},{3,24},{3,25}, //97->100
					     {4,20},{4,21},{4,22},{4,23}       ,{4,25}, //101->105
					     {5,20},{5,21},{5,22}              ,{5,25}, //106->109
					            {6,21},{6,22},{6,23}       ,{6,25}, //110->113
					                   {7,22},{7,23}       ,{7,25}, //114->116
					     {0,26},{1,26},{2,26},{3,26},{4,26},{5,26},{6,26},{7,26} //117->124
  }; 
  
  
  
  const unsigned hex_lag::fine2CoarseVertexMapping[8][8]={
    {0,8,24,11,16,20,26,23},
    {8,1,9,24,20,17,21,26},
    {24,9,2,10,26,21,18,22},
    {11,24,10,3,23,26,22,19},
    {16,20,26,23,4,12,25,15},
    {20,17,21,26,12,5,13,25},
    {26,21,18,22,25,13,6,14},
    {23,26,22,19,15,25,14,7} };
    
  const unsigned hex_lag::faceDofs[6][9] = { 
    {0, 1, 5, 4, 8, 17, 12, 16, 20},
    {1, 2, 6, 5, 9, 18, 13, 17, 21},
    {2, 3, 7, 6, 10, 19, 14, 18, 22},
    {3, 0, 4, 7, 11, 16, 15, 19, 23},
    {0, 3, 2, 1, 11, 10, 9, 8, 24},
    {4, 5, 6, 7, 12, 13, 14, 15, 25}
  };
  
  //******************************************************************
  
  const double hex_const::X[32][3]= {
    {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
    {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
    {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
    {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
    {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
    {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
    {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
    {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5}
  };

  const int hex_const::IND[4][3]= { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} };

  const int hex_const::KVERT_IND[32][2]= {
    {0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0},
    {0,1},{1,1},{2,1},{3,1},{4,1},{5,1},{6,1},{7,1},
    {0,2},{1,2},{2,2},{3,2},{4,2},{5,2},{6,2},{7,2},
    {0,3},{1,3},{2,3},{3,3},{4,3},{5,3},{6,3},{7,3}
  };
  
  
  double hexpwLinear::eval_phi(const int *I,const double* x) const {
    return (1.-I[0])*(1.-I[1])*(1.-I[2]) + 
	    x[0]*eval_dphidx(I,x) + 
	    x[1]*eval_dphidy(I,x) + 
	    x[2]*eval_dphidz(I,x);
  }

  double hexpwLinear::eval_dphidx(const int *I,const double* x) const {
    return I[0];
  }

  double hexpwLinear::eval_dphidy(const int *I,const double* x) const {
    return I[1];
  }

  double hexpwLinear::eval_dphidz(const int *I,const double* x) const {
    return I[2];
  }

  // ************************************************************

  double HexLinear::eval_phi(const int *I,const double* x) const {
    return lagLinear(x[0],I[0])*lagLinear(x[1],I[1])*lagLinear(x[2],I[2]);
  }

  double HexLinear::eval_dphidx(const int *I,const double* x) const {
    return dlagLinear(x[0],I[0])*lagLinear(x[1],I[1])*lagLinear(x[2],I[2]);
  }

  double HexLinear::eval_dphidy(const int *I,const double* x) const {
    return lagLinear(x[0],I[0])*dlagLinear(x[1],I[1])*lagLinear(x[2],I[2]);
  }

  double HexLinear::eval_dphidz(const int *I,const double* x) const {
    return lagLinear(x[0],I[0])*lagLinear(x[1],I[1])*dlagLinear(x[2],I[2]);
  }

  double HexLinear::eval_d2phidxdy(const int *I,const double* x) const {
    return dlagLinear(x[0],I[0])*dlagLinear(x[1],I[1])*lagLinear(x[2],I[2]);
  }

  double HexLinear::eval_d2phidydz(const int *I,const double* x) const {
    return lagLinear(x[0],I[0])*dlagLinear(x[1],I[1])*dlagLinear(x[2],I[2]);
  }

  double HexLinear::eval_d2phidzdx(const int *I,const double* x) const {
    return dlagLinear(x[0],I[0])*lagLinear(x[1],I[1])*dlagLinear(x[2],I[2]);
  }

  // ************************************************************

  double HexBiquadratic::eval_phi(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_dphidx(const int *I,const double* x) const {
    return dlagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_dphidy(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*dlagBiquadratic(x[1],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_dphidz(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1])*dlagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_d2phidx2(const int *I,const double* x) const {
    return d2lagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_d2phidy2(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*d2lagBiquadratic(x[1],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_d2phidz2(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1])*d2lagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_d2phidxdy(const int *I,const double* x) const {
    return dlagBiquadratic(x[0],I[0])*dlagBiquadratic(x[1],I[1])*lagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_d2phidydz(const int *I,const double* x) const {
    return lagBiquadratic(x[0],I[0])*dlagBiquadratic(x[1],I[1])*dlagBiquadratic(x[2],I[2]);
  }

  double HexBiquadratic::eval_d2phidzdx(const int *I,const double* x) const {
    return dlagBiquadratic(x[0],I[0])*lagBiquadratic(x[1],I[1])*dlagBiquadratic(x[2],I[2]);
  }

  // ************************************************************

  double HexQuadratic::eval_phi(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)? 
      lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2])
      : 
      (-2.+ix*x[0]+jx*x[1]+kx*x[2])*lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2]);
  }

  double HexQuadratic::eval_dphidx(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      dlagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2])*(ix*lagQuadratic(x[0],I[0]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dlagQuadratic(x[0],I[0]));
  }

  double HexQuadratic::eval_dphidy(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      lagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[0],I[0])*lagQuadratic(x[2],I[2])*(jx*lagQuadratic(x[1],I[1]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dlagQuadratic(x[1],I[1]));
  }

  double HexQuadratic::eval_dphidz(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*dlagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*(kx*lagQuadratic(x[2],I[2]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dlagQuadratic(x[2],I[2]));
  }

  double HexQuadratic::eval_d2phidx2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      d2lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2])*( 2.*ix*dlagQuadratic(x[0],I[0]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*d2lagQuadratic(x[0],I[0]) );
  }

  double HexQuadratic::eval_d2phidy2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      lagQuadratic(x[0],I[0])*d2lagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[2],I[2])*lagQuadratic(x[0],I[0])*( 2.*jx*dlagQuadratic(x[1],I[1]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*d2lagQuadratic(x[1],I[1]) );
  }

  double HexQuadratic::eval_d2phidz2(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*d2lagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*( 2.*kx*dlagQuadratic(x[2],I[2]) + (-2.+ix*x[0]+jx*x[1]+kx*x[2])*d2lagQuadratic(x[2],I[2]) );
  }

  double HexQuadratic::eval_d2phidxdy(const int *I,const double* x) const {
    const double ix=I[0]-1., jx=I[1]-1., kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      dlagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1])*lagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[2],I[2])*( ix*lagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1]) + 
		       jx*lagQuadratic(x[1],I[1])*dlagQuadratic(x[0],I[0]) +
		       (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dlagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1]));
	      
  }

  double HexQuadratic::eval_d2phidydz(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      lagQuadratic(x[0],I[0])*dlagQuadratic(x[1],I[1])*dlagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[0],I[0])*( jx*lagQuadratic(x[1],I[1])*dlagQuadratic(x[2],I[2]) + 
		       kx*lagQuadratic(x[2],I[2])*dlagQuadratic(x[1],I[1]) +
		       (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dlagQuadratic(x[1],I[1])*dlagQuadratic(x[2],I[2]));
	    
  }

  double HexQuadratic::eval_d2phidzdx(const int *I,const double* x) const {
    const double ix=I[0]-1.,jx=I[1]-1.,kx=I[2]-1.;
    return (fabs(ix*jx*kx)==0)?
      dlagQuadratic(x[0],I[0])*lagQuadratic(x[1],I[1])*dlagQuadratic(x[2],I[2]) 
      :
      lagQuadratic(x[1],I[1])*( kx*lagQuadratic(x[2],I[2])*dlagQuadratic(x[0],I[0]) + 
		       ix*lagQuadratic(x[0],I[0])*dlagQuadratic(x[2],I[2]) +
		       (-2.+ix*x[0]+jx*x[1]+kx*x[2])*dlagQuadratic(x[2],I[2])*dlagQuadratic(x[0],I[0]));
  }

} //end namespace femus


