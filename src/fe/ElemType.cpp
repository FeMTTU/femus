/*=========================================================================

 Program: FEMUS
 Module: ElemType
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

#include "GaussPoints.hpp"
#include "ElemType.hpp"
#include "Elem.hpp"


namespace femus {



unsigned elem_type::_refindex=1;

//hex const vectors
const double HEX_X[125][3]= {{-1, -1, -1},{1, -1, -1},{1, 1, -1},{-1, 1, -1},{-1, -1, 1},{1, -1, 1},{1, 1, 1},{-1, 1, 1},{0, -1, -1},
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

const int HEX_IND[27][3]= {{0, 0, 0},{2, 0, 0},{2, 2, 0},{0, 2, 0},
  {0, 0, 2},{2, 0, 2},{2, 2, 2},{0, 2, 2},
  {1, 0, 0},{2, 1, 0},{1, 2, 0},{0, 1, 0},
  {1, 0, 2},{2, 1, 2},{1, 2, 2},{0, 1, 2},
  {0, 0, 1},{2, 0, 1},{2, 2, 1},{0, 2, 1},
  {1, 0, 1},{2, 1, 1},{1, 2, 1},{0, 1, 1},{1, 1, 0},{1, 1, 2},{1, 1, 1}
};

const int HEX_KVERT_IND[125][2]= {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7}
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

const double HEX_X0[32][3]= {{-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
  {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
  {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
  {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
  {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
  {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5},
  {-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5, 0.5,-0.5},{-0.5, 0.5,-0.5},
  {-0.5,-0.5, 0.5},{0.5,-0.5, 0.5},{0.5, 0.5, 0.5},{-0.5, 0.5, 0.5}
};

const int HEX_IND0[4][3]= {{1,0,0},{0, 1,0},{0, 0,1},{0,0,0}};

const int HEX_KVERT_IND0[32][2]= {{0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0},
  {0,1},{1,1},{2,1},{3,1},{4,1},{5,1},{6,1},{7,1},
  {0,2},{1,2},{2,2},{3,2},{4,2},{5,2},{6,2},{7,2},
  {0,3},{1,3},{2,3},{3,3},{4,3},{5,3},{6,3},{7,3}
};

// wedge const vectors
const double WEDGE_X[75][3]= {{0, 0, -1},      {1, 0, -1},     {0, 1, -1},
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

const int WEDGE_IND[18][3]= {{0, 0, 0},{2, 0, 0},{0, 2, 0},
  {0, 0, 2},{2, 0, 2},{0, 2, 2},
  {1, 0, 0},{1, 1, 0},{0, 1, 0},
  {1, 0, 2},{1, 1, 2},{0, 1, 2},
  {0, 0, 1},{2, 0, 1},{0, 2, 1},
  {1, 0, 1},{1, 1, 1},{0, 1, 1}
};

const int WEDGE_KVERT_IND[75][2]= {{0,0},{1,0},{2,0},{4,3},{5,3},{6,3},
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

// tet const vectors
const double TET_X[35][3]= {{0, 0, 0},      {1, 0, 0},       {0, 1, 0},   {0, 0, 1}, //0->4
  {0.5, 0, 0},    {0.5, 0.5, 0},   {0, 0.5, 0},
  {0.,  0, 0.5},  {0.5, 0., 0.5},  {0, 0.5, 0.5}, //5->9
  {0.25,0,0},     {0.25,0.25,0.},  {0,0.25,0},
  {0,0,0.25},     {0.25,0,0.25},   {0,0.25,0.25}, //10->15
  {0.75,0.25,0},  {0.5,0.25,0},    {0.75,0,0},
  {0.75,0,0.25},  {0.5,0.25,0.25}, {0.5,0,0.25},  //16->21
  {0,0.75,0},     {0.25,0.5,0},    {0.25,0.75,0},
  {0,0.75,0.25},  {0,0.5,0.25},    {0.25,0.5,0.25},//22->27
  {0.25,0,0.5},   {0.25,0.25,0.5}, {0,0.25,0.5},
  {0,0,0.75},     {0.25,0,0.75},   {0,0.25,0.75}, //28->33
  {0.25,0.25,0.25}
};//34

const int TET_IND[10][3]= {{0, 0, 0},{2, 0, 0},{0, 2, 0},{0, 0, 2},
  {1, 0, 0},{1, 1, 0},{0, 1, 0},
  {0, 0, 1},{1, 0, 1},{0, 1, 1}
};

const int TET_KVERT_IND[35][2]= {{0,0},{1,0},{2,0},{3,3},
  {0,1},{1,1},{2,1},
  {0,3},{1,3},{2,3},
  {0,4},{0,5},{0,6},{0,7},{0,8},{0,9},
  {1,4},{1,5},{1,6},{1,7},{1,8},{1,9},
  {2,4},{2,5},{2,6},{2,7},{2,8},{2,9},
  {3,4},{3,5},{3,6},{3,7},{3,8},{3,9},
  {4,5}
};

// square const vectors
const double QUAD_X[25][3]= {{-1,-1},{1,-1},{1, 1},{-1, 1},
  { 0,-1},{1, 0},{0, 1},{-1, 0},{0, 0},
  {-0.5,-1},{0,-0.5},{-0.5,0},{-1,-0.5},
  { 0.5,-1},{1,-0.5},{ 0.5,0},
  { 1, 0.5},{0.5, 1},{0, 0.5},
  {-0.5, 1},{-1,0.5},
  {-0.5,-0.5},{0.5, -0.5},{0.5, 0.5},{-0.5, 0.5}
};


const int QUAD_IND[9][3]= {{0, 0},{2, 0},{2, 2},{0, 2},
  {1, 0},{2, 1},{1, 2},{0, 1},
  {1, 1}
};


const int QUAD_KVERT_IND[25][2]= {
  {0,0},{1,1},{2,2},{3,3},
  {0,1},{1,2},{2,3},{3,0},{0,2},
  {0,4},{0,5},{0,6},{0,7},
  {1,4},{1,5},{1,6},
  {2,5},{2,6},{2,7},
  {3,6},{3,7},
  {0,8},{1,8},{2,8},{3,8}
};


const double QUAD_X0[12][3]= {{-0.5,-0.5},{0.5, -0.5},{0.5, 0.5},{-0.5, 0.5},
  {-0.5,-0.5},{0.5, -0.5},{0.5, 0.5},{-0.5, 0.5},
  {-0.5,-0.5},{0.5, -0.5},{0.5, 0.5},{-0.5, 0.5}
};
const int QUAD_IND0[3][3]= {{1, 0},{0, 1},{0, 0}};

const int QUAD_KVERT_IND0[12][2]= {{0,0},{1,0},{2,0},{3,0},
  {0,1},{1,1},{2,1},{3,1},
  {0,2},{1,2},{2,2},{3,2}
};

// triangle const vectors
const double TRI_X[15][3]= {{0, 0},      {1, 0},      {0, 1},
  {0.5, 0},    {0.5, 0.5},  {0, 0.5},
  {0.25,0},    {0.25,0.25}, {0,0.25},
  {0.75,0.25}, {0.5,0.25},  {0.75,0},
  {0,0.75},    {0.25,0.5},  {0.25,0.75}
};

const int TRI_IND[6][3]= {{0, 0},{2, 0},{0, 2},
  {1, 0},{1, 1},{0, 1}
};


const int TRI_KVERT_IND[15][2]= {{0,0},{1,0},{2,0},
  {0,1},{1,1},{2,1},
  {0,3},{0,4},{0,5},
  {1,3},{1,4},{1,5},
  {2,3},{2,4},{2,5}
};

// line const vectors
const double LINE_X[5][3]= {{-1},{1},{0},{-0.5},{0.5}};


const int LINE_IND[3][3]= {{0},{2},{1}};


const int LINE_KVERT_IND[5][2]= {{0,0},{1,1},{0,1},{0,2},{1,2}};

const double *GaussLine[5]= {
  GaussLine0[0],
  GaussLine1[0],
  GaussLine2[0],
  GaussLine3[0],
  GaussLine4[0]
};

const double *GaussSquare[5]= {
  GaussSquare0[0],
  GaussSquare1[0],
  GaussSquare2[0],
  GaussSquare3[0],
  GaussSquare4[0]
};

const double *GaussTriangle[5]= {
  GaussTriangle0[0],
  GaussTriangle1[0],
  GaussTriangle2[0],
  GaussTriangle3[0],
  GaussTriangle4[0]
};

const double *GaussCube[5]= {
  GaussCube0[0],
  GaussCube1[0],
  GaussCube2[0],
  GaussCube3[0],
  GaussCube4[0]
};

const double *GaussWedge[5]= {
  GaussWedge0[0],
  GaussWedge1[0],
  GaussWedge2[0],
  GaussWedge3[0],
  GaussWedge4[0]
};

const double *GaussTetrahedra[5]= {
  GaussTetrahedra0[0],
  GaussTetrahedra1[0],
  GaussTetrahedra2[0],
  GaussTetrahedra3[0],
  GaussTetrahedra4[0]
};


using std::cout;
using std::endl;

elem_type::~elem_type() {
  delete [] X;
  delete [] KVERT_IND;
  delete [] IND;

  delete [] prol_val;
  delete [] prol_ind;
  delete [] mem_prol_val;
  delete [] mem_prol_ind;

  delete [] phi;
  delete [] phi_memory;
  delete [] dphidxi;
  delete [] dphidxi_memory;
  delete [] dphideta;
  delete [] dphideta_memory;
  delete [] dphidzeta;
  delete [] dphidzeta_memory;
};

elem_type::elem_type(const char *solid, const char *order, const char *order_gauss) {
  //Gauss Point order
  int gauss_order;
  if (!strcmp(order_gauss,"zero")  || !strcmp(order_gauss,"first")) {
    gauss_order=0;
  } else if (!strcmp(order_gauss,"second") || !strcmp(order_gauss,"third") ) {
    gauss_order=1;
  } else if (!strcmp(order_gauss,"fourth") || !strcmp(order_gauss,"fifth") ) {
    gauss_order=2;
  } else if (!strcmp(order_gauss,"sixth")  || !strcmp(order_gauss,"seventh") ) {
    gauss_order=3;
  } else if (!strcmp(order_gauss,"eighth") || !strcmp(order_gauss,"ninth") ) {
    gauss_order=4;
  } else {
    cout<<order_gauss<<"is not a valid option for the Gauss points of"<<solid<<endl;
    exit(0);
  }

  if (!strcmp(order,"linear")) {
    SolType_=0;   
  } else if (!strcmp(order,"quadratic")) {
    SolType_=1;   
  } else if (!strcmp(order,"biquadratic")) {
    SolType_=2;   
  } else if (!strcmp(order,"constant")) {
    SolType_=3;   
  } else if (!strcmp(order,"disc_linear")) {
    SolType_=4;   
  }
  
  if (!strcmp(solid,"hex")) {//HEX
    Jacobian_ptr=&elem_type::Jacobian3D;
    
    ncf_[0]=8;
    ncf_[1]=20;
    ncf_[2]=27;
    
    if (!strcmp(order,"linear")) {
      type_=0;
      nc_=8;
      nf_=27;
      pt_basis = & hex_1;
    } else if (!strcmp(order,"quadratic")) {
      type_=1;
      nc_=20;
      nf_=81;
      pt_basis = & hex_th;
    } else if (!strcmp(order,"biquadratic")) {
      type_=2;
      nc_=27;
      nf_=125;
      pt_basis = & hex_2;
    } else if (!strcmp(order,"constant")) {
      type_=17;
      nc_=1;
      nf_=8;
      pt_basis = & hex_0;
    } else if (!strcmp(order,"disc_linear")) {
      type_=18;
      nc_=4;
      nf_=32;
      pt_basis = & hex_pwl;
    } else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }
    if (type_<15) {
      IND=new const int * [nc_];
      for (int i=0; i<nc_; i++)
        IND[i]=HEX_IND[i];
      KVERT_IND=new const int * [nf_];
      X=new const double * [nf_];
      for (int i=0; i<nf_; i++) {
        KVERT_IND[i]=HEX_KVERT_IND[i];
        X[i]=HEX_X[i];
      }
    } else {
      IND=new const int * [nc_];
      for (int i=0; i<nc_; i++)
        IND[i]=HEX_IND0[i];
      KVERT_IND=new const int * [nf_];
      X=new const double * [nf_];
      for (int i=0; i<nf_; i++) {
        KVERT_IND[i]=HEX_KVERT_IND0[i];
        X[i]=HEX_X0[i];
      }
    }
    GaussWeight=GaussCube[gauss_order];
    GaussPoints=GaussPointsCube[gauss_order];
  } else if (!strcmp(solid,"wedge")) { //WEDGE
    Jacobian_ptr=&elem_type::Jacobian3D;
         
    ncf_[0]=6;
    ncf_[1]=15;
    ncf_[2]=18;
       
    if (!strcmp(order,"linear")) {
      type_=3;
      nc_=6;
      nf_=18;
      pt_basis = & wedge_1;
    } else if (!strcmp(order,"quadratic")) {
      type_=4;
      nc_=15;
      nf_=57;
      pt_basis = & wedge_th;
    } else if (!strcmp(order,"biquadratic")) {
      type_=5;
      nc_=18;
      nf_=75;
      pt_basis = & wedge_2;
    } else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }
    //
    IND=new const int * [nc_];
    for (int i=0; i<nc_; i++)
      IND[i]=WEDGE_IND[i];

    KVERT_IND=new const int * [nf_];
    X=new const double * [nf_];
    for (int i=0; i<nf_; i++) {
      KVERT_IND[i]=WEDGE_KVERT_IND[i];
      X[i]=WEDGE_X[i];
    }
    GaussWeight=GaussWedge[gauss_order];
    GaussPoints=GaussPointsWedge[gauss_order];
  } else if (!strcmp(solid,"tet")) { //TETRAHEDRA
    Jacobian_ptr=&elem_type::Jacobian3D;
         
    ncf_[0]=4;
    ncf_[1]=10;
    ncf_[2]=10;
    
    if (!strcmp(order,"linear")) {
      type_=6;
      nc_=4;
      nf_=10;
      pt_basis = & tet_1;
    } else if (!strcmp(order,"biquadratic")) {
      type_=7;
      nc_=10;
      nf_=35;
      pt_basis = & tet_2;
    } else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }

    IND=new const int * [nc_];
    for (int i=0; i<nc_; i++)
      IND[i]=TET_IND[i];

    KVERT_IND=new const int * [nf_];
    X=new const double * [nf_];
    for (int i=0; i<nf_; i++) {
      KVERT_IND[i]=TET_KVERT_IND[i];
      X[i]=TET_X[i];
    }
    GaussWeight=GaussTetrahedra[gauss_order];
    GaussPoints=GaussPointsTetrahedra[gauss_order];
  } else if (!strcmp(solid,"quad")) { //QUAD
    Jacobian_ptr=&elem_type::Jacobian2D;
    Jacobian_sur_ptr=&elem_type::JacobianSur2D;
    ncf_[0]=4;
    ncf_[1]=8;
    ncf_[2]=9;
    if (!strcmp(order,"linear")) {
      type_=8;
      nc_=4;
      nf_=9;
      pt_basis = & quad_1;
    } else if (!strcmp(order,"quadratic")) {
      type_=9;
      nc_=8;
      nf_=21;
      pt_basis = & quad_th;
    } else if (!strcmp(order,"biquadratic")) {
      type_=10;
      nc_=9;
      nf_=25;
      pt_basis = & quad_2;
    } else if (!strcmp(order,"constant")) {
      type_=15;
      nc_=1;
      nf_=4;
      pt_basis = & quad_0;
    } else if (!strcmp(order,"disc_linear")) {
      type_=16;
      nc_=3;
      nf_=12;
      pt_basis = & quad_pwl;
    } else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }
    if (type_<15) {
      IND=new const int * [nc_];
      for (int i=0; i<nc_; i++)
        IND[i]=QUAD_IND[i];
      KVERT_IND=new const int * [nf_];
      X=new const double * [nf_];
      for (int i=0; i<nf_; i++) {
        KVERT_IND[i]=QUAD_KVERT_IND[i];
        X[i]=QUAD_X[i];
      }
    } else {
      IND=new const int * [nc_];
      for (int i=0; i<nc_; i++)
        IND[i]=QUAD_IND0[i];
      KVERT_IND=new const int * [nf_];
      X=new const double * [nf_];
      for (int i=0; i<nf_; i++) {
        KVERT_IND[i]=QUAD_KVERT_IND0[i];
        X[i]=QUAD_X0[i];
      }
    }
    GaussWeight=GaussSquare[gauss_order];
    GaussPoints=GaussPointsSquare[gauss_order];
  } else if (!strcmp(solid,"tri")) { //TRIANGLE
    Jacobian_ptr=&elem_type::Jacobian2D;
    Jacobian_sur_ptr=&elem_type::JacobianSur2D;
    
    ncf_[0]=3;
    ncf_[1]=6;
    ncf_[2]=6;
    
    if (!strcmp(order,"linear")) {
      type_=11;
      nc_=3;
      nf_=6;
      pt_basis = & tri_1;
    } else if (!strcmp(order,"biquadratic")) {
      type_=12;
      nc_=6;
      nf_=15;
      pt_basis = & tri_2;
    } else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }

    IND=new const int * [nc_];
    for (int i=0; i<nc_; i++)
      IND[i]=TRI_IND[i];

    KVERT_IND=new const int * [nf_];
    X=new const double * [nf_];
    for (int i=0; i<nf_; i++) {
      KVERT_IND[i]=TRI_KVERT_IND[i];
      X[i]=TRI_X[i];
    }
    GaussWeight=GaussTriangle[gauss_order];
    GaussPoints=GaussPointsTriangle[gauss_order];
  }

  else if (!strcmp(solid,"line")) { //line
    Jacobian_ptr=&elem_type::Jacobian1D;
    Jacobian_sur_ptr=&elem_type::JacobianSur1D;
    
    ncf_[0]=2;
    ncf_[1]=3;
    ncf_[2]=3;
    
    if (!strcmp(order,"linear")) {
      type_=13;
      nc_=2;
      nf_=3;
      pt_basis = & line_1;
    } else if (!strcmp(order,"biquadratic")) {
      type_=14;
      nc_=3;
      nf_=5;
      pt_basis = & line_2;
    } else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }

    IND=new const int * [nc_];
    for (int i=0; i<nc_; i++)
      IND[i]=LINE_IND[i];

    KVERT_IND=new const int * [nf_];
    X=new const double * [nf_];
    for (int i=0; i<nf_; i++) {
      KVERT_IND[i]=LINE_KVERT_IND[i];
      X[i]=LINE_X[i];
    }
    GaussWeight=GaussLine[gauss_order];
    GaussPoints=GaussPointsLine[gauss_order];
  } else {
    cout<<solid<<" is not a valid option"<<endl;
    exit(0);
  }

  int counter=0;
  for (int i=0; i<nf_; i++) {
    for (int j=0; j<nc_; j++) {
      double phi=pt_basis->eval_phi(IND[j],X[i]);
      if (type_==16) {
        if (i/4==1) phi=pt_basis->eval_dphidx(IND[j],X[i]);
        else if (i/4==2) phi=pt_basis->eval_dphidy(IND[j],X[i]);
      }
      if (phi!=0)
        counter++;
    }
  }
  double *pt_d;
  int *pt_i;

  prol_val=new double * [nf_+1];
  prol_ind=new int * [nf_+1];
  mem_prol_val=new double [counter];
  mem_prol_ind=new int [counter];

  pt_d=mem_prol_val;
  pt_i= mem_prol_ind;
  for (int i=0; i<nf_; i++) {
    prol_val[i]=pt_d;
    prol_ind[i]=pt_i;
    for (int j=0; j<nc_; j++) {
      double phi=pt_basis->eval_phi(IND[j],X[i]);
      if (type_==16) {
        if (i/4==1)
          phi=pt_basis->eval_dphidx(IND[j],X[i])/2.;
        else if (i/4==2)
          phi=pt_basis->eval_dphidy(IND[j],X[i])/2.;
      } else if (type_==18) {
        if (i/8==1)
          phi=pt_basis->eval_dphidx(IND[j],X[i])/2.;
        else if (i/8==2)
          phi=pt_basis->eval_dphidy(IND[j],X[i])/2.;
        else if (i/8==3)
          phi=pt_basis->eval_dphidz(IND[j],X[i])/2.;
      }
      if (phi!=0) {
        *(pt_d++)=phi;
        *(pt_i++)=j;
      }
    }
  }
  prol_val[nf_]=pt_d;
  prol_ind[nf_]=pt_i;


//   rest_val=new double * [nc_+1];
//   rest_ind=new int * [nc_+1];
//   mem_rest_val=new double [counter];
//   mem_rest_ind=new int [counter];
//
//   pt_d=mem_rest_val;
//   pt_i=mem_rest_ind;
//
//   for(int i=0;i<nc_;i++){
//     rest_val[i]=pt_d;
//     rest_ind[i]=pt_i;
//     for (int j=0;j<nf_;j++){
//       double phi=pt_basis->eval_phi(IND[i],X[j]);
//       if(type_==16){
// 	if(j/4==1) phi=pt_basis->eval_dphidx(IND[i],X[j]);
// 	else if(j/4==2) phi=pt_basis->eval_dphidy(IND[i],X[j]);
//       }
//       if(phi!=0){
// 	*(pt_d++)=phi;
// 	*(pt_i++)=j;
//       }
//     }
//   }
//   rest_val[nc_]=pt_d;
//   rest_ind[nc_]=pt_i;

  phi= new double*[GaussPoints];
  dphidxi  = new double*[GaussPoints];
  dphideta = new double*[GaussPoints];
  dphidzeta= new double*[GaussPoints];

  phi_memory=new double [GaussPoints*nc_];
  dphidxi_memory  =new double [GaussPoints*nc_];
  dphideta_memory =new double [GaussPoints*nc_];
  dphidzeta_memory=new double [GaussPoints*nc_];

  for (unsigned i=0; i<GaussPoints; i++) {
    phi[i]=&phi_memory[i*nc_];
    dphidxi[i]  =&dphidxi_memory[i*nc_];
    dphideta[i] =&dphideta_memory[i*nc_];
    dphidzeta[i]=&dphidzeta_memory[i*nc_];
  }

  const double *ptx=GaussWeight+GaussPoints,*pty=GaussWeight+2*GaussPoints,*ptz=GaussWeight+3*GaussPoints;
  for (unsigned i=0; i<GaussPoints; i++,ptx++,pty++,ptz++) {
    double x[3];
    x[0]=*ptx;
    x[1]=*pty;
    x[2]=*ptz;
    for (int j=0; j<nc_; j++) {
      phi[i][j]=pt_basis->eval_phi(IND[j],x);
      dphidxi[i][j]=pt_basis->eval_dphidx(IND[j],x);
      dphideta[i][j]=pt_basis->eval_dphidy(IND[j],x);
      dphidzeta[i][j]=pt_basis->eval_dphidz(IND[j],x);
    }
  }

}


//----------------------------------------------------------------------------------------------------
// prolungator for LsysPde  Matrix 
//----------------------------------------------------------------------------------------------------
void elem_type::BuildProlongation(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat, 
				  const unsigned &index_sol, const unsigned &kkindex_sol) const {
  vector<int> cols(27);
   
  for (int i=0; i<nf_; i++) {
    int i0=KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=lspdec._msh->el->GetChildElement(ielc,i0);
    int i1=KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=lspdef._msh->el->GetDof(ielf,i1,type_);
    int irow=lspdef.GetKKDof(index_sol,kkindex_sol,iadd);  //  local-id to dof 
    int ncols=prol_ind[i+1]-prol_ind[i];
    cols.assign(ncols,0);
    for (int k=0; k<ncols; k++) {
      int j=prol_ind[i][k]; 
      int jadd=lspdec._msh->el->GetDof(ielc,j,type_);
      int jj=lspdec.GetKKDof(index_sol,kkindex_sol,jadd); 
      cols[k]=jj;
    }
    Projmat->insert_row(irow,ncols,cols,prol_val[i]);
  }
}



void elem_type::BuildRestrictionTranspose(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat, 
					  const unsigned &index_sol, const unsigned &kkindex_sol, const bool &TestDisp) const {
  vector<int> cols(27);
  bool fluid_region = (2==lspdec._msh->el->GetElementMaterial(ielc))?1:0;
  
  vector <double> copy_prol_val;
  copy_prol_val.reserve(27); 
  for (int i=0; i<nf_; i++) {
    int i0=KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=lspdec._msh->el->GetChildElement(ielc,i0);
    int i1=KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=lspdef._msh->el->GetDof(ielf,i1,type_);
        
    int irow=lspdef.GetKKDof(index_sol,kkindex_sol,iadd);  //  local-id to dof 
    int ncols=prol_ind[i+1]-prol_ind[i];
    
    bool isolidmark=lspdef._msh->el->GetNodeRegion(iadd);
    
    cols.assign(ncols,0);
    copy_prol_val.resize(ncols);
    for (int k=0; k<ncols; k++) {
      int j=prol_ind[i][k]; 
      int jadd=lspdec._msh->el->GetDof(ielc,j,type_);
      int jj=lspdec.GetKKDof(index_sol,kkindex_sol,jadd); 
      cols[k]=jj;
      
      bool jsolidmark=lspdef._msh->el->GetNodeRegion(jadd); 
      
      copy_prol_val[k]=(!TestDisp || !fluid_region || isolidmark==jsolidmark)?prol_val[i][k]:0.;
    }
      //Projmat->insert_row(irow,ncols,cols,prol_val[i]);
    Projmat->insert_row(irow,ncols,cols,&copy_prol_val[0]);
  }
}



//----------------------------------------------------------------------------------------------------
// prolungator for single solution
//-----------------------------------------------------------------------------------------------------
void elem_type::prolongation(const mesh &meshf,const mesh &meshc, const int& ielc,
			    SparseMatrix* Projmat) const {
  vector<int> cols(27);
  for (int i=0; i<nf_; i++) {
    int i0=KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=meshc.el->GetChildElement(ielc,i0);
    int i1=KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=meshf.el->GetDof(ielf,i1,type_);
    int irow=meshf.GetMetisDof(iadd,SolType_);  //  local-id to dof 
    int ncols=prol_ind[i+1]-prol_ind[i];
    cols.assign(ncols,0);
    for (int k=0; k<ncols; k++) {
      int j=prol_ind[i][k]; 
      int jadd=meshc.el->GetDof(ielc,j,type_);
      int jj=meshc.GetMetisDof(jadd,SolType_);
      cols[k]=jj;
    }
    Projmat->insert_row(irow,ncols,cols,prol_val[i]);
  }
}

//----------------------------------------------------------------------------------------------------
// prolungator for solution printing
//-----------------------------------------------------------------------------------------------------
void elem_type::ProlQitoQj(const mesh& mymesh,const int& iel, SparseMatrix* Projmat, 
			  bool testnode[],const unsigned &itype) const{
  vector<int> cols(27);
  for (int i=0; i<ncf_[itype]; i++) {
    int inode=mymesh.el->GetDof(iel,i,type_);
    int irow=mymesh.GetMetisDof(inode,itype);
    if (testnode[irow]==0) {
      testnode[irow]=1;
      int ncols=prol_ind[i+1]-prol_ind[i];
      cols.assign(ncols,0);
      for (int k=0; k<ncols; k++) {
        int jj=prol_ind[i][k];
        int jnode=mymesh.el->GetDof(iel,jj,type_);
	int jadd=mymesh.GetMetisDof(jnode,SolType_);
        cols[k]=jadd;
      }
      Projmat->insert_row(irow,ncols,cols,prol_val[i]);
    }
  }
}

//---------------------------------------------------------------------------------------------------------
void elem_type::JacobianSur2D(const double vt[][27],const unsigned &ig,
                              double &Weight, double *other_phi, double gradphi[][3], double normal[3]) const {

    double Jac[3][3];
    double JacI[3][3];

    Jac[0][0] = 0.;
    Jac[1][0] = 0.;
    Jac[2][0] = 0.;
    
    Jac[0][1] = 0.;
    Jac[1][1] = 0.;
    Jac[2][1] = 0.;
 
//  for(double *pt_d=Jac[0]; pt_d<Jac[0]+9; pt_d++) *pt_d=0.;
    const double *dfx=dphidxi[ig];
    const double *dfy=dphideta[ig];

    const double *vx=vt[0];
    const double *vy=vt[1];
    const double *vz=vt[2];
  
   for(int inode=0;inode<nc_;inode++,dfx++,dfy++,vx++,vy++,vz++){
       
       Jac[0][0] += (*dfx)*(*vx);
       Jac[1][0] += (*dfx)*(*vy);
       Jac[2][0] += (*dfx)*(*vz);
    
       Jac[0][1] += (*dfy)*(*vx);
       Jac[1][1] += (*dfy)*(*vy);
       Jac[2][1] += (*dfy)*(*vz);
    }
    
    //   normal module
    double nx = Jac[1][0]*Jac[2][1] - Jac[1][1]*Jac[2][0];
    double ny = Jac[0][1]*Jac[2][0] - Jac[2][1]*Jac[0][0]; 
    double nz = Jac[0][0]*Jac[1][1] - Jac[0][1]*Jac[1][0];
    double modn = sqrt(nx*nx + ny*ny + nz*nz);

    normal[0] =  (nx)/modn;
    normal[1] =  (ny)/modn;
    normal[2] =  (nz)/modn;
    
    Jac[0][2] += normal[0];
    Jac[1][2] += normal[1];
    Jac[2][2] += normal[2];
    
    
   //the determinant of the matrix is the area 
   double det=(Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])+
               Jac[0][1]*(Jac[1][2]*Jac[2][0]-Jac[1][0]*Jac[2][2])+
               Jac[0][2]*(Jac[1][0]*Jac[2][1]-Jac[1][1]*Jac[2][0]));

   Weight=det*GaussWeight[ig];

   double *fi=phi[ig];

   for(int inode=0;inode<nc_;inode++,other_phi++,fi++){
     *other_phi=*fi;
   }
}

//---------------------------------------------------------------------------------------------------------
void elem_type::JacobianSur1D(const double vt[][27],const unsigned &ig,
                              double &Weight, double *other_phi, double gradphi[][3], double normal[3]) const {

//      cout << "Calling 1d surface jacobian" << endl;
  double Jac[2][2];
  double JacI[2][2];

  Jac[0][0] = 0.;
  Jac[1][0] = 0.;

  const double *dfeta=dphidxi[ig];
  const double *vx=vt[0];
  const double *vy=vt[1];

  for (int inode=0; inode<nc_; inode++,dfeta++,vx++,vy++) {
    Jac[0][0] += (*dfeta)*(*vx);
    Jac[1][0] += (*dfeta)*(*vy);
  }

//   normal module
  double modn = sqrt(Jac[0][0]*Jac[0][0] + Jac[1][0]*Jac[1][0]);

  normal[0] =  Jac[1][0]/modn;
  normal[1] = -Jac[0][0]/modn;
  normal[2] =  0.;

  //The derivative of x with respect to eta (dx/deta) has the opposite sign with respect to the normal
  //obtained as cross product between (dx/deta , dy/deta, 0) x (0,0,1)
  //The Jacobian has the structure
  // |dx/deta  -nx|
  // |dy/deta  -ny|
  Jac[0][1] = -normal[0];
  Jac[1][1] = -normal[1];

  //The determinant of that matrix is the area
  double det= (Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]);

  JacI[0][0] =  Jac[1][1]/det;
  JacI[0][1] = -Jac[0][1]/det;
  JacI[1][0] = -Jac[1][0]/det;
  JacI[1][1] =  Jac[0][0]/det;

  Weight = det*GaussWeight[ig];

  double *gradf=gradphi[0];
  double *fi=phi[ig];
  dfeta=dphidxi[ig];
  for (int inode=0; inode<nc_; inode++,other_phi++,fi++,dfeta++) {
    *other_phi=*fi;
    //the derivative computed in this way are wrong
    *(gradf++)=(*dfeta)*JacI[0][0];
    *(gradf++)=(*dfeta)*JacI[0][1];
    *(gradf++)=0.;

  }

}

//---------------------------------------------------------------------------------------------------------
void elem_type::Jacobian3D(const vector < vector < double > > &vt,const unsigned &ig,
                           double &Weight, vector < double > &other_phi, vector < double > &gradphi) const {
  double Jac[3][3];
  double JacI[3][3];
  for (double *pt_d=Jac[0]; pt_d<Jac[0]+9; pt_d++) *pt_d=0.;
  const double *dfx=dphidxi[ig];
  const double *dfy=dphideta[ig];
  const double *dfz=dphidzeta[ig];
  const double *vx=&vt[0][0];
  const double *vy=&vt[1][0];
  const double *vz=&vt[2][0];
  for (int inode=0; inode<nc_; inode++,dfx++,dfy++,dfz++,vx++,vy++,vz++) {
    double *pt_d=Jac[0];
    *(pt_d++)+=(*dfx)*(*vx);
    *(pt_d++)+=(*dfx)*(*vy);
    *(pt_d++)+=(*dfx)*(*vz);
    *(pt_d++)+=(*dfy)*(*vx);
    *(pt_d++)+=(*dfy)*(*vy);
    *(pt_d++)+=(*dfy)*(*vz);
    *(pt_d++)+=(*dfz)*(*vx);
    *(pt_d++)+=(*dfz)*(*vy);
    *(pt_d++)+=(*dfz)*(*vz);
  }
  double det=(Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])+
              Jac[0][1]*(Jac[1][2]*Jac[2][0]-Jac[1][0]*Jac[2][2])+
              Jac[0][2]*(Jac[1][0]*Jac[2][1]-Jac[1][1]*Jac[2][0]));

  JacI[0][0]= (-Jac[1][2]*Jac[2][1] + Jac[1][1]*Jac[2][2])/det;
  JacI[0][1]= ( Jac[0][2]*Jac[2][1] - Jac[0][1]*Jac[2][2])/det;
  JacI[0][2]= (-Jac[0][2]*Jac[1][1] + Jac[0][1]*Jac[1][2])/det;
  JacI[1][0]= ( Jac[1][2]*Jac[2][0] - Jac[1][0]*Jac[2][2])/det;
  JacI[1][1]= (-Jac[0][2]*Jac[2][0] + Jac[0][0]*Jac[2][2])/det;
  JacI[1][2]= ( Jac[0][2]*Jac[1][0] - Jac[0][0]*Jac[1][2])/det;
  JacI[2][0]= (-Jac[1][1]*Jac[2][0] + Jac[1][0]*Jac[2][1])/det;
  JacI[2][1]= ( Jac[0][1]*Jac[2][0] - Jac[0][0]*Jac[2][1])/det;
  JacI[2][2]= (-Jac[0][1]*Jac[1][0] + Jac[0][0]*Jac[1][1])/det;

  Weight=det*GaussWeight[ig];

  double *other_f=&other_phi[0];
  double *gradf=&gradphi[0];
  double *fi=phi[ig];
  dfx=dphidxi[ig];
  dfy=dphideta[ig];
  dfz=dphidzeta[ig];
  for (int inode=0; inode<nc_; inode++,other_f++,fi++,dfx++,dfy++,dfz++) {
    *other_f=*fi;
    double *pt_d=JacI[0];
    for (int j=0; j<3; j++) {
      *(gradf++)=(*dfx)*(*pt_d)+(*dfy)*(*(pt_d+1))+(*dfz)*(*(pt_d+2));
      pt_d+=3;
    }
  }
}

//---------------------------------------------------------------------------------------------------------
void elem_type::Jacobian2D(const vector < vector < double > > &vt,const unsigned &ig,
                           double &Weight, vector < double > &other_phi, vector < double > &gradphi) const {

  double Jac[2][2];
  double JacI[2][2];
  for (double *pt_d=Jac[0]; pt_d<Jac[0]+4; pt_d++) *pt_d=0.;
  const double *dfx=dphidxi[ig];
  const double *dfy=dphideta[ig];
  const double *vx=&vt[0][0];
  const double *vy=&vt[1][0];
  for (int inode=0; inode<nc_; inode++,dfx++,dfy++,vx++,vy++) {
    double *pt_d=Jac[0];
    *(pt_d++)+=(*dfx)*(*vx);
    *(pt_d++)+=(*dfx)*(*vy);
    *(pt_d++)+=(*dfy)*(*vx);
    *(pt_d++)+=(*dfy)*(*vy);
  }
  double det=(Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]);

  JacI[0][0]= Jac[1][1]/det;
  JacI[0][1]=-Jac[0][1]/det;
  JacI[1][0]=-Jac[1][0]/det;
  JacI[1][1]= Jac[0][0]/det;

  Weight=det*GaussWeight[ig];

  double *other_f=&other_phi[0];
  double *gradf=&gradphi[0];
  double *fi=phi[ig];
  dfx=dphidxi[ig];
  dfy=dphideta[ig];
  for (int inode=0; inode<nc_; inode++,other_f++,fi++,dfx++,dfy++) {
    *other_f=*fi;
    double *pt_d=JacI[0];
    for (int j=0; j<2; j++) {
      *(gradf++)=(*dfx)*(*pt_d)+(*dfy)*(*(pt_d+1));
      pt_d+=2;
    }
    //*(gradf++)=0.;
  }
}

//---------------------------------------------------------------------------------------------------------
void elem_type::Jacobian1D(const vector < vector < double > > &vt,const unsigned &ig,
                           double &Weight, vector < double > &other_phi, vector < double > &gradphi) const {

  double Jac=0.;
  double h;

  const double *dfx=dphidxi[ig];
  const double *vx=&vt[0][0];
//   Jac = fabs((*(vx++)) - (*vx));
  
   for (int inode=0; inode<nc_; inode++,dfx++,vx++) {
     Jac+=(*dfx)*(*vx);
   }
  Weight=Jac*GaussWeight[ig];

  double *other_f=&other_phi[0];
  double *gradf=&gradphi[0];
  double *fi=phi[ig];
  dfx=dphidxi[ig];
  for (int inode=0; inode<nc_; inode++,other_f++,fi++,dfx++) {
    *other_f=*fi;
    *(gradf++)=(*dfx)*(1./Jac);
  }
}

//---------------------------------------------------------------------------------------------------------
double* elem_type::GetPhi(const unsigned &ig) const {
  return phi[ig];
}

double* elem_type::GetDPhiDXi(const unsigned &ig) const {
  return dphidxi[ig];
}

double* elem_type::GetDPhiDEta(const unsigned &ig) const {
  return dphideta[ig];
}

double* elem_type::GetDPhiDZeta(const unsigned &ig) const {
  return dphidzeta[ig];
}


//---------------------------------------------------------------------------------------------------------
void elem_type::GetArea(const double *vtx,const double *vty, const double *vtz, const unsigned &ig,
                        double &Weight, double *other_phi) const {

  double Jac[2][3];
  for (double *pt_d=Jac[0]; pt_d<Jac[0]+6; pt_d++) *pt_d=0.;
  const double *dfx=dphidxi[ig];
  const double *dfy=dphideta[ig];
  const double *vx=vtx;
  const double *vy=vty;
  const double *vz=vtz;
  for (int inode=0; inode<nc_; inode++,dfx++,dfy++,vx++,vy++,vz++) {
    double *pt_d=Jac[0];
    *(pt_d++)+=(*dfx)*(*vx);
    *(pt_d++)+=(*dfx)*(*vy);
    *(pt_d++)+=(*dfx)*(*vz);
    *(pt_d++)+=(*dfy)*(*vx);
    *(pt_d++)+=(*dfy)*(*vy);
    *(pt_d++)+=(*dfy)*(*vz);
  }

  double det1=Jac[0][1]*Jac[1][2]-Jac[1][1]*Jac[0][2];
  double det2=Jac[0][0]*Jac[1][2]-Jac[1][0]*Jac[0][2];
  double det3=Jac[0][0]*Jac[1][1]-Jac[1][0]*Jac[0][1];
  double det=sqrt(det1*det1+det2*det2+det3*det3);

  Weight=det*GaussWeight[ig];
  double *fi=phi[ig];
  for (int inode=0; inode<nc_; inode++,other_phi++,fi++) {
    *other_phi=*fi;
  }
}


} //end namespace femus


