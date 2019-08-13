/*=========================================================================

 Program: FEMuS
 Module: Marker
 Authors: Eugenio Aulisa and Giacomo Capodaglio

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <stdlib.h>
#include<iostream>
#include<cmath>

#include "PolynomialBases.hpp"
#include "GeomElTypeEnum.hpp"

namespace femus {

//BEGIN Interface
  void ProjectNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector< std::vector < double > > &aN,
      const short unsigned &ielType, const unsigned &solType) {

    if(ielType == QUAD) {
      ProjectQuadNodalToPolynomialCoefficients(aP, aN, solType);
    }

    else if(ielType == TRI) {
      ProjectTriNodalToPolynomialCoefficients(aP, aN, solType);
    }

    else if(ielType == HEX) {
      ProjectHexNodalToPolynomialCoefficients(aP, aN, solType);
    }

    else if(ielType == TET) {
      ProjectTetNodalToPolynomialCoefficients(aP, aN, solType);
    }

    else if(ielType == WEDGE) {
      ProjectWedgeNodalToPolynomialCoefficients(aP, aN, solType);
    }
  }

  void InterpolatePolynomialCoefficients(std::vector<std::vector < std::vector <double > > > &aXs, const std::vector<std::vector < std::vector <double > > > &aX0,
                                         const std::vector<std::vector < std::vector <double > > > aX1, const double &s){
     
    aXs.resize(aX0.size());
    for(unsigned i=0; i<aX0.size(); i++){
      aXs[i].resize(aX0[i].size());
      for(unsigned j=0; j<aX0[i].size(); j++){
	aXs[i][j].resize(aX0[i][j].size());
	for(unsigned k=0; k<aX0[i][j].size(); k++){
	  aXs[i][j][k] = (1.-s)*aX0[i][j][k] + s * aX1[i][j][k];
	}
      }
    }
  }


  void GetPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi,
                                  short unsigned &ielType, const unsigned & solType) {

    if(ielType == QUAD) {
      GetQuadPolynomialShapeFunction(phi,  xi,  solType);
    }

    else if(ielType == TRI) {
      GetTriPolynomialShapeFunction(phi, xi, solType);
    }

    else if(ielType == HEX) {
      GetHexPolynomialShapeFunction(phi, xi, solType);
    }

    else if(ielType == TET) {
      GetTetPolynomialShapeFunction(phi, xi, solType);
    }

    else if(ielType == WEDGE) {
      GetWedgePolynomialShapeFunction(phi, xi, solType);
    }
  }

  void GetPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
                                          const std::vector < double >& xi,  short unsigned &ielType, const unsigned & solType) {

    if(ielType == QUAD) {
      GetQuadPolynomialShapeFunctionGradient(phi, gradPhi,  xi,  solType);
    }

    else if(ielType == TRI) {
      GetTriPolynomialShapeFunctionGradient(phi, gradPhi, xi, solType);
    }

    else if(ielType == HEX) {
      GetHexPolynomialShapeFunctionGradient(phi, gradPhi, xi, solType);
    }

    else if(ielType == TET) {
      GetTetPolynomialShapeFunctionGradient(phi, gradPhi, xi, solType);
    }

    else if(ielType == WEDGE) {
      GetWedgePolynomialShapeFunctionGradient(phi, gradPhi, xi, solType);
    }
  }

  void GetPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi,
      short unsigned &ielType, const unsigned & solType) {

    if(ielType == QUAD) {
      GetQuadPolynomialShapeFunctionGradientHessian(phi, gradPhi,  hessPhi,  xi,  solType);
    }

    else if(ielType == TRI) {
      GetTriPolynomialShapeFunctionGradientHessian(phi, gradPhi, hessPhi, xi, solType);
    }

    else if(ielType == HEX) {
      GetHexPolynomialShapeFunctionGradientHessian(phi, gradPhi, hessPhi, xi, solType);
    }

    else if(ielType == TET) {
      GetTetPolynomialShapeFunctionGradientHessian(phi, gradPhi, hessPhi, xi, solType);
    }

    else if(ielType == WEDGE) {
      GetWedgePolynomialShapeFunctionGradientHessian(phi, gradPhi, hessPhi, xi, solType);
    }
  }
//END Interface

//BEGIN QUAD specialized functions
  const unsigned quadNumberOfDofs[3] = {4, 8, 9};

  void ProjectQuadNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solType) {

    unsigned dim =  aN.size();
    aP.resize(dim);

    unsigned nDofs = quadNumberOfDofs[solType];
    if(nDofs > aN[0].size()) {
      std::cout << "Error in ProjectQuadNodalToPolynomialCoefficients(...) the number of Dofs is inconsistent" << std::endl;
      abort();
    }

    for(unsigned k = 0; k < dim; k++) {
      aP[k].resize(nDofs);
    }

    if(solType == 0) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = 0.25 * (aN[k][0] + aN[k][1] + aN[k][2] + aN[k][3]);
        aP[k][1] = 0.25 * (-aN[k][0] + aN[k][1] + aN[k][2] - aN[k][3]);
        aP[k][2] = 0.25 * (-aN[k][0] - aN[k][1] + aN[k][2] + aN[k][3]);
        aP[k][3] = 0.25 * (aN[k][0] - aN[k][1] + aN[k][2] - aN[k][3]);
      }
    }
    else if(solType == 1) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = -0.25 * (aN[k][0] + aN[k][1] + aN[k][2] + aN[k][3]) +
                   0.5 * (aN[k][4] + aN[k][5] + aN[k][6] + aN[k][7]);
        aP[k][1] = 0.5 * (aN[k][5] - aN[k][7]);
        aP[k][2] = 0.5 * (aN[k][6] - aN[k][4]);
        aP[k][3] = 0.25 * (aN[k][0] - aN[k][1] + aN[k][2] - aN[k][3]);
        aP[k][4] = 0.25 * (aN[k][0] + aN[k][1] + aN[k][2] + aN[k][3]) - 0.5 * (aN[k][4] + aN[k][6]);
        aP[k][5] = 0.25 * (aN[k][0] + aN[k][1] + aN[k][2] + aN[k][3]) - 0.5 * (aN[k][5] + aN[k][7]);
        aP[k][6] = 0.25 * (-aN[k][0] - aN[k][1] + aN[k][2] + aN[k][3]) + 0.5 * (aN[k][4] - aN[k][6]);
        aP[k][7] = 0.25 * (-aN[k][0] + aN[k][1] + aN[k][2] - aN[k][3]) + 0.5 * (aN[k][7] - aN[k][5]);
      }
    }
    else if(solType == 2) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][8];
        aP[k][1] = 0.5 * (aN[k][5] - aN[k][7]);
        aP[k][2] = 0.5 * (aN[k][6] - aN[k][4]);
        aP[k][3] = 0.25 * (aN[k][0] - aN[k][1] + aN[k][2] - aN[k][3]);
        aP[k][4] =  0.5 * (aN[k][5] + aN[k][7]) - aN[k][8];
        aP[k][5] =  0.5 * (aN[k][4] + aN[k][6]) - aN[k][8];
        aP[k][6] = 0.25 * (-aN[k][0] - aN[k][1] + aN[k][2] + aN[k][3]) +
                   0.5 * (aN[k][4] - aN[k][6]);
        aP[k][7] = 0.25 * (-aN[k][0] + aN[k][1] + aN[k][2] - aN[k][3]) +
                   0.5 * (-aN[k][5] + aN[k][7]);
        aP[k][8] = 0.25 * (aN[k][0] + aN[k][1] + aN[k][2] + aN[k][3]) -
                   0.5 * (aN[k][4] + aN[k][5] + aN[k][6] + aN[k][7]) + aN[k][8];
      }
    }
  }


  void GetQuadPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) {

    const unsigned nDofs = quadNumberOfDofs[solType];

    phi.resize(nDofs);

    phi[0] = 1.;
    phi[1] = xi[0]; // x
    phi[2] = xi[1];  // y
    phi[3] = xi[0] * xi[1];  // x y

    if(solType > 0) {
      phi[4] = xi[0] * xi[0]; // x x
      phi[5] = xi[1] * xi[1]; // y y
      phi[6] = phi[4] * xi[1]; // xx y
      phi[7] = phi[5] * xi[0]; // x yy

      if(solType > 1) {
        phi[8] = phi[3] * phi[3]; // xx yy
      }
    }

  }

  void GetQuadPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      const std::vector < double >& xi, const unsigned & solType) {

    GetQuadPolynomialShapeFunction(phi,  xi, solType);

    const unsigned dim = 2;

    const unsigned nDofs = quadNumberOfDofs[solType];

    gradPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      gradPhi[i].assign(dim, 0.);
    }

    //phi_x
    gradPhi[1][0] = 1.; // 1
    gradPhi[3][0] = xi[1];  //  y
    //phi_y
    gradPhi[2][1] = 1.;  // 1
    gradPhi[3][1] = xi[0];  // x

    if(solType > 0) {
      //phi_x
      gradPhi[4][0] = 2. * xi[0]; // 2 x
      gradPhi[6][0] = 2. * phi[3]; // 2 x y
      gradPhi[7][0] = phi[5]; //  yy
      //phi_y
      gradPhi[5][1] = 2. * xi[1]; // 2 y
      gradPhi[6][1] = phi[4]; // xx
      gradPhi[7][1] = 2. * phi[3]; // 2 x y

      if(solType > 1) {
        //phi_x
        gradPhi[8][0] = 2. * phi[7]; // 2 x yy
        //phi_y
        gradPhi[8][1] = 2. * phi[6]; // 2 xx y
      }
    }
  }

  void GetQuadPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) {

    GetQuadPolynomialShapeFunctionGradient(phi,   gradPhi,  xi, solType);

    const unsigned dim = 2;
    const unsigned nDofs = quadNumberOfDofs[solType];

    hessPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      hessPhi[i].resize(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessPhi[i][i1].assign(dim, 0.);
      }
    }

    //phi_xy
    hessPhi[3][1][0] = hessPhi[3][0][1] = 1.; //1

    if(solType > 0) {
      //phi_xx
      hessPhi[4][0][0] = 2.; // 2
      hessPhi[6][0][0] = 2. * xi[1]; // 2 y
      //phi_xy
      hessPhi[6][1][0] = hessPhi[6][0][1] = 2. * xi[0]; // 2 x
      hessPhi[7][1][0] = hessPhi[7][0][1] = 2. * xi[1]; // 2 y
      //phi_yy
      hessPhi[5][1][1] = 2.; // 2
      hessPhi[7][1][1] = 2. * xi[0]; // 2 x

      if(solType > 1) {
        //phi_xx
        hessPhi[8][0][0] = 2. * phi[5]; // 2 yy
        //phi_xy
        hessPhi[8][1][0] = hessPhi[8][0][1] = 4. * phi[3]; // 4 x y
        //phi_yy
        hessPhi[8][1][1] = 2. * phi[4]; // 2 xx
      }
    }
  }
//END QUAD

//BEGIN TRI specialized functions
  const unsigned triNumberOfDofs[3] = {3, 6, 7};

  void ProjectTriNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solType) {

    unsigned dim =  aN.size();
    aP.resize(dim);

    unsigned nDofs = triNumberOfDofs[solType];
    if(nDofs > aN[0].size()) {
      std::cout << "Error in ProjectTriNodalToPolynomialCoefficients(...) the number of Dofs is inconsistent" << std::endl;
      abort();
    }

    for(unsigned k = 0; k < dim; k++) {
      aP[k].resize(nDofs);
    }

    if(solType == 0) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][0];
        aP[k][1] = - aN[k][0] + aN[k][1];
        aP[k][2] = - aN[k][0] + aN[k][2];
      }
    }
    else if(solType == 1) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][0];
        aP[k][1] = - 3 * aN[k][0] - aN[k][1] + 4 * aN[k][3];
        aP[k][2] = - 3 * aN[k][0] - aN[k][2] + 4 * aN[k][5];
        aP[k][3] = 4 * aN[k][0] - 4 * aN[k][3] + 4 * aN[k][4] - 4 * aN[k][5];
        aP[k][4] = 2 * aN[k][0] + 2 * aN[k][1] - 4 * aN[k][3];
        aP[k][5] = 2 * aN[k][0] + 2 * aN[k][2] - 4 * aN[k][5];
      }
    }
    else if(solType == 2) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][0];
        aP[k][1] = - 3 * aN[k][0] - aN[k][1] + 4 * aN[k][3];
        aP[k][2] = - 3 * aN[k][0] - aN[k][2] + 4 * aN[k][5];
        aP[k][3] = 7 * aN[k][0] + 3 * aN[k][1] + 3 * aN[k][2] - 16 * aN[k][3] -
                   8 * aN[k][4] - 16 * aN[k][5] + 27 * aN[k][6];
        aP[k][4] = 2 * aN[k][0] + 2 * aN[k][1] - 4 * aN[k][3];
        aP[k][5] = 2 * aN[k][0] + 2 * aN[k][2] - 4 * aN[k][5];
        aP[k][6] = - 3 * aN[k][0] - 3 * aN[k][1] - 3 * aN[k][2] + 12 * aN[k][3] +
                   12 * aN[k][4] + 12 * aN[k][5] - 27 * aN[k][6];
      }
    }
  }
  

  void GetTriPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) {

    const unsigned nDofs = triNumberOfDofs[solType];

    phi.resize(nDofs);

    phi[0] = 1.;
    phi[1] = xi[0]; // x
    phi[2] = xi[1]; // y

    if(solType > 0) {
      phi[3] = xi[0] * xi[1]; // x y
      phi[4] = xi[0] * xi[0];  // x x
      phi[5] = xi[1] * xi[1]; // y y

      if(solType > 1) {
        phi[6] = phi[4] * xi[1] + phi[5] * xi[0]; // xx y + x yy
      }
    }

  }

  void GetTriPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      const std::vector < double >& xi, const unsigned & solType) {

    GetTriPolynomialShapeFunction(phi,  xi, solType);

    const unsigned dim = 2;

    const unsigned nDofs = triNumberOfDofs[solType];

    gradPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      gradPhi[i].assign(dim, 0.);
    }

    //phi_x
    gradPhi[1][0] = 1.; // 1
    //phi_y
    gradPhi[2][1] = 1.;  // 1

    if(solType > 0) {
      //phi_x
      gradPhi[3][0] = xi[1]; // y
      gradPhi[4][0] = 2.*xi[0] ;  // 2 x
      //phi_y
      gradPhi[3][1] = xi[0]; // x
      gradPhi[5][1] = 2.*xi[1]; // 2*y


      if(solType > 1) {
        //phi_x
        gradPhi[6][0] = gradPhi[4][0] * xi[1] + phi[5]; // 2 x y + y y
        //phi_y
        gradPhi[6][1] = phi[4] + gradPhi[5][1] * xi[0]; // xx  + 2 x y
      }
    }
  }

  void GetTriPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) {

    GetTriPolynomialShapeFunctionGradient(phi,   gradPhi,  xi, solType);

    const unsigned dim = 2;
    const unsigned nDofs = triNumberOfDofs[solType];

    hessPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      hessPhi[i].resize(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessPhi[i][i1].assign(dim, 0.);
      }
    }

    if(solType > 0) {
      //phi_xx
      hessPhi[4][0][0] = 2.; // 2
      //phi_xy
      hessPhi[3][1][0] = hessPhi[3][0][1] = 1.; // 1
      //phi_yy
      hessPhi[5][1][1] = 2.; // 2

      if(solType > 1) {
        //phi_xx
        hessPhi[6][0][0] = 2. * xi[1]; // 2 y
        //phi_xy
        hessPhi[6][1][0] = hessPhi[6][0][1] = 2. * xi[0] + 2. * xi[1]; // 2. x + 2 * y
        //phi_yy
        hessPhi[6][1][1] = 2. * xi[0]; // 2 * x
      }
    }
  }
//END TRI

//BEGIN HEX specialized functions
  const unsigned hexNumberOfDofs[3] = {8, 20, 27};

  void ProjectHexNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solType) {

    unsigned dim =  aN.size();
    aP.resize(dim);

    unsigned nDofs = hexNumberOfDofs[solType];
    if(nDofs > aN[0].size()) {
      std::cout << "Error in ProjectHexNodalToPolynomialCoefficients(...) the number of Dofs is inconsistent" << std::endl;
      abort();
    }

    for(unsigned k = 0; k < dim; k++) {
      aP[k].resize(nDofs);
    }

    if(solType == 0) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = 0.125 * (aN[k][0] + aN[k][1] + aN[k][2] + aN[k][3] + aN[k][4] +
                            aN[k][5] + aN[k][6] + aN[k][7]) ;
        aP[k][1] = 0.125 * (- aN[k][0] + aN[k][1] + aN[k][2] - aN[k][3] - aN[k][4] +
                            aN[k][5] + aN[k][6] - aN[k][7]) ;
        aP[k][2] = 0.125 * (- aN[k][0] - aN[k][1] + aN[k][2] + aN[k][3] - aN[k][4] -
                            aN[k][5] + aN[k][6] + aN[k][7]) ;
        aP[k][3] = 0.125 * (- aN[k][0] - aN[k][1] - aN[k][2] - aN[k][3] + aN[k][4] +
                            aN[k][5] + aN[k][6] + aN[k][7]) ;
        aP[k][4] = 0.125 * (aN[k][0] - aN[k][1] + aN[k][2] - aN[k][3] + aN[k][4] -
                            aN[k][5] + aN[k][6] - aN[k][7]) ;
        aP[k][5] = 0.125 * (aN[k][0] - aN[k][1] - aN[k][2] + aN[k][3] - aN[k][4] +
                            aN[k][5] + aN[k][6] - aN[k][7]) ;
        aP[k][6] = 0.125 * (aN[k][0] + aN[k][1] - aN[k][2] - aN[k][3] - aN[k][4] -
                            aN[k][5] + aN[k][6] + aN[k][7]) ;
        aP[k][7] = 0.125 * (- aN[k][0] + aN[k][1] - aN[k][2] + aN[k][3] + aN[k][4] -
                            aN[k][5] + aN[k][6] - aN[k][7]) ;
      }
    }
    else if(solType == 1) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = 0.25 * (- aN[k][0] + aN[k][9] + aN[k][10] + aN[k][11] + aN[k][12] + aN[k][13] + aN[k][14] +
                           aN[k][15] + aN[k][16] + aN[k][17] + aN[k][18] - aN[k][1] + aN[k][19] - aN[k][2] -
                           aN[k][3] - aN[k][4] - aN[k][5] - aN[k][6] - aN[k][7] + aN[k][8]);
        aP[k][1] = 0.125 * (aN[k][0] + 2 * aN[k][9] - 2 * aN[k][11] + 2 * aN[k][13] - 2 * aN[k][15] - 2 * aN[k][16] +
                            2 * aN[k][17] + 2 * aN[k][18] - aN[k][1] - 2 * aN[k][19] - aN[k][2] + aN[k][3] +
                            aN[k][4] - aN[k][5] - aN[k][6] + aN[k][7]);
        aP[k][2] = 0.125 * (aN[k][0] + 2 * aN[k][10] - 2 * aN[k][12] + 2 * aN[k][14] - 2 * aN[k][16] - 2 * aN[k][17] +
                            2 * aN[k][18] + aN[k][1] + 2 * aN[k][19] - aN[k][2] - aN[k][3] + aN[k][4] + aN[k][5] - aN[k][6] - aN[k][7] - 2 * aN[k][8]);
        aP[k][3] = 0.125 * (aN[k][0] - 2 * aN[k][9] - 2 * aN[k][10] - 2 * aN[k][11] + 2 * aN[k][12] + 2 * aN[k][13] + 2 * aN[k][14] +
                            2 * aN[k][15] + aN[k][1] + aN[k][2] + aN[k][3] - aN[k][4] - aN[k][5] - aN[k][6] - aN[k][7] - 2 * aN[k][8]);
        aP[k][4] = 0.25 * (aN[k][16] - aN[k][17] + aN[k][18] - aN[k][19]);
        aP[k][5] = 0.25 * (- aN[k][9] + aN[k][11] + aN[k][13] - aN[k][15]);
        aP[k][6] = 0.25 * (- aN[k][10] - aN[k][12] + aN[k][14] + aN[k][8]);
        aP[k][7] = 0.125 * (aN[k][0] - 2 * aN[k][10] - 2 * aN[k][12] - 2 * aN[k][14] + aN[k][1] + aN[k][2] + aN[k][3] + aN[k][4] +
                            aN[k][5] + aN[k][6] + aN[k][7] - 2 * aN[k][8]);
        aP[k][8] = 0.125 * (aN[k][0] - 2 * aN[k][9] - 2 * aN[k][11] - 2 * aN[k][13] - 2 * aN[k][15] + aN[k][1] + aN[k][2] + aN[k][3] +
                            aN[k][4] + aN[k][5] + aN[k][6] + aN[k][7]);
        aP[k][9] = 0.125 * (aN[k][0] - 2  * aN[k][16] - 2 * aN[k][17] - 2 *  aN[k][18] + aN[k][1] - 2 * aN[k][19] + aN[k][2] + aN[k][3] +
                            aN[k][4] + aN[k][5] + aN[k][6] + aN[k][7]);
        aP[k][10] = 0.125 * (- aN[k][0] + aN[k][1] - aN[k][2] + aN[k][3] + aN[k][4] - aN[k][5] + aN[k][6] - aN[k][7]);
        aP[k][11] = 0.125 * (- aN[k][0] - 2 * aN[k][10] + 2 *  aN[k][12] - 2 * aN[k][14] - aN[k][1] + aN[k][2] + aN[k][3] - aN[k][4] -
                             aN[k][5] + aN[k][6] + aN[k][7] + 2 * aN[k][8]);
        aP[k][12] = 0.125 * (- aN[k][0] + 2 * aN[k][10] - 2 * aN[k][12] - 2 * aN[k][14] - aN[k][1] - aN[k][2] - aN[k][3] + aN[k][4] +
                             aN[k][5] + aN[k][6] + aN[k][7] + 2 * aN[k][8]);
        aP[k][13] = 0.125 * (- aN[k][0] + 2 * aN[k][9] + 2 * aN[k][11] - 2 * aN[k][13] - 2 * aN[k][15] - aN[k][1] - aN[k][2] - aN[k][3] +
                             aN[k][4] + aN[k][5] + aN[k][6] + aN[k][7]);
        aP[k][14] = 0.125 * (- aN[k][0] - 2 * aN[k][9] + 2 * aN[k][11] - 2 * aN[k][13] + 2 * aN[k][15] + aN[k][1] + aN[k][2] - aN[k][3] -
                             aN[k][4] + aN[k][5] + aN[k][6] - aN[k][7]);
        aP[k][15] = 0.125 * (- aN[k][0] + 2 * aN[k][16] - 2 * aN[k][17] - 2 * aN[k][18] + aN[k][1] + 2 * aN[k][19] + aN[k][2] - aN[k][3] - aN[k][4] +
                             aN[k][5] + aN[k][6] - aN[k][7]);
        aP[k][16] = 0.125 * (- aN[k][0] + 2 * aN[k][16] + 2 * aN[k][17] - 2 * aN[k][18] - aN[k][1] - 2 * aN[k][19] + aN[k][2] + aN[k][3] - aN[k][4] -
                             aN[k][5] + aN[k][6] + aN[k][7]);
        aP[k][17] = 0.125 * (aN[k][0] + 2 * aN[k][10] + 2 * aN[k][12] - 2 * aN[k][14] + aN[k][1] - aN[k][2] - aN[k][3] - aN[k][4] - aN[k][5] + aN[k][6] +
                             aN[k][7] - 2 * aN[k][8]);
        aP[k][18] = 0.125 * (aN[k][0] + 2 * aN[k][9] - 2 * aN[k][11] - 2 * aN[k][13] + 2 * aN[k][15] - aN[k][1] - aN[k][2] + aN[k][3] - aN[k][4] +
                             aN[k][5] + aN[k][6] - aN[k][7]);
        aP[k][19] = 0.125 * (aN[k][0] - 2 * aN[k][16] + 2 * aN[k][17] - 2 * aN[k][18] - aN[k][1] + 2 * aN[k][19] + aN[k][2] - aN[k][3] + aN[k][4] -
                             aN[k][5] + aN[k][6] - aN[k][7]);
      }
    }
    else if(solType == 2) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][26];
        aP[k][1] = 0.5 * (aN[k][21] - aN[k][23]);
        aP[k][2] = 0.5 * (aN[k][22] - aN[k][20]);
        aP[k][3] = 0.5 * (aN[k][25] - aN[k][24]);
        aP[k][4] = 0.25 * (aN[k][16] - aN[k][17] + aN[k][18] - aN[k][19]);
        aP[k][5] = 0.25 * (- aN[k][9] + aN[k][11] + aN[k][13] - aN[k][15]);
        aP[k][6] = 0.25 * (- aN[k][10] - aN[k][12] + aN[k][14] + aN[k][8]);
        aP[k][7] = 0.5 * (aN[k][21] + aN[k][23] - 2 * aN[k][26]);
        aP[k][8] = 0.5 * (aN[k][20] + aN[k][22] - 2 * aN[k][26]);
        aP[k][9] = 0.5 * (aN[k][24] + aN[k][25] - 2 * aN[k][26]);
        aP[k][10] = 0.125 * (- aN[k][0] + aN[k][1] - aN[k][2] + aN[k][3] + aN[k][4] - aN[k][5] + aN[k][6] - aN[k][7]);
        aP[k][11] = 0.25 * (- aN[k][16] - aN[k][17] + aN[k][18] + aN[k][19] + 2 * aN[k][20] - 2 * aN[k][22]);
        aP[k][12] = 0.25 * (- aN[k][9] - aN[k][11] + aN[k][13] + aN[k][15] + 2 * aN[k][24] - 2 * aN[k][25]);
        aP[k][13] = 0.25 * (- aN[k][10] + aN[k][12] + aN[k][14] + 2 * aN[k][24] - 2 * aN[k][25] - aN[k][8]);
        aP[k][14] = 0.25 * (- aN[k][16] + aN[k][17] + aN[k][18] - aN[k][19] - 2 * aN[k][21] + 2 * aN[k][23]);
        aP[k][15] = 0.25 * (aN[k][9] - aN[k][11] + aN[k][13] - aN[k][15] - 2 * aN[k][21] + 2 * aN[k][23]);
        aP[k][16] = 0.25 * (aN[k][10] - aN[k][12] + aN[k][14] + 2 * aN[k][20] - 2 * aN[k][22] - aN[k][8]);
        aP[k][17] = 0.125 * (aN[k][0] + 2 * aN[k][10] + 2 * aN[k][12] - 2 * aN[k][14] + aN[k][1] - aN[k][2] - aN[k][3] - aN[k][4] - aN[k][5] + aN[k][6] + aN[k][7] - 2 * aN[k][8]);
        aP[k][18] = 0.125 * (aN[k][0] + 2 * aN[k][9] - 2 * aN[k][11] - 2 * aN[k][13] + 2 * aN[k][15] - aN[k][1] - aN[k][2] + aN[k][3] - aN[k][4] + aN[k][5] + aN[k][6] - aN[k][7]);
        aP[k][19] = 0.125 * (aN[k][0] - 2 * aN[k][16] + 2 * aN[k][17] - 2 * aN[k][18] - aN[k][1] + 2 * aN[k][19] + aN[k][2] - aN[k][3] + aN[k][4] - aN[k][5] + aN[k][6] - aN[k][7]);
        aP[k][20] = 0.25 * (aN[k][16] + aN[k][17] + aN[k][18] + aN[k][19] - 2 * (aN[k][20] + aN[k][21] + aN[k][22] + aN[k][23])) + aN[k][26];
        aP[k][21] = 0.25 * (aN[k][9] + aN[k][11] + aN[k][13] + aN[k][15] - 2 * (aN[k][21] + aN[k][23] + aN[k][24] + aN[k][25])) + aN[k][26];
        aP[k][22] = 0.25 * (aN[k][10] + aN[k][12] + aN[k][14] - 2 * (aN[k][20] + aN[k][22] + aN[k][24] + aN[k][25] - 2 * aN[k][26]) + aN[k][8]);
        aP[k][23] = 0.125 * (- aN[k][0] - 2 * aN[k][10] + 2 * aN[k][12] - 2 * aN[k][14] + 2 * aN[k][16] + 2 * aN[k][17] - 2 * aN[k][18] - aN[k][1] -
                             2 * aN[k][19] - 4 * aN[k][20] + 4 * aN[k][22] + aN[k][2] + aN[k][3] - aN[k][4] - aN[k][5] + aN[k][6] + aN[k][7] + 2 * aN[k][8]);
        aP[k][24] = 0.125 * (- aN[k][0] + 2 * aN[k][9] + 2 * aN[k][10] + 2 * aN[k][11] - 2 * aN[k][12] - 2 * aN[k][13] - 2 * aN[k][14] -
                             2 * aN[k][15] - aN[k][1] - 4 * aN[k][24] + 4 * aN[k][25] - aN[k][2] - aN[k][3] + aN[k][4] + aN[k][5] + aN[k][6] + aN[k][7] + 2 * aN[k][8]);
        aP[k][25] = 0.125 * (- aN[k][0] - 2 * aN[k][9] + 2 * aN[k][11] - 2 * aN[k][13] + 2 * aN[k][15] + 2 * aN[k][16] - 2 * aN[k][17] -
                             2 * aN[k][18] + aN[k][1] + 2 * aN[k][19] + 4 * aN[k][21] - 4 * aN[k][23] + aN[k][2] - aN[k][3] - aN[k][4] + aN[k][5] + aN[k][6] - aN[k][7]);
        aP[k][26] = 0.125 * (aN[k][0] - 2 * aN[k][9] - 2 * aN[k][10] - 2 * aN[k][11] - 2 * aN[k][12] - 2 * aN[k][13] - 2 * aN[k][14] - 2 * aN[k][15] -
                             2 * aN[k][16] - 2 * aN[k][17] - 2 * aN[k][18] + aN[k][1] - 2 * aN[k][19] + 4 * aN[k][20] + 4 * aN[k][21] + 4 * aN[k][22] + 4 * aN[k][23] +
                             4 * aN[k][24] + 4 * aN[k][25] - 8 * aN[k][26] + aN[k][2] + aN[k][3] + aN[k][4] + aN[k][5] + aN[k][6] + aN[k][7] - 2 * aN[k][8]);
      }
    }
  }
  

  void GetHexPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) {

    const unsigned nDofs = hexNumberOfDofs[solType];

    phi.resize(nDofs);


    //common for linear, quadratic and biquadratic
    phi[0] = 1.;
    phi[1] = xi[0]; // x
    phi[2] = xi[1]; // y
    phi[3] = xi[2]; // z
    phi[4] = phi[1] * phi[2]; // x y
    phi[5] = phi[1] * phi[3]; // x z
    phi[6] = phi[2] * phi[3]; // y z
    if(solType < 1) { //only linear
      phi[7] = phi[4] * phi[3];  // x y z
    }
    else { //only quadratic and biquadratic
      phi[7] = phi[1] * phi[1];   // x x
      phi[8] = phi[2] * phi[2];   // y y
      phi[9] = phi[3] * phi[3];   // z z
      phi[10] = phi[4] * phi[3];  // x y z
      phi[11] = phi[7] * phi[2];  // xx y
      phi[12] = phi[7] * phi[3];  // xx z
      phi[13] = phi[8] * phi[3];  // yy z
      phi[14] = phi[1] * phi[8];  // x yy
      phi[15] = phi[1] * phi[9];  // x zz
      phi[16] = phi[2] * phi[9];  // y zz
      phi[17] = phi[11] * phi[3]; // xx y z
      phi[18] = phi[1] * phi[13]; // x yy z
      phi[19] = phi[1] * phi[16]; // x y zz

      if(solType > 1) { //only biquadratic
        phi[20] = phi[7] * phi[8];  // xx yy
        phi[21] = phi[7] * phi[9];  // xx zz
        phi[22] = phi[8] * phi[9];  // yy zz
        phi[23] = phi[11] * phi[9]; // xx y zz
        phi[24] = phi[7] * phi[13]; // xx yy z
        phi[25] = phi[14] * phi[9]; // x yy zz
        phi[26] = phi[20] * phi[9]; // xx yy zz
      }
    }

  }

  void GetHexPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      const std::vector < double >& xi, const unsigned & solType) {

    GetHexPolynomialShapeFunction(phi,  xi, solType);

    const unsigned dim = 3;

    const unsigned nDofs = hexNumberOfDofs[solType];

    gradPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      gradPhi[i].assign(dim, 0.);
    }

    //common for linear, quadratic and biquadratic
    //phi_x
    gradPhi[1][0] = 1.; // 1
    gradPhi[4][0] = xi[1] ; // y
    gradPhi[5][0] = xi[2] ; // z
    //phi_y
    gradPhi[2][1] = 1.; // 1
    gradPhi[4][1] = xi[0]; // x
    gradPhi[6][1] = xi[2]; // z

    //phi_z
    gradPhi[3][2] = 1.; // 1
    gradPhi[5][2] = xi[0]; // x
    gradPhi[6][2] = xi[1]; // y

    if(solType < 1) {  //only linear
      gradPhi[7][0] = phi[6];  // y z
      gradPhi[7][1] = phi[5];  // x z
      gradPhi[7][2] = phi[4];  // x y
    }
    else { //only quadratic and biquadratic
      //phi_x
      gradPhi[7][0] = 2 * xi[0];   // 2 x
      gradPhi[10][0] = phi[6];  // y z
      gradPhi[11][0] = 2 * phi[4];  // 2 x y
      gradPhi[12][0] = 2 * phi[5];  // 2 x z
      gradPhi[14][0] = phi[8];  // yy
      gradPhi[15][0] = phi[9];  // zz
      gradPhi[17][0] = 2 * phi[10]; // 2 x y z
      gradPhi[18][0] = phi[13]; //  yy z
      gradPhi[19][0] = phi[16]; // y zz
      //phi_y
      gradPhi[8][1] = 2 * xi[1];   // 2 y
      gradPhi[10][1] = phi[5];  // x z
      gradPhi[11][1] = phi[7];  // xx
      gradPhi[13][1] = 2 * phi[6];  // 2 y z
      gradPhi[14][1] = 2 * phi[4];  // 2 x y
      gradPhi[16][1] = phi[9];  // zz
      gradPhi[17][1] = phi[12]; // xx  z
      gradPhi[18][1] = 2 * phi[10]; // 2 x y z
      gradPhi[19][1] = phi[15]; // x zz
      //phi_z
      gradPhi[9][2] = 2 * xi[2];   // 2 z
      gradPhi[10][2] = phi[4];  // x y
      gradPhi[12][2] = phi[7];  // xx
      gradPhi[13][2] = phi[8];  // yy
      gradPhi[15][2] = 2 * phi[5];  // 2 x z
      gradPhi[16][2] = 2 * phi[6];  // 2 y z
      gradPhi[17][2] = phi[11]; // xx y
      gradPhi[18][2] = phi[14]; // x yy
      gradPhi[19][2] = 2 * phi[10]; // 2 x y z

      if(solType > 1) { //only biquadratic
        //phi_x
        gradPhi[20][0] = 2 * phi[14];  // 2 x yy
        gradPhi[21][0] = 2 * phi[15];  // 2 x zz
        gradPhi[23][0] = 2 * phi[19]; // 2 x y zz
        gradPhi[24][0] = 2 * phi[18]; // 2 x yy z
        gradPhi[25][0] = phi[22]; //  yy zz
        gradPhi[26][0] = 2 * phi[25]; // 2 x yy zz
        //phi_y
        gradPhi[20][1] = 2 * phi[11];  // 2 xx y
        gradPhi[22][1] = 2 * phi[16];  // 2 y zz
        gradPhi[23][1] = phi[21]; // xx zz
        gradPhi[24][1] = 2 * phi[17]; // 2 xx y z
        gradPhi[25][1] = 2 * phi[19]; // 2 x y zz
        gradPhi[26][1] = 2 * phi[23]; // 2 xx y zz
        //phi_z
        gradPhi[21][2] = 2 * phi[12];  // 2 xx z
        gradPhi[22][2] = 2 * phi[13];  // 2 yy z
        gradPhi[23][2] = 2 * phi[17]; // 2 xx y z
        gradPhi[24][2] = phi[20]; // xx yy
        gradPhi[25][2] = 2 * phi[18]; // 2 x yy z
        gradPhi[26][2] = 2 * phi[24]; // 2 xx yy z
      }
    }



  }

  void GetHexPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) {

    GetHexPolynomialShapeFunctionGradient(phi,   gradPhi,  xi, solType);

    const unsigned dim = 3;
    const unsigned nDofs = hexNumberOfDofs[solType];

    hessPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      hessPhi[i].resize(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessPhi[i][i1].assign(dim, 0.);
      }
    }

    // common for linear, quadratic and biquadratic
    //phi_xy
    hessPhi[4][1][0] = hessPhi[4][0][1] = 1.; // 1
    //phi_xz
    hessPhi[5][2][0] = hessPhi[5][0][2] = 1.; // 1
    //phi_yz
    hessPhi[6][2][1] = hessPhi[6][1][2] = 1.; // 1
    if(solType < 1) { //only linear
      //phi_xy
      hessPhi[7][1][0] = hessPhi[7][0][1] = xi[2];  //  z
      //phi_xz
      hessPhi[7][2][0] = hessPhi[7][0][2] = xi[1];  //  y
      //phi_yz
      hessPhi[7][2][1] = hessPhi[7][1][2] = xi[0] ; //  x
    }
    else { //only quadraic and biquadratic
      //phi_xx
      hessPhi[7][0][0] = 2.;   // 2
      hessPhi[11][0][0] = 2 * xi[1];  // 2 y
      hessPhi[12][0][0] = 2 * xi[2];  // 2 z
      hessPhi[17][0][0] = 2 * phi[6]; // 2 y z
      //phi_xy
      hessPhi[10][1][0] = hessPhi[10][0][1] = xi[2];  // z
      hessPhi[11][1][0] = hessPhi[11][0][1] = 2 * xi[0];  // 2 x
      hessPhi[14][1][0] = hessPhi[14][0][1] = 2 * xi[1];  // 2 y
      hessPhi[17][1][0] = hessPhi[17][0][1] = 2 * phi[5]; // 2 x z
      hessPhi[18][1][0] = hessPhi[18][0][1] = 2 * phi[6]; //  2 y z
      hessPhi[19][1][0] = hessPhi[19][0][1] = phi[9]; // zz
      //phi_xz
      hessPhi[10][2][0] = hessPhi[10][0][2] = xi[1];  // y
      hessPhi[12][2][0] = hessPhi[12][0][2] = 2 * xi[0];  // 2 x
      hessPhi[15][2][0] = hessPhi[15][0][2] = 2 * xi[2];  // 2 z
      hessPhi[17][2][0] = hessPhi[17][0][2] = 2 * phi[4]; // 2 x y
      hessPhi[18][2][0] = hessPhi[18][0][2] = phi[8]; //  yy
      hessPhi[19][2][0] = hessPhi[19][0][2] = 2 * phi[6]; // 2 y z
      //phi_yy
      hessPhi[8][1][1] = 2.;   // 2
      hessPhi[13][1][1] = 2. * xi[2];  // 2 z
      hessPhi[14][1][1] = 2. * xi[0];  // 2 x
      hessPhi[18][1][1] = 2. * phi[5]; // 2 x z
      //phi_yz
      hessPhi[10][2][1] = hessPhi[10][1][2] = xi[0];  // x
      hessPhi[13][2][1] = hessPhi[13][1][2] = 2 * xi[1];  // 2 y
      hessPhi[16][2][1] = hessPhi[16][1][2] = 2 * xi[2];  // 2 z
      hessPhi[17][2][1] = hessPhi[17][1][2] = phi[7]; // xx
      hessPhi[18][2][1] = hessPhi[18][1][2] = 2 * phi[4]; // 2 x y
      hessPhi[19][2][1] = hessPhi[19][1][2] = 2 * phi[5]; // 2 x z
      //phi_zz
      hessPhi[9][2][2] = 2.;   // 2
      hessPhi[15][2][2] = 2 * xi[0];  // 2 x
      hessPhi[16][2][2] = 2 * xi[1];  // 2 y
      hessPhi[19][2][2] = 2 * phi[4]; // 2 x y

      if(solType > 1) { // only biquadratic
        //phi_xx
        hessPhi[20][0][0] = 2 * phi[8];  // 2 yy
        hessPhi[21][0][0] = 2 * phi[9];  // 2 zz
        hessPhi[23][0][0] = 2 * phi[16]; // 2  y zz
        hessPhi[24][0][0] = 2 * phi[13]; // 2  yy z
        hessPhi[26][0][0] = 2 * phi[22]; // 2  yy zz
        //phi_xy
        hessPhi[20][1][0] = hessPhi[20][0][1] = 4 * phi[4];  // 4 x y
        hessPhi[23][1][0] = hessPhi[23][0][1] = 2 * phi[15]; // 2 x zz
        hessPhi[24][1][0] = hessPhi[24][0][1] = 4 * phi[10]; // 4 x y z
        hessPhi[25][1][0] = hessPhi[25][0][1] = 2 * phi[16]; //  2 y zz
        hessPhi[26][1][0] = hessPhi[26][0][1] = 4 * phi[19]; // 4 x y zz
        //phi_xz
        hessPhi[21][2][0] = hessPhi[21][0][2] = 4 * phi[5];  // 4 x z
        hessPhi[23][2][0] = hessPhi[23][0][2] = 4 * phi[10]; // 4 x y z
        hessPhi[24][2][0] = hessPhi[24][0][2] = 2 * phi[14]; // 2 x yy
        hessPhi[25][2][0] = hessPhi[25][0][2] = 2 * phi[13]; //  2 yy z
        hessPhi[26][2][0] = hessPhi[26][0][2] = 4 * phi[18]; // 4 x yy z
        //phi_yy
        hessPhi[20][1][1] = 2 * phi[7];  // 2 xx
        hessPhi[22][1][1] = 2 * phi[9];  // 2 zz
        hessPhi[24][1][1] = 2 * phi[12]; // 2 xx z
        hessPhi[25][1][1] = 2 * phi[15]; // 2 x zz
        hessPhi[26][1][1] = 2 * phi[21]; // 2 xx zz
        //phi_yz
        hessPhi[22][2][1] = hessPhi[22][1][2] = 4 * phi[6];  // 4 y z
        hessPhi[23][2][1] = hessPhi[23][1][2] = 2 * phi[12]; // 2 xx z
        hessPhi[24][2][1] = hessPhi[24][1][2] = 2 * phi[11]; // 2 xx y
        hessPhi[25][2][1] = hessPhi[25][1][2] = 4 * phi[10]; // 4 x y z
        hessPhi[26][2][1] = hessPhi[26][1][2] = 4 * phi[17]; // 4 xx y z
        //phi_zz
        hessPhi[21][2][2] = 2 * phi[7];  // 2 xx
        hessPhi[22][2][2] = 2 * phi[8];  // 2 yy
        hessPhi[23][2][2] = 2 * phi[11]; // 2 xx y
        hessPhi[25][2][2] = 2 * phi[14]; // 2 x yy
        hessPhi[26][2][2] = 2 * phi[20]; // 2 xx yy
      }
    }

  }
//END HEX

//BEGIN TET specialized functions
  const unsigned tetNumberOfDofs[3] = {4, 10, 15};

  void ProjectTetNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solType) {

    unsigned dim =  aN.size();
    aP.resize(dim);

    unsigned nDofs = tetNumberOfDofs[solType];
    if(nDofs > aN[0].size()) {
      std::cout << "Error in ProjectTetNodalToPolynomialCoefficients(...) the number of Dofs is inconsistent" << std::endl;
      abort();
    }

    for(unsigned k = 0; k < dim; k++) {
      aP[k].resize(nDofs);
    }

    if(solType == 0) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][0] ;
        aP[k][1] = - aN[k][0] + aN[k][1] ;
        aP[k][2] = - aN[k][0] + aN[k][2] ;
        aP[k][3] = - aN[k][0] + aN[k][3] ;
      }
    }
    else if(solType == 1) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][0];
        aP[k][1] =  - 3 * aN[k][0] - aN[k][1] + 4 * aN[k][4];
        aP[k][2] = - 3 * aN[k][0] - aN[k][2] + 4 * aN[k][6];
        aP[k][3] =  - 3 * aN[k][0] - aN[k][3] + 4 * aN[k][7];
        aP[k][4] = 4 * aN[k][0] - 4 * aN[k][4] + 4 * aN[k][5] - 4 * aN[k][6];
        aP[k][5] = 4 * aN[k][0] - 4 * aN[k][4] - 4 * aN[k][7] + 4 * aN[k][8];
        aP[k][6] = 4 * aN[k][0] + 4 * aN[k][9] - 4 * aN[k][6] - 4 * aN[k][7];
        aP[k][7] =  2 * aN[k][0] + 2 * aN[k][1] - 4 * aN[k][4];
        aP[k][8] = 2 * aN[k][0] + 2 * aN[k][2] - 4 * aN[k][6];
        aP[k][9] = 2 * aN[k][0] + 2 * aN[k][3] - 4 * aN[k][7];
      }
    }
    else if(solType == 2) {
      for(int k = 0; k < dim; k++) {
        aP[k][14] = 4 * (aN[k][0] - 8 * aN[k][9] + 27 * aN[k][10] + 27 * aN[k][11] + 27 * aN[k][12] + 27 * aN[k][13] - 64 * aN[k][14] +
                         aN[k][1] + aN[k][2] + aN[k][3] - 8 * (aN[k][4] + aN[k][5] + aN[k][6] + aN[k][7] + aN[k][8]));
        aP[k][13] = -3 * (aN[k][0] - 4 * aN[k][9] + 9 * aN[k][13] + aN[k][2] + aN[k][3] - 4 * (aN[k][6] + aN[k][7]));
        aP[k][12] = -3  * (aN[k][0] + 9 * aN[k][11] + aN[k][1] + aN[k][3] - 4 * (aN[k][4] + aN[k][7] + aN[k][8]));
        aP[k][11] = -3 * (aN[k][0] + 9  * aN[k][10] + aN[k][1] + aN[k][2] - 4 * (aN[k][4] + aN[k][5] + aN[k][6]));
        aP[k][10] = -13 * aN[k][0] + 32  * aN[k][9] - 135 * aN[k][10] - 135 * aN[k][11] - 81  * aN[k][12] - 135 * aN[k][13] +
                    256 * aN[k][14] - 7 * aN[k][1] - 7 * aN[k][2] - 7 * aN[k][3] + 8 * (7 * aN[k][4] + 4 * aN[k][5] + 7 * (aN[k][6] + aN[k][7]) + 4 * aN[k][8]);
        aP[k][9] = 2 * (aN[k][0] + aN[k][3] - 2 * aN[k][7]);
        aP[k][8] = 2 * (aN[k][0] + aN[k][2] - 2 * aN[k][6]);
        aP[k][7] = 2 * (aN[k][0] + aN[k][1] - 2 * aN[k][4]);
        aP[k][6] = 7 * aN[k][0] - 8 * aN[k][9] + 27 * aN[k][13] + 3 * aN[k][2] + 3 * aN[k][3] - 16 * (aN[k][6] + aN[k][7]);
        aP[k][5] = 7 * aN[k][0] + 27 * aN[k][11] + 3 * aN[k][1] + 3 * aN[k][3] - 8 * (2 * (aN[k][4] + aN[k][7]) + aN[k][8]);
        aP[k][4] = 7 * aN[k][0] + 27 * aN[k][10] + 3 * aN[k][1] + 3 * aN[k][2] - 8 * (2 * aN[k][4] + aN[k][5] + 2 * aN[k][6]);
        aP[k][3] = -3 * aN[k][0] - aN[k][3] + 4 * aN[k][7];
        aP[k][2] = -3 * aN[k][0] - aN[k][2] + 4 * aN[k][6];
        aP[k][1] = -3 * aN[k][0] - aN[k][1] + 4 * aN[k][4];
        aP[k][0] = aN[k][0];
      }
    }
  }
  

  void GetTetPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) {

    const unsigned nDofs = tetNumberOfDofs[solType];

    phi.resize(nDofs);

    phi[0] = 1.;
    phi[1] = xi[0]; // x
    phi[2] = xi[1];  // y
    phi[3] = xi[2]; // z

    if(solType > 0) {
      phi[4] = phi[1] * phi[2];  // x y
      phi[5] = phi[1] * phi[3];  // x z
      phi[6] = phi[2] * phi[3];  // y z
      phi[7] = phi[1] * phi[1]; // x x
      phi[8] = phi[2] * phi[2]; // y y
      phi[9] = phi[3] * phi[3]; // z z


      if(solType > 1) {
        phi[10] = phi[4] * phi[3]; // x y z
        phi[11] = phi[1] * phi[8] + phi[7] * phi[2] ; // x y y + x x y
        phi[12] = phi[1] * phi[9] + phi[7] * phi[3] ; // x z z + x x z
        phi[13] = phi[2] * phi[9] + phi[8] * phi[3] ; // y z z + y y z
        phi[14] = phi[4] * phi[9] + phi[4] * phi[6] + phi[7] * phi[6]; // x y z z + x y y z + x x y z
      }
    }

  }


  void GetTetPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      const std::vector < double >& xi, const unsigned & solType) {

    GetTetPolynomialShapeFunction(phi,  xi, solType);

    const unsigned dim = 3;

    const unsigned nDofs = tetNumberOfDofs[solType];

    gradPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      gradPhi[i].assign(dim, 0.);
    }

    //phi_x
    gradPhi[1][0] = 1.; // 1
    //phi_y
    gradPhi[2][1] = 1.;  // 1
    //phi_z
    gradPhi[3][2] = 1.; // 1

    if(solType > 0) {
      //phi_x
      gradPhi[4][0] = xi[1] ;  //  y
      gradPhi[5][0] = xi[2] ;  //  z
      gradPhi[7][0] = 2 * xi[0] ; // 2 x
      //phi_y
      gradPhi[4][1] = xi[0] ;  // x
      gradPhi[6][1] = xi[2] ;  //  z
      gradPhi[8][1] = 2 * xi[1]; // 2 y
      //phi_z
      gradPhi[5][2] = xi[0];  // x
      gradPhi[6][2] = xi[1];  // y
      gradPhi[9][2] = 2 * xi[2]; // 2 z

      if(solType > 1) {
        //phi_x
        gradPhi[10][0] = phi[6]; //  y z
        gradPhi[11][0] = phi[8] + 2 * phi[4] ; //  y y + 2 x y
        gradPhi[12][0] = phi[9] + 2 * phi[5] ; // z z + 2 x z
        gradPhi[14][0] = xi[1] * (gradPhi[12][0] + phi[6]) ; //  y z z +  y y z + 2 x y z
        //phi_y
        gradPhi[10][1] = phi[5]; // x z
        gradPhi[11][1] = 2 * phi[4] + phi[7] ; //  2 x y + x x
        gradPhi[13][1] = phi[9] + 2 * phi[6]; //  z z + 2 y z
        gradPhi[14][1] = xi[0] * (gradPhi[13][1] + phi[5]); // x z z + 2 x y z + x x z
        //phi_z
        gradPhi[10][2] = phi[4]; // x y
        gradPhi[12][2] = 2 * phi[5] + phi[7]; // 2 x z  + x x
        gradPhi[13][2] = 2 * phi[6] + phi[8]; // 2 y z  + y y
        gradPhi[14][2] = xi[0] * (gradPhi[13][2] + phi[4]); //  2 x y z  + x y y  + x x y
      }
    }
  }

  void GetTetPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) {

    GetTetPolynomialShapeFunctionGradient(phi,   gradPhi,  xi, solType);

    const unsigned dim = 3;
    const unsigned nDofs = tetNumberOfDofs[solType];

    hessPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      hessPhi[i].resize(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessPhi[i][i1].assign(dim, 0.);
      }
    }

    if(solType > 0) {
      //phi_xx
      hessPhi[7][0][0] = 2.; // 2
      //phi_xy
      hessPhi[4][1][0] = hessPhi[4][0][1] = 1.;  //  1
      //phi_xz
      hessPhi[5][2][0] = hessPhi[5][0][2] = 1.;  //  1
      //phi_yy
      hessPhi[8][1][1] = 2.; // 2
      //phi_yz
      hessPhi[6][2][1] = hessPhi[6][1][2] = 1.;  //  1
      //phi_zz
      hessPhi[9][2][2] = 2.; // 2

      if(solType > 1) {
        //phi_xx
        hessPhi[11][0][0] = 2 * xi[1]; //  2 y
        hessPhi[12][0][0] = 2 * xi[2]; // 2 z
        hessPhi[14][0][0] = 2 * phi[6]; //  2 y z
        //phi_xy
        hessPhi[10][1][0] = hessPhi[10][0][1] = xi[2] ; //  z
        hessPhi[11][1][0] = hessPhi[11][0][1] = 2 * xi[1] + 2 * xi[0] ; //  2 y + 2 x
        hessPhi[14][1][0] = hessPhi[14][0][1] = gradPhi[13][1] + 2 * phi[5] ; //  z z +  2 y z + 2 x z
        //phi_xz
        hessPhi[10][2][0] = hessPhi[10][0][2] = xi[1] ; //  y
        hessPhi[12][2][0] = hessPhi[12][0][2] = 2 * (xi[0] + xi[2]) ; // 2 z + 2 x
        hessPhi[14][2][0] = hessPhi[14][0][2] = gradPhi[13][2] +  2 * phi[4] ; //  2 y z +  y y  + 2 x y
        //phi_yy
        hessPhi[11][1][1] = 2 * xi[0] ; //  2 x
        hessPhi[13][1][1] = 2 * xi[2]; //  2 z
        hessPhi[14][1][1] = 2 * phi[5]; // 2 x z
        //phi_yz
        hessPhi[10][2][1] = hessPhi[10][1][2] = xi[0] ; // x
        hessPhi[13][2][1] = hessPhi[13][1][2] = 2 * xi[2] + 2 * xi[1]; //  2 z + 2 y
        hessPhi[14][2][1] = hessPhi[14][1][2] = 2 * phi[5] + gradPhi[11][1] ; // 2 x z  + 2 x y  + x x
        //phi_zz
        hessPhi[12][2][2] = 2 * xi[0] ; // 2 x
        hessPhi[13][2][2] = 2 * xi[1] ; // 2 y
        hessPhi[14][2][2] = 2 * phi[4]; //  2 x y
      }
    }
  }
//END TET

//BEGIN WEDGE specialized functions
  const unsigned wedgeNumberOfDofs[3] = {6, 15, 21};

  void ProjectWedgeNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solType) {

    unsigned dim =  aN.size();
    aP.resize(dim);

    unsigned nDofs = wedgeNumberOfDofs[solType];
    if(nDofs > aN[0].size()) {
      std::cout << "Error in ProjectWedgeNodalToPolynomialCoefficients(...) the number of Dofs is inconsistent" << std::endl;
      abort();
    }

    for(unsigned k = 0; k < dim; k++) {
      aP[k].resize(nDofs);
    }

    if(solType == 0) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = 0.5 * (aN[k][0] + aN[k][3]);
        aP[k][1] = 0.5 * (- aN[k][0] + aN[k][1] - aN[k][3] + aN[k][4]);
        aP[k][2] = 0.5 * (- aN[k][0] + aN[k][2] - aN[k][3] + aN[k][5]);
        aP[k][3] = 0.5 * (- aN[k][0] + aN[k][3]);
        aP[k][4] = 0.5 * (aN[k][0] - aN[k][1] - aN[k][3] + aN[k][4]);
        aP[k][5] = 0.5 * (aN[k][0] - aN[k][2] - aN[k][3] + aN[k][5]);
      }
    }
    else if(solType == 1) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][12];
        aP[k][1] = - aN[k][0] + 2 * aN[k][9] - aN[k][12] + aN[k][13] - aN[k][1] - aN[k][3] - aN[k][4] + 2 * aN[k][6];
        aP[k][2] = - aN[k][0] + 2 * aN[k][11] - aN[k][12] + aN[k][14] - aN[k][2] - aN[k][3] - aN[k][5] + 2 * aN[k][8];
        aP[k][3] = 0.5 * (- aN[k][0] + aN[k][3]);
        aP[k][4] = 2 * (aN[k][0] - aN[k][9] + aN[k][10] - aN[k][11] + aN[k][3] - aN[k][6] + aN[k][7] - aN[k][8]);
        aP[k][5] = 0.5 * (3 * aN[k][0] + 4 * aN[k][9] + aN[k][1] - 3 * aN[k][3] - aN[k][4] - 4 * aN[k][6]);
        aP[k][6] = 0.5 * (3 * aN[k][0] + 4 * aN[k][11] + aN[k][2] - 3 * aN[k][3] - aN[k][5] - 4 * aN[k][8]);
        aP[k][7] = aN[k][0] - 2 * aN[k][9] + aN[k][1] + aN[k][3] + aN[k][4] - 2 * aN[k][6];
        aP[k][8] = aN[k][0] - 2 * aN[k][11] + aN[k][2] + aN[k][3] + aN[k][5] - 2 * aN[k][8];
        aP[k][9] = 0.5 * (aN[k][0] - 2 * aN[k][12] + aN[k][3]);
        aP[k][10] = -2 * (aN[k][0] + aN[k][9] - aN[k][10] + aN[k][11] - aN[k][3] - aN[k][6] + aN[k][7] - aN[k][8]);
        aP[k][11] = - aN[k][0] - 2 * aN[k][9] - aN[k][1] + aN[k][3] + aN[k][4] + 2 * aN[k][6];
        aP[k][12] = - aN[k][0] - 2 * aN[k][11] - aN[k][2] + aN[k][3] + aN[k][5] + 2 * aN[k][8];
        aP[k][13] = 0.5 * (- aN[k][0] + 2 * aN[k][12] - 2 * aN[k][13] + aN[k][1] - aN[k][3] + aN[k][4]);
        aP[k][14] = 0.5 * (- aN[k][0] + 2 * aN[k][12] - 2 * aN[k][14] + aN[k][2] - aN[k][3] + aN[k][5]);
      }
    }
    else if(solType == 2) {
      for(int k = 0; k < dim; k++) {
        aP[k][0] = aN[k][12];
        aP[k][1] = - 3 * aN[k][12] - aN[k][13] + 4 * aN[k][15];
        aP[k][2] = - 3 * aN[k][12] - aN[k][14] + 4 * aN[k][17];
        aP[k][3] = 0.5 * (- aN[k][0] + aN[k][3]);
        aP[k][4] = 7 * aN[k][12] + 3 * aN[k][13] + 3 * aN[k][14] - 8 * (2 * aN[k][15] + aN[k][16] + 2 * aN[k][17]) + 27 * aN[k][20];
        aP[k][5] = 0.5 * (3 * aN[k][0] + 4 * aN[k][9] + aN[k][1] - 3 * aN[k][3] - aN[k][4] - 4 * aN[k][6]);
        aP[k][6] = 0.5 * (3 * aN[k][0] + 4 * aN[k][11] + aN[k][2] - 3 * aN[k][3] - aN[k][5] - 4 * aN[k][8]);
        aP[k][7] = 2 * (aN[k][12] + aN[k][13] - 2 * aN[k][15]);
        aP[k][8] = 2 * (aN[k][12] + aN[k][14] - 2 * aN[k][17]);
        aP[k][9] = 0.5 * (aN[k][0] - 2 * aN[k][12] + aN[k][3]);
        aP[k][10] = 0.5 * (- 7 * aN[k][0] - 16 * aN[k][9] - 8 * aN[k][10] - 16 * aN[k][11] - 27 * aN[k][18] -
                           3 * aN[k][1] + 27 * aN[k][19] - 3 * aN[k][2] + 7 * aN[k][3] + 3 * aN[k][4] + 3 * aN[k][5] +
                           8 * (2 * aN[k][6] + aN[k][7] + 2 * aN[k][8])) ;
        aP[k][11] = - 3 * (aN[k][12] + aN[k][13] + aN[k][14] - 4 * (aN[k][15] + aN[k][16] + aN[k][17]) + 9 * aN[k][20]);
        aP[k][12] = - aN[k][0] - 2 * aN[k][9] - aN[k][1] + aN[k][3] + aN[k][4] + 2 * aN[k][6];
        aP[k][13] = - aN[k][0] - 2 * aN[k][11] - aN[k][2] + aN[k][3] + aN[k][5] + 2 * aN[k][8];
        aP[k][14] = 0.5 * (-3 * aN[k][0] + 4 * aN[k][9] + 6 * aN[k][12] + 2 * aN[k][13] - 8 * aN[k][15] - aN[k][1] - 3 * aN[k][3] - aN[k][4] + 4 * aN[k][6]);
        aP[k][15] = 0.5 * (-3 * aN[k][0] + 4 * aN[k][11] + 6 * aN[k][12] + 2 * aN[k][14] - 8 * aN[k][17] - aN[k][2] - 3 * aN[k][3] - aN[k][5] + 4 * aN[k][8]);
        aP[k][16] = aN[k][0] - 2 * aN[k][9] - 2 * aN[k][12] - 2 * aN[k][13] + 4 * aN[k][15] + aN[k][1] + aN[k][3] + aN[k][4] - 2 * aN[k][6];
        aP[k][17] = aN[k][0] - 2 * aN[k][11] - 2 * aN[k][12] - 2 * aN[k][14] + 4 * aN[k][17] + aN[k][2] + aN[k][3] + aN[k][5] - 2 * aN[k][8];
        aP[k][18] = 1.5 * (aN[k][0] + 4 * aN[k][9] + 4 * aN[k][10] + 4 * aN[k][11] + 9 * aN[k][18] + aN[k][1] - 9 * aN[k][19] + aN[k][2] - aN[k][3] - aN[k][4] - aN[k][5] - 4 * (aN[k][6] + aN[k][7] + aN[k][8]));
        aP[k][19] = 0.5 * (7 * aN[k][0] - 16 * aN[k][9] - 8 * aN[k][10] - 16 * aN[k][11] - 14 * aN[k][12] - 6 * aN[k][13] - 6 * aN[k][14] + 32 * aN[k][15] +
                           16 * aN[k][16] + 32 * aN[k][17] + 27 * aN[k][18] + 3 * aN[k][1] + 27 * aN[k][19] - 54 * aN[k][20] + 3 * aN[k][2] + 7 * aN[k][3] +
                           3 * aN[k][4] + 3 * aN[k][5] - 8 * (2 * aN[k][6] + aN[k][7] + 2 * aN[k][8]));
        aP[k][20] = -1.5 * (aN[k][0] - 4 * aN[k][9] - 4 * aN[k][10] - 4 * aN[k][11] - 2 * aN[k][12] - 2 * aN[k][13] - 2 * aN[k][14] + 8 * aN[k][15] + 8 * aN[k][16] + 8 * aN[k][17] + 9 * aN[k][18] + aN[k][1] + 9 * aN[k][19] - 18 * aN[k][20] +
                            aN[k][2] + aN[k][3] + aN[k][4] + aN[k][5] - 4 * (aN[k][6] + aN[k][7] + aN[k][8]));
      }
    }
  }
  

  void GetWedgePolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) {

    const unsigned nDofs = wedgeNumberOfDofs[solType];

    phi.resize(nDofs);

    //common for linear, quadratic and biquadratic
    phi[0] = 1.;
    phi[1] = xi[0]; // x
    phi[2] = xi[1];  // y
    phi[3] = xi[2]; // z
    if(solType < 1) {  //only linear
      phi[4] = phi[1] * phi[3];  // x z
      phi[5] = phi[2] * phi[3];  // y z
    }
    else {  //only quadratic and biquadratic
      phi[4] = xi[0] * xi[1];  // x y
      phi[5] = xi[0] * xi[2];  // x z
      phi[6] = xi[1] * xi[2];  // y z
      phi[7] = xi[0] * xi[0]; // x x
      phi[8] = xi[1] * xi[1]; // y y
      phi[9] = xi[2] * xi[2]; // z z
      phi[10] = phi[4] * xi[2]; // x y z
      if(solType < 2) {  //only quadratic
        phi[11] = phi[7] * xi[2]; // x x z
        phi[12] = phi[8] * xi[2]; // y y z
        phi[13] = xi[0] * phi[9]; // x z z
        phi[14] = xi[1] * phi[9]; // y z z;
      }
      else { //only biquadratic
        phi[11] = phi[7] * xi[1] + xi[0] * phi[8]; // xx y + x yy
        phi[12] = phi[7] * xi[2]; // x x z
        phi[13] = phi[8] * xi[2]; // y y z
        phi[14] = xi[0]  * phi[9]; // x z z
        phi[15] = xi[1] * phi[9]; // y z z
        phi[16] = phi[7] * phi[9]; // xx zz
        phi[17] = phi[8] * phi[9]; // yy zz
        phi[18] = phi[11] * xi[2]; // x yy z + xx y z
        phi[19] = phi[10] * xi[2]; // x y zz
        phi[20] = phi[18] * xi[2]; // x yy zz + xx y zz
      }
    }
  }

  void GetWedgePolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      const std::vector < double >& xi, const unsigned & solType) {

    GetWedgePolynomialShapeFunction(phi,  xi, solType);

    const unsigned dim = 3;

    const unsigned nDofs = wedgeNumberOfDofs[solType];

    gradPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      gradPhi[i].assign(dim, 0.);
    }

    //common for linear, quadratic and biquadratic
    //phi_x
    gradPhi[1][0] = 1.; // 1
    //phi_y
    gradPhi[2][1] = 1.;  // 1
    //phi_z
    gradPhi[3][2] = 1.; // 1
    if(solType < 1) { //only linear
      //phi_x
      gradPhi[4][0] = xi[2] ;  //  z
      //phi_y
      gradPhi[5][1] = xi[2] ;  // z
      //phi_z
      gradPhi[4][2] = xi[0] ;  // x
      gradPhi[5][2] = xi[1] ;  // y

    }
    else { //only quadratic and biquadratic
      //phi_x
      gradPhi[4][0] = xi[1] ;  //  y
      gradPhi[5][0] = xi[2] ;  //  z
      gradPhi[7][0] = 2 * xi[0] ; // 2 x
      gradPhi[10][0] = phi[6]; //  y z

      //phi_y
      gradPhi[4][1] = xi[0] ;  //  x
      gradPhi[6][1] = xi[2] ;  // z
      gradPhi[8][1] = 2 * xi[1]; // 2 y
      gradPhi[10][1] = phi[5]; // x z

      //phi_z
      gradPhi[5][2] = xi[0] ;  // x
      gradPhi[6][2] = xi[1] ;  // y
      gradPhi[9][2] = 2 * xi[2]; // 2 z
      gradPhi[10][2] = phi[4]; // x y

      if(solType < 2) {  //only quadratic
        //phi_x
        gradPhi[11][0] = 2 * phi[5]; // 2 x z
        gradPhi[13][0] = phi[9]; //  z z

        //phi_y
        gradPhi[12][1] = 2 * phi[6]; // 2 y z
        gradPhi[14][1] = phi[9]; //  z z;


        //phi_z
        gradPhi[11][2] = phi[7]; // x x
        gradPhi[12][2] = phi[8]; // y y

        gradPhi[13][2] = 2 * phi[5]; // 2 x z
        gradPhi[14][2] = 2 * phi[6]; // 2 y z;
      }
      else { //only biquadratic
        //phi_x
        gradPhi[11][0] = 2 * phi[4] + phi[8]; // 2 x y +  yy
        gradPhi[12][0] = 2 * phi[5]; // 2 x z
        gradPhi[14][0] = phi[9]; //  z z
        gradPhi[16][0] = 2 * phi[14]; // 2 x zz
        gradPhi[18][0] = phi[13] + 2 * phi[10]; //  yy z + 2 x y z
        gradPhi[19][0] = phi[15]; // y zz
        gradPhi[20][0] = gradPhi[18][0] * xi[2]; //  yy zz + 2 x y zz
        //phi_y
        gradPhi[11][1] = phi[7] + 2 * phi[4]; // xx + 2 x y
        gradPhi[13][1] = 2 * phi[6]; // 2 y z
        gradPhi[15][1] = phi[9]; //  z z
        gradPhi[17][1] = 2 * phi[15]; // 2 y zz
        gradPhi[18][1] = 2 * phi[10] + phi[12]; // 2 x y z + xx z
        gradPhi[19][1] = phi[14]; // x zz
        gradPhi[20][1] = gradPhi[18][1] * xi[2]; // 2 x y zz + xx zz
        //phi_z
        gradPhi[12][2] = phi[7]; // x x
        gradPhi[13][2] = phi[8]; // y y
        gradPhi[14][2] = 2  * phi[5]; // 2 x z
        gradPhi[15][2] = 2 * phi[6]; // 2 y z
        gradPhi[16][2] = 2 * phi[12]; // 2 xx z
        gradPhi[17][2] = 2 * phi[13]; // 2 yy z
        gradPhi[18][2] = xi[1] * (phi[4] + phi[7]); // x yy + xx y
        gradPhi[19][2] = 2 * phi[10]; // 2 x y z
        gradPhi[20][2] = 2 * xi[2] * gradPhi[18][2]; // 2 x yy z + 2 xx y z
      }
    }
  }

  void GetWedgePolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) {

    GetWedgePolynomialShapeFunctionGradient(phi,   gradPhi,  xi, solType);

    const unsigned dim = 3;
    const unsigned nDofs = wedgeNumberOfDofs[solType];

    hessPhi.resize(nDofs);

    for(int i = 0; i < nDofs; i++) {
      hessPhi[i].resize(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessPhi[i][i1].assign(dim, 0.);
      }
    }


    if(solType < 1) { //only linear
      //phi_xz
      hessPhi[4][0][2] = hessPhi[4][2][0] = 1.;  //  1
      //phi_yz
      hessPhi[5][2][1] = hessPhi[5][1][2] = 1.;  // 1
    }
    else {  //only quadratic and biquadratic
      //phi_xx
      hessPhi[7][0][0] = 2.; // 2
      //phi_xy
      hessPhi[4][1][0] = hessPhi[4][0][1] = 1.;  //  1
      hessPhi[10][1][0] = hessPhi[10][0][1] = xi[2] ; //  z
      //phi_xz
      hessPhi[5][0][2] = hessPhi[5][2][0] = 1.;  //  1
      hessPhi[10][0][2] = hessPhi[10][2][0] =  xi[1] ; //  y
      //phi_yy
      hessPhi[8][1][1] = 2.; // 2
      //phi_yz
      hessPhi[6][2][1] = hessPhi[6][1][2] = 1.;  // 1
      hessPhi[10][2][1] = hessPhi[10][1][2] = xi[0] ; // x
      //phi_zz
      hessPhi[9][2][2] = 2.; // 2

      if(solType < 2) { //only quadratic
        //phi_xx
        hessPhi[11][0][0] = 2 * xi[2]; // 2 z
        //phi_xz
        hessPhi[11][0][2] = hessPhi[11][2][0] = 2 * xi[0]; // 2 x
        hessPhi[13][0][2] = hessPhi[13][2][0] = 2 * xi[2]; //  2 z
        //phi_yy
        hessPhi[12][1][1] = 2 * xi[2]; // 2 z
        //phi_yz
        hessPhi[12][2][1] = hessPhi[12][1][2] = 2 * xi[1]; // 2 y
        hessPhi[14][2][1] = hessPhi[14][1][2] = 2 * xi[2]; //  2 z;
        //phi_zz
        hessPhi[13][2][2] = 2 * xi[0]; // 2 x
        hessPhi[14][2][2] = 2 * xi[1]; // 2 y ;
      }
      else { //only biquadratic
        //phi_xx
        hessPhi[11][0][0] = 2 * xi[1]; // 2 y
        hessPhi[12][0][0] = 2 * xi[2]; // 2 z
        hessPhi[16][0][0] = 2 * phi[9]; // 2 zz
        hessPhi[18][0][0] = 2 * phi[6]; // 2 y z
        hessPhi[20][0][0] = 2 * phi[15]; //  2 y zz
        //phi_xy
        hessPhi[11][1][0] = hessPhi[11][0][1] = 2 * xi[0] + 2 * xi[1]; // 2 x + 2 y
        hessPhi[18][1][0] = hessPhi[18][0][1] = 2 * (phi[6] + phi[5]); // 2 y z + 2 x z
        hessPhi[19][1][0] = hessPhi[19][0][1] = phi[9]; // zz
        hessPhi[20][1][0] = hessPhi[20][0][1] = hessPhi[18][0][1] * xi[2]; // 2 y zz + 2 x zz
        //phi_xz
        hessPhi[12][0][2] = hessPhi[12][2][0] = 2 * xi[0]; // 2 x
        hessPhi[14][0][2] = hessPhi[14][2][0] = 2 * xi[2]; //  2 z
        hessPhi[16][0][2] = hessPhi[16][2][0] = 4 * phi[5]; // 4 x z
        hessPhi[18][0][2] = hessPhi[18][2][0] = phi[8] + 2 * phi[4]; //  yy  + 2 x y
        hessPhi[19][0][2] = hessPhi[19][2][0] = 2 * phi[6]; // 2 y z
        hessPhi[20][0][2] = hessPhi[20][2][0] = 2 * xi[2] * hessPhi[18][0][3]; //  2 yy z + 4 x y z
        //phi_yy
        hessPhi[11][1][1] = 2 * xi[0] ; // 2 x
        hessPhi[13][1][1] = 2 * xi[2]; // 2 z
        hessPhi[17][1][1] = 2 * phi[9]; // 2 zz
        hessPhi[18][1][1] = 2 * phi[5]; // 2 x z
        hessPhi[20][1][1] = 2 * phi[14]; // 2 x zz
        //phi_yz
        hessPhi[13][2][1] = hessPhi[13][1][2] = 2 * xi[1]; // 2 y
        hessPhi[15][2][1] = hessPhi[15][1][2] = 2 * xi[2]; //  2 z
        hessPhi[17][2][1] = hessPhi[17][1][2] = 4 * phi[6]; // 4 y z
        hessPhi[18][2][1] = hessPhi[18][1][2] = 2 * phi[4] + phi[7]; // 2 x y  + xx
        hessPhi[19][2][1] = hessPhi[19][1][2] = 2 * phi[5]; // 2 x z
        hessPhi[20][2][1] = hessPhi[20][1][2] = 2 * xi[2] * hessPhi[18][1][2]; // 4 x y z + 2 xx z
        //phi_zz
        hessPhi[14][2][2] = 2  * xi[0]; // 2 x
        hessPhi[15][2][2] = 2 * xi[1]; // 2 y
        hessPhi[16][2][2] = 2 * phi[7]; // 2 xx
        hessPhi[17][2][2] = 2 * phi[8]; // 2 yy
        hessPhi[19][2][2] = 2 * phi[4]; // 2 x y
        hessPhi[20][2][2] = xi[1] * (hessPhi[19][2][2] + hessPhi[16][2][2]) ; // 2 x yy  + 2 xx y
      }
    }
  }
//END WEDGE



  bool CheckIfPointIsInsideReferenceDomain(std::vector<double> &xi, const short unsigned &ielType, const double &eps) {
    if(ielType == 0) {
      return CheckIfPointIsInsideReferenceDomainHex(xi, eps);
    }
    else if(ielType == 1) {
      return CheckIfPointIsInsideReferenceDomainTet(xi, eps);
    }
    else if(ielType == 2) {
      return CheckIfPointIsInsideReferenceDomainWedge(xi, eps);
    }
    else if(ielType == 3) {
      return CheckIfPointIsInsideReferenceDomainQuad(xi, eps);
    }
    else if(ielType == 4) {
      return CheckIfPointIsInsideReferenceDomainTri(xi, eps);
    }
    else{// if(ielType == 5) {
      return CheckIfPointIsInsideReferenceDomainLine(xi, eps);
    }
    
  }
  bool CheckIfPointIsInsideReferenceDomainHex(std::vector<double> &xi, const double &eps) {
    double threshold = 1. + eps;
    return (fabs(xi[0]) < threshold && fabs(xi[1]) < threshold && fabs(xi[2]) < threshold) ? true : false;
  }
  bool CheckIfPointIsInsideReferenceDomainTet(std::vector<double> &xi, const double &eps) {
    return (xi[0] > - eps && xi[1] > -eps &&  xi[2] > -eps && xi[0] + xi[1] + xi[2] < 1. + eps) ? true : false;
  }
  bool CheckIfPointIsInsideReferenceDomainWedge(std::vector<double> &xi, const double &eps) {
    double threshold = 1. + eps;
    return (xi[0] > -eps && xi[1] > -eps && xi[0] + xi[1] < threshold && fabs(xi[2]) < threshold) ? true : false;
  }
  bool CheckIfPointIsInsideReferenceDomainQuad(std::vector<double> &xi, const double &eps) {
    double threshold = 1. + eps;
    return (fabs(xi[0]) < threshold && fabs(xi[1]) < threshold) ? true : false;
  }
  bool CheckIfPointIsInsideReferenceDomainTri(std::vector<double> &xi, const double &eps) {
    return (xi[0] > -eps && xi[1] > -eps && xi[0] + xi[1] < 1. + eps) ? true : false;
  }
  bool CheckIfPointIsInsideReferenceDomainLine(std::vector<double> &xi, const double &eps) {
    double threshold = 1. + eps;
    return (fabs(xi[0]) < threshold) ? true : false;
  }


  bool SPDCheck2D(const std::vector< std::vector <double> > &A) {
    bool SPD = true;

    if(A[0][1] != A[1][0]) {
      SPD = false;
      std::cout << "The 2D matrix is not symmetric" << std::endl;
    }

    else if(A[0][0] < 1.0e-8) {
      SPD = false;
    }
    else {
      double lambda = A[1][1] - ((A[0][1] * A[0][1]) / A[0][0]);
      if(lambda < 0 || fabs(lambda) < 1.0e-8) {
        SPD = false;
      }
    }

    if(SPD == true) {
      std::cout << "The 2D matrix is SPD" << std::endl;
    }

    else {
      std::cout << "The 2D matrix is not SPD" << std::endl;
    }

    return SPD;

  }




  bool SPDCheck3D(const std::vector< std::vector <double> > &A) {

    bool SPD = true;
    bool notSymm = false;

    if(A[0][2] != A[2][0] || A[1][2] != A[2][1] || A[0][1] != A[1][0]) {
      notSymm = true;
      std::cout << "The 3D matrix is not symmetric" << std::endl;
    }

    std::vector < std::vector <double> > B(2);
    for(int i = 0; i < 2; i++) {
      B[i].resize(2);
    }

    B[0][0] = A[0][0];
    B[0][1] = A[0][1];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1];



    if(SPDCheck2D(B) == false || notSymm == true) {
      SPD = false;
    }

    else {
      double l00 = sqrt(A[0][0]);
      double l11 = sqrt(A[1][1] - ((A[0][1] * A[0][1]) / A[0][0]));
      double l01 = A[0][1] / l00;

      double detL =  l00 * l11;
      std::vector < std::vector < double > > Lm1(2);
      for(int i = 0; i < 2; i++) {
        Lm1[i].resize(2);
      }

      Lm1[0][0] = (1 / detL) * l11;
      Lm1[1][0] = -(1 / detL) * l01;
      Lm1[0][1] = 0. ;
      Lm1[1][1] = (1 / detL) * l00 ;

      std::vector < double > K(2);

      K[0] = Lm1[0][0] * A[0][2] + Lm1[0][1] * A[1][2];
      K[1] = Lm1[1][0] * A[0][2] + Lm1[1][1] * A[1][2];

      double KK = sqrt(K[0] * K[0] + K[1] * K[1]) ;

      if(A[2][2] - KK < 0 || fabs(A[2][2] - KK) < 1.0e-8) {
        SPD = false;
      }
    }

    if(SPD == true) {
      std::cout << "The 3D matrix is SPD" << std::endl;
    }

    else {
      std::cout << "The 3D matrix is not SPD" << std::endl;
    }

    return SPD;

  }



  bool GetNewLocalCoordinates(std::vector <double> &xi, const std::vector< double > &x, const std::vector <double> &phi,
                              const std::vector < std::vector <double > > &gradPhi,
                              const std::vector < std::vector <double > > &a) {

    const unsigned dim = gradPhi[0].size();
    const unsigned  nDofs = phi.size();

    bool convergence = false;
    std::vector < double > F(dim, 0.);
    std::vector < std::vector < double > > J(dim);

    for(int k = 0; k < dim; k++) {
      J[k].assign(dim, 0.);
    }

    for(int k = 0; k < dim; k++) {
      for(int i = 0; i < nDofs; i++) {
        F[k] += a[k][i] * phi[i];

        for(int i1 = 0; i1 < dim; i1++) {
          J[k][i1] += a[k][i] * gradPhi[i][i1];
        }
      }
      F[k] -= x[k];
    }


    std::vector < std::vector < double > >  Jm1;
    InverseMatrix(J, Jm1);

    double delta2 = 0.;

    for(int i1 = 0; i1 < dim; i1++) {
      double deltak = 0.;

      for(int i2 = 0; i2 < dim; i2++) {
        deltak -= Jm1[i1][i2] * F[i2];
      }

      xi[i1] += deltak;
      delta2 += deltak * deltak;
    }

    if(delta2 < 1.0e-9) {
      convergence = true;
    }

    return convergence;
  }


  bool GetNewLocalCoordinatesHess(std::vector <double> &xi, const std::vector< double > &x, const std::vector <double> &phi,
                                  const std::vector < std::vector <double > > &gradPhi, const std::vector < std::vector < std::vector <double> > > hessPhi,
                                  const std::vector < std::vector <double > > &a) {

    const unsigned dim = gradPhi[0].size();
    const unsigned  nDofs = phi.size();

    bool convergence = false;
    std::vector < double > xp(dim, 0.);
    std::vector < std::vector < double > > gradXp(dim);
    std::vector < std::vector < std::vector < double > > > hessXp(dim);

    for(int k = 0; k < dim; k++) {
      gradXp[k].assign(dim, 0.);
      hessXp[k].resize(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessXp[k][i1].assign(dim, 0.);
      }
    }

    for(int k = 0; k < dim; k++) {
      for(int i = 0; i < nDofs; i++) {
        xp[k] += a[k][i] * phi[i];

        for(int i1 = 0; i1 < dim; i1++) {
          gradXp[k][i1] += a[k][i] * gradPhi[i][i1];

          for(int i2 = 0; i2 < dim; i2++) {
            hessXp[k][i1][i2] += a[k][i] * hessPhi[i][i1][i2];
          }
        }
      }
    }

    std::vector < double > gradF(dim, 0.);
    std::vector < std::vector < double > >  hessF(dim);

    for(int i1 = 0; i1 < dim; i1++) {
      hessF[i1].assign(dim, 0.);
    }

    for(int k = 0; k < dim; k++) {
      for(int i1 = 0; i1 < dim; i1++) {
        gradF[i1] += -2. * (x[k] - xp[k]) * gradXp[k][i1];

        for(int i2 = 0; i2 < dim; i2++) {
          hessF[i1][i2] += -2. * (x[k] - xp[k]) * hessXp[k][i1][i2] + 2. * gradXp[k][i1] * gradXp[k][i2];
        }
      }
    }

    std::vector < std::vector < double > >  hessFm1;
    InverseMatrix(hessF, hessFm1);

    double delta2 = 0.;

    for(int i1 = 0; i1 < dim; i1++) {
      double deltak = 0.;

      for(int i2 = 0; i2 < dim; i2++) {
        deltak -= hessFm1[i1][i2] * gradF[i2];
      }

      xi[i1] += deltak;
      delta2 += deltak * deltak;
    }

    if(delta2 < 1.0e-9) {
      convergence = true;
    }

    return convergence;
  }



  void InverseMatrix(const std::vector< std::vector <double> > &A, std::vector< std::vector <double> > &invA) {

    unsigned dim = A.size();
    invA.resize(dim);

    for(int i = 0; i < dim; i++) {
      invA[i].resize(dim);
    }

    double detA;
    if(dim == 1) {

      if(A[0][0] == 0) {
        std::cout << " ERROR: the matrix is singular " << std::endl;
        abort();
      }
      invA[0][0] = 1. / A[0][0];
    }
    else if(dim == 2) {

      detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

      if(detA == 0) {
        std::cout << " ERROR: the matrix is singular " << std::endl;
        abort();
      }
      else {
        invA[0][0] = A[1][1] / detA;
        invA[0][1] = -A[0][1] / detA;
        invA[1][0] = -A[1][0] / detA;
        invA[1][1] = A[0][0] / detA;
      }
    }
    else if(dim == 3) {

      detA = (A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1])
             - (A[2][0] * A[1][1] * A[0][2] + A[2][1] * A[1][2] * A[0][0] + A[2][2] * A[1][0] * A[0][1]) ;

      if(detA == 0) {
        std::cout << " ERROR: the matrix is singular " << std::endl;
        abort();
      }
      else {

        invA[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) / detA ;
        invA[0][1] = (A[0][2] * A[2][1] - A[2][2] * A[0][1]) / detA ;
        invA[0][2] = (A[0][1] * A[1][2] - A[1][1] * A[0][2]) / detA ;
        invA[1][0] = (A[1][2] * A[2][0] - A[2][2] * A[1][0]) / detA ;
        invA[1][1] = (A[0][0] * A[2][2] - A[2][0] * A[0][2]) / detA ;
        invA[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / detA ;
        invA[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) / detA ;
        invA[2][1] = (A[0][1] * A[2][0] - A[2][1] * A[0][0]) / detA ;
        invA[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) / detA ;

      }
    }
    else {
      std::cout << " ERROR: the matrix is neither 2x2 nor 3x3 so we cannot use this function " << std::endl;
      abort();
    }
  }

  void GetConvexHullSphere(const std::vector< std::vector < double > > &xv, std::vector <double> &xc, double & r, const double tolerance) {
    unsigned dim = xv.size();
    unsigned ndofs = xv[0].size();
    xc.resize(dim, 0.);
    for(int d = 0; d < dim; d++) {
      for(int i = 0; i < ndofs; i++) {
        xc[d] += xv[d][i];
      }
      xc[d] /= ndofs;
    }
    double r2 = 0.;
    for(unsigned j = 0; j < ndofs; j++) {
      double d2 = 0.;
      for(int d = 0; d < dim; d++) {
        d2 += (xv[d][j] - xc[d]) * (xv[d][j] - xc[d]);
      }
      r2 = (r2 > d2) ? r2 : d2;
    }

    r = (1. + tolerance) * sqrt(r2);
  }

  void GetBoundingBox(const std::vector< std::vector < double > > &xv, std::vector< std::vector < double > > &xe, const double tolerance) {
    unsigned dim = xv.size();
    unsigned ndofs = xv[0].size();
    xe.resize(dim);
    for(int d = 0; d < dim; d++) {
      xe[d].resize(2);
      xe[d][0] = xv[d][0];
      xe[d][1] = xv[d][0];
    }
    for(int d = 0; d < dim; d++) {
      for(int i = 1; i < ndofs; i++) {
        xe[d][0] = (xv[d][i] < xe[d][0]) ? xv[d][i] : xe[d][0];
        xe[d][1] = (xv[d][i] > xe[d][1]) ? xv[d][i] : xe[d][1];
      }
    }
    for(int d = 0; d < dim; d++) {
      double epsilon = tolerance * (xe[d][1] - xe[d][0]);
      xe[d][0] -= epsilon;
      xe[d][1] += epsilon;
    }
  }

  void GetInverseMapping(const unsigned &solType, short unsigned &ielType, const std::vector < std::vector < std::vector <double > > > &aP,
                         const std::vector <double > &xl, std::vector <double > &xi) {

    for(short unsigned jtype = 0; jtype < solType + 1; jtype++) {
      std::vector < double > phi;
      std::vector < std::vector < double > > gradPhi;
      bool convergence = false;
      while(!convergence) {
        GetPolynomialShapeFunctionGradient(phi, gradPhi, xi, ielType, jtype);
        convergence = GetNewLocalCoordinates(xi, xl, phi, gradPhi, aP[jtype]);
      }
    }
  }

  const double XI[6][27][3] = {{
      { -1, -1, -1}, {1, -1, -1}, {1, 1, -1}, { -1, 1, -1}, { -1, -1, 1}, {1, -1, 1}, {1, 1, 1}, { -1, 1, 1}, {0, -1, -1},
      {1, 0, -1}, {0, 1, -1}, { -1, 0, -1}, {0, -1, 1}, {1, 0, 1}, {0, 1, 1}, { -1, 0, 1}, { -1, -1, 0},
      {1, -1, 0}, {1, 1, 0}, { -1, 1, 0}, {0, -1, 0}, {1, 0, 0}, {0, 1, 0}, { -1, 0, 0}, {0, 0, -1},
      {0, 0, 1}, {0, 0, 0}
    },
    {
      {0, 0, 0},        {1, 0, 0},       {0, 1, 0},   {0, 0, 1},         		//0->4
      {0.5, 0, 0},      {0.5, 0.5, 0},   {0, 0.5, 0},
      {0.,  0, 0.5},    {0.5, 0., 0.5},  {0, 0.5, 0.5},                  		//5->9
      {1. / 3., 1. / 3., 0.}, {1. / 3., 0., 1. / 3.}, {1. / 3., 1. / 3., 1. / 3.}, {0., 1. / 3., 1. / 3.}, //external faces of internal tetrahedra
      {0.25, 0.25, 0.25} //34
    },
    {
      {0, 0, -1},      {1, 0, -1},      {0, 1, -1},                     //vertici triangoli
      {0, 0, 1},       {1, 0, 1},       {0, 1, 1},
      {0.5, 0, -1},    {0.5, 0.5, -1},  {0, 0.5, -1},                   //midpoints triangoli
      {0.5, 0, 1},     {0.5, 0.5, 1},   {0, 0.5, 1},
      {0, 0, 0},       {1, 0, 0} ,      {0, 1, 0},                      //midpoints quadrati
      {0.5, 0, 0},     {0.5, 0.5, 0},    {0, 0.5, 0}, //0->17           //facce quadrati
      {1. / 3., 1. / 3., -1}, {1. / 3., 1. / 3., 1}, {1. / 3., 1. / 3., 0}
    },
    {
      { -1, -1}, {1, -1}, {1, 1}, { -1, 1},
      { 0, -1}, {1, 0}, {0, 1}, { -1, 0}, {0, 0}
    },
    {
      {0, 0},         {1, 0},         {0, 1},
      {0.5, 0},       {0.5, 0.5},     {0, 0.5},
      {1. / 3., 1. / 3.}
    },
    {
      { -1}, {1}, {0}
    }
  };

  void GetClosestPointInReferenceElement(const std::vector< std::vector < double > > &xv, const std::vector <double> &x,
                                         const short unsigned &ieltype, std::vector < double > &xi) {
    unsigned dim = xv.size();
    unsigned ndofs = xv[0].size();
    unsigned jmin = ndofs;
    double d2min = 1.0e100;
    for(unsigned j = 0; j < ndofs; j++) {
      double d2 = 0;
      for(int d = 0; d < dim; d++) {
        d2 += (xv[d][j] - x[d]) * (xv[d][j] - x[d]);
      }
      if(d2 < d2min) {
        d2min = d2;
        jmin = j;
      }
    }

    xi.resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      xi[k] = XI[ieltype][jmin][k];
    }
  }

  void PrintLine(const std::string output_path, const std::vector < std::vector< std::vector<double> > > &xn, const bool &streamline, const unsigned &step) {

    // *********** open vtu files *************
    std::ofstream fout;

    std::string dirnamePVTK = "";
    Files files;
    files.CheckDir(output_path, dirnamePVTK);
    
    
    

    std::string filename_prefix = (streamline) ? "streamline" : "line";

    std::ostringstream filename;
    filename << output_path << "/./" << filename_prefix << "." << step << ".vtu";

    fout.open(filename.str().c_str());
    
    if(!fout.is_open()) {
      std::cout << std::endl << " The output file " << filename.str() << " cannot be opened.\n";
      abort();
    }

    unsigned nvt = 0;
    unsigned nel = 0;

    for(unsigned l = 0; l < xn.size(); l++) {
      nvt += xn[l].size();
      nel += xn[l].size() - 1;
    }

    const unsigned dim_array_coord [] = { nvt * 3 * static_cast< unsigned int >( sizeof(float) )};
    const unsigned dim_array_conn[]   = { nel * 2 * static_cast< unsigned int >( sizeof(int) ) };
    const unsigned dim_array_off []   = { nel * static_cast< unsigned int >( sizeof(int) ) };
    const unsigned dim_array_type []  = { nel * static_cast< unsigned int >( sizeof(short unsigned) ) };

    unsigned buffer_size = (dim_array_coord[0] > dim_array_conn[0]) ? dim_array_coord[0] : dim_array_conn[0];
    void* buffer_void = new char [buffer_size];
    char* buffer_char = static_cast <char*>(buffer_void);

    size_t cch;
    cch = b64::b64_encode(&buffer_char[0], buffer_size , NULL, 0);
    std::vector <char> enc;
    enc.resize(cch);
    char* pt_char;

    // *********** write vtu header ************
    fout << "<?xml version=\"1.0\"?>" << std::endl;
    fout << "<VTKFile type = \"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    fout << "  <UnstructuredGrid>" << std::endl;


    fout  << "    <Piece NumberOfPoints= \"" << nvt
          << "\" NumberOfCells= \"" << nel
          << "\" >" << std::endl;

    //-----------------------------------------------------------------------------------------------
    // print coordinates *********************************************Solu*******************************************
    fout  << "      <Points>" << std::endl;
    fout  << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;

    // point pointer to common mamory area buffer of void type;
    float* var_coord = static_cast< float* >(buffer_void);

    unsigned counter = 0;
    for(unsigned l = 0; l < xn.size(); l++) {
      for(unsigned i = 0; i < xn[l].size(); i++) {
        for(unsigned d = 0; d < 3; d++) {
          var_coord[counter * 3 + d] = (xn[l][i].size() > d) ? xn[l][i][d] : 0.;
        }
        counter++;
      }
    }

    cch = b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), NULL, 0);
    b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), &enc[0], cch);
    pt_char = &enc[0];
    for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;

    //print coordinates array
    cch = b64::b64_encode(&var_coord[0], dim_array_coord[0] , NULL, 0);
    b64::b64_encode(&var_coord[0], dim_array_coord[0], &enc[0], cch);
    pt_char = &enc[0];
    for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;
    fout << std::endl;

    fout  << "        </DataArray>" << std::endl;
    fout  << "      </Points>" << std::endl;
    //-----------------------------------------------------------------------------------------------

    // Printing of element connectivity - offset - format type  *
    fout  << "      <Cells>" << std::endl;
    //-----------------------------------------------------------------------------------------------
    //print connectivity
    fout  << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">" << std::endl;

    // point pointer to common mamory area buffer of void type;
    int* var_conn = static_cast <int*>(buffer_void);
    unsigned icount = 0;
    unsigned skip = 0;
    unsigned restart = 0;
    for(unsigned iel = 0; iel < nel; iel++) {
      for(unsigned j = 0; j < 2; j++) {
        var_conn[icount] = iel + j + skip;
        icount++;
      }
      restart++;
      if(restart == xn[skip].size() - 1) {
        skip++;
        restart = 0;
      }
    }

    //print connectivity dimension
    cch = b64::b64_encode(&dim_array_conn[0], sizeof(dim_array_conn), NULL, 0);
    b64::b64_encode(&dim_array_conn[0], sizeof(dim_array_conn), &enc[0], cch);
    pt_char = &enc[0];
    for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;

    //print connectivity array
    cch = b64::b64_encode(&var_conn[0], dim_array_conn[0] , NULL, 0);
    b64::b64_encode(&var_conn[0], dim_array_conn[0], &enc[0], cch);
    pt_char = &enc[0];
    for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    //------------------------------------------------------------------------------------------------
    fout  << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">" << std::endl;
    // point pointer to common memory area buffer of void type;
    int* var_off = static_cast <int*>(buffer_void);
    // print offset array
    for(int iel = 0; iel < nel; iel++) {
      var_off[iel] = (iel + 1) * 2;
    }

    //print offset dimension
    cch = b64::b64_encode(&dim_array_off[0], sizeof(dim_array_off), NULL, 0);
    b64::b64_encode(&dim_array_off[0], sizeof(dim_array_off), &enc[0], cch);
    pt_char = &enc[0];
    for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;

    //print offset array
    cch = b64::b64_encode(&var_off[0], dim_array_off[0] , NULL, 0);
    b64::b64_encode(&var_off[0], dim_array_off[0], &enc[0], cch);
    pt_char = &enc[0];
    for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;

    fout  << std::endl;

    fout  << "        </DataArray>" << std::endl;



    //--------------------------------------------------------------------------------------------------

    //Element format type : 23:Serendipity(8-nodes)  28:Quad9-Biquadratic
    fout  << "        <DataArray type=\"UInt16\" Name=\"types\" format=\"binary\">" << std::endl;

    // point pointer to common mamory area buffer of void type;
    unsigned short* var_type = static_cast <unsigned short*>(buffer_void);

    for(unsigned iel = 0; iel < nel; iel++) {
      var_type[iel] = 3;
    }

    //print element format dimension
    cch = b64::b64_encode(&dim_array_type[0], sizeof(dim_array_type), NULL, 0);
    b64::b64_encode(&dim_array_type[0], sizeof(dim_array_type), &enc[0], cch);
    pt_char = &enc[0];
    for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;

    //print element format array
    cch = b64::b64_encode(&var_type[0], dim_array_type[0] , NULL, 0);
    b64::b64_encode(&var_type[0], dim_array_type[0], &enc[0], cch);
    pt_char = &enc[0];
    for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;

    fout  << std::endl;
    fout  << "        </DataArray>" << std::endl;



    //----------------------------------------------------------------------------------------------------
//
    fout  << "      </Cells>" << std::endl;



    //-----------------------------------------------------------------------------------------------
    // Printing of element connectivity - offset - format type  *
//   fout  << "      <Cells>" << std::endl;
//   fout  << "      </Cells>" << std::endl;

    fout << "    </Piece>" << std::endl;
    fout << "  </UnstructuredGrid>" << std::endl;
    fout << "</VTKFile>" << std::endl;
    fout.close();

    delete [] var_coord;
    //--------------------------------------------------------------------------------------------------------
    return;
  }



}
