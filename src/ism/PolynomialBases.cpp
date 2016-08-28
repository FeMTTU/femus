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

#include "PolynomialBases.hpp"
#include "GeomElTypeEnum.hpp"

namespace femus {

  //BEGIN Interface
  void ProjectNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector< std::vector < double > > &aN,
      const short unsigned &ielType, const unsigned &solType) {

    if (ielType == QUAD) {
      ProjectQuadNodalToPolynomialCoefficients(aP, aN, solType);
    }
  }

  void GetPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi,
                                  short unsigned &ielType, const unsigned & solType) {

    if (ielType == QUAD) {
      GetQuadPolynomialShapeFunction(phi,  xi,  solType);
    }
  }

  void GetPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
                                          const std::vector < double >& xi,  short unsigned &ielType, const unsigned & solType) {

    if (ielType == QUAD) {
      GetQuadPolynomialShapeFunctionGradient(phi, gradPhi,  xi,  solType);
    }
  }

  void GetPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
      std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi,
      short unsigned &ielType, const unsigned & solType) {

    if (ielType == QUAD) {
      GetQuadPolynomialShapeFunctionGradientHessian(phi, gradPhi,  hessPhi,  xi,  solType);
    }
  }
  //END Interface

  //BEGIN QUAD specialized functions
   const unsigned quadNumberOfDofs[3] = {4, 8, 9};
  
  void ProjectQuadNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solType) {

    unsigned dim =  aN.size();
    aP.resize(dim);
    unsigned nDofs = aN[0].size();
    
    if( nDofs != quadNumberOfDofs[solType]){
      std::cout << "Error in ProjectQuadNodalToPolynomialCoefficients(...) the number of Dofs is inconsistent"<<std::endl;
      abort();
    }

    for (unsigned k = 0; k < dim; k++) {
      aP[k].resize(nDofs);
    }

    if (solType == 0) {
      for (int k = 0; k < dim; k++) {
        aP[k][0] = 0.25 * (aN[k][0] + aN[k][1] + aN[k][2] + aN[k][3]);
        aP[k][1] = 0.25 * (-aN[k][0] + aN[k][1] + aN[k][2] - aN[k][3]);
        aP[k][2] = 0.25 * (-aN[k][0] - aN[k][1] + aN[k][2] + aN[k][3]);
        aP[k][3] = 0.25 * (aN[k][0] - aN[k][1] + aN[k][2] - aN[k][3]);
      }
    }
    else if (solType == 1) {
      for (int k = 0; k < dim; k++) {
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
    else if (solType == 2) {
      for (int k = 0; k < dim; k++) {
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

    if (solType > 0) {
      phi[4] = xi[0] * xi[0]; // x x
      phi[5] = xi[1] * xi[1]; // y y
      phi[6] = phi[4] * xi[1]; // xx y
      phi[7] = phi[5] * xi[0]; // x yy

      if (solType > 1) {
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

    for (int i = 0; i < nDofs; i++) {
      gradPhi[i].assign(dim, 0.);
    }

    //phi_x
    gradPhi[1][0] = 1.; // 1
    gradPhi[3][0] = xi[1];  //  y
    //phi_y
    gradPhi[2][1] = 1.;  // 1
    gradPhi[3][1] = xi[0];  // x

    if (solType > 0) {
      //phi_x
      gradPhi[4][0] = 2. * xi[0]; // 2 x
      gradPhi[6][0] = 2. * phi[3]; // 2 x y
      gradPhi[7][0] = phi[5]; //  yy
      //phi_y
      gradPhi[5][1] = 2. * xi[1]; // 2 y
      gradPhi[6][1] = phi[4]; // xx
      gradPhi[7][1] = 2. * phi[3]; // 2 x y

      if (solType > 1) {
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

    for (int i = 0; i < nDofs; i++) {
      hessPhi[i].resize(dim);

      for (int i1 = 0; i1 < dim; i1++) {
        hessPhi[i][i1].assign(dim, 0.);
      }
    }

    //phi_xy
    hessPhi[3][1][0] = hessPhi[3][0][1] = 1.; //1

    if (solType > 0) {
      //phi_xx
      hessPhi[4][0][0] = 2.; // 2
      hessPhi[6][0][0] = 2. * xi[1]; // 2 y
      //phi_xy
      hessPhi[6][1][0] = hessPhi[6][0][1] = 2. * xi[0]; // 2 x
      hessPhi[7][1][0] = hessPhi[7][0][1] = 2. * xi[1]; // 2 y
      //phi_yy
      hessPhi[5][1][1] = 2.; // 2
      hessPhi[7][1][1] = 2. * xi[0]; // 2 x

      if (solType > 1) {
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

}
