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

    else if (ielType == TRI) {
        ProjectTriNodalToPolynomialCoefficients(aP, aN, solType);
    }
}

void GetPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi,
                                short unsigned &ielType, const unsigned & solType) {

    if (ielType == QUAD) {
        GetQuadPolynomialShapeFunction(phi,  xi,  solType);
    }

    else if (ielType == TRI) {
        GetTriPolynomialShapeFunction(phi, xi, solType);
    }
}

void GetPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
                                        const std::vector < double >& xi,  short unsigned &ielType, const unsigned & solType) {

    if (ielType == QUAD) {
        GetQuadPolynomialShapeFunctionGradient(phi, gradPhi,  xi,  solType);
    }

    else if(ielType == TRI) {
        GetTriPolynomialShapeFunctionGradient(phi, gradPhi, xi, solType);
    }
}

void GetPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi,
        std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi,
        short unsigned &ielType, const unsigned & solType) {

    if (ielType == QUAD) {
        GetQuadPolynomialShapeFunctionGradientHessian(phi, gradPhi,  hessPhi,  xi,  solType);
    }

    else if(ielType == TRI) {
        GetTriPolynomialShapeFunctionGradientHessian(phi, gradPhi, hessPhi, xi, solType);
    }
}
//END Interface

//BEGIN QUAD specialized functions
const unsigned quadNumberOfDofs[3] = {4, 8, 9};

void ProjectQuadNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solType) {

    unsigned dim =  aN.size();
    aP.resize(dim);
    unsigned nDofs = aN[0].size();

    if( nDofs != quadNumberOfDofs[solType]) {
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

//BEGIN TRI specialized functions
const unsigned triNumberOfDofs[3] = {3, 6, 7};

void ProjectTriNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solType) {

    unsigned dim =  aN.size();
    aP.resize(dim);
    unsigned nDofs = aN[0].size();

    if( nDofs != triNumberOfDofs[solType]) {
        std::cout << "Error in ProjectTriNodalToPolynomialCoefficients(...) the number of Dofs is inconsistent"<<std::endl;
        abort();
    }

    for (unsigned k = 0; k < dim; k++) {
        aP[k].resize(nDofs);
    }

    if (solType == 0) {
        for (int k = 0; k < dim; k++) {
            aP[k][0] = aN[k][0];
            aP[k][1] = - aN[k][0] + aN[k][1];
            aP[k][2] = - aN[k][0] + aN[k][2];
        }
    }
    else if (solType == 1) {
        for (int k = 0; k < dim; k++) {
            aP[k][0] = aN[k][0];
            aP[k][1] = - 3 * aN[k][0] - aN[k][1] + 4 * aN[k][3];
            aP[k][2] = - 3 * aN[k][0] - aN[k][2] + 4 * aN[k][5];
            aP[k][3] = 4 * aN[k][0] - 4 * aN[k][3] + 4 * aN[k][4] - 4 * aN[k][5];
            aP[k][4] = 2 * aN[k][0] + 2 * aN[k][1] - 4 * aN[k][3];
            aP[k][5] = 2 * aN[k][0] + 2 * aN[k][2] - 4 * aN[k][5];
        }
    }
    else if (solType == 2) {
        for (int k = 0; k < dim; k++) {
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

    if (solType > 0) {
        phi[3] = xi[0] * xi[1]; // x y
        phi[4] = xi[0] * xi[0];  // x x
        phi[5] = xi[1] * xi[1]; // y y

        if (solType > 1) {
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

    for (int i = 0; i < nDofs; i++) {
        gradPhi[i].assign(dim, 0.);
    }

    //phi_x
    gradPhi[1][0] = 1.; // 1
    //phi_y
    gradPhi[2][1] = 1.;  // 1

    if (solType > 0) {
        //phi_x
        gradPhi[3][0] = xi[1]; // y
        gradPhi[4][0] = 2.*xi[0] ;  // 2 x
        //phi_y
        gradPhi[3][1] = xi[0]; // x
        gradPhi[5][1] = 2.*xi[1]; // 2*y


        if (solType > 1) {
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

    for (int i = 0; i < nDofs; i++) {
        hessPhi[i].resize(dim);

        for (int i1 = 0; i1 < dim; i1++) {
            hessPhi[i][i1].assign(dim, 0.);
        }
    }

    if (solType > 0) {
        //phi_xx
        hessPhi[4][0][0] = 2.; // 2
        //phi_xy
        hessPhi[3][1][0] = hessPhi[3][0][1] = 1.; // 1
        //phi_yy
        hessPhi[5][1][1] = 2.; // 2

        if (solType > 1) {
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
    unsigned nDofs = aN[0].size();

    if( nDofs != hexNumberOfDofs[solType]) {
        std::cout << "Error in ProjectHexNodalToPolynomialCoefficients(...) the number of Dofs is inconsistent"<<std::endl;
        abort();
    }

    for (unsigned k = 0; k < dim; k++) {
        aP[k].resize(nDofs);
    }

    if (solType == 0) {
        for (int k = 0; k < dim; k++) {
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
    else if (solType == 1) {
        for (int k = 0; k < dim; k++) {
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
    else if (solType == 2) {
        for (int k = 0; k < dim; k++) {
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
//END HEX
}
