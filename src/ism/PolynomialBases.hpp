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

#ifndef __femus_ism_PolynomialBases_hpp__
#define __femus_ism_PolynomialBases_hpp__

#include "vector"

namespace femus {

  // interface
  void ProjectNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector< std::vector < double > > &aN, const short unsigned &ielType, const unsigned &solType) ;
  void GetPolynomialShapeFunction(std::vector < double >& phi,  const std::vector < double >& xi, short unsigned &ielType, const unsigned & solType) ;
  void GetPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, short unsigned &ielType, const unsigned & solType) ;
  void GetPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, short unsigned &ielType, const unsigned & solType) ;

  // QUAD specialized functions
  void ProjectQuadNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetQuadPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetQuadPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetQuadPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;

  // TRI specialized functions
  void ProjectTriNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetTriPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetTriPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetTriPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;
  
  // HEX specialized functions
  void ProjectHexNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetHexPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetHexPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetHexPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;
  
  // TET specialized functions
  void ProjectTetNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetTetPolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetTetPolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetTetPolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;
  
  // WEDGE specialized functions
  void ProjectWedgeNodalToPolynomialCoefficients(std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetWedgePolynomialShapeFunction(std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetWedgePolynomialShapeFunctionGradient(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetWedgePolynomialShapeFunctionGradientHessian(std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;
  
  
}
#endif
