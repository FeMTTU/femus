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

#include <vector>
#include "Files.hpp"
#include <b64/b64.h>

namespace femus {
  // interface
  void ProjectNodalToPolynomialCoefficients (std::vector < std::vector <double > > &aP, const std::vector< std::vector < double > > &aN, const short unsigned &ielType, const unsigned &solType) ;
  void InterpolatePolynomialCoefficients (std::vector<std::vector < std::vector <double > > > &aXs, const std::vector<std::vector < std::vector <double > > > &aX0,
                                          const std::vector<std::vector < std::vector <double > > >aX1, const double &s);
  void GetPolynomialShapeFunction (std::vector < double >& phi,  const std::vector < double >& xi, short unsigned &ielType, const unsigned & solType) ;
  void GetPolynomialShapeFunctionGradient (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, short unsigned &ielType, const unsigned & solType) ;
  void GetPolynomialShapeFunctionGradientHessian (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, short unsigned &ielType, const unsigned & solType) ;

  // LINE specialized functions
  void ProjectLineNodalToPolynomialCoefficients (std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetLinePolynomialShapeFunction (std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetLinePolynomialShapeFunctionGradient (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetLinePolynomialShapeFunctionGradientHessian (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;
    
  // QUAD specialized functions
  void ProjectQuadNodalToPolynomialCoefficients (std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetQuadPolynomialShapeFunction (std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetQuadPolynomialShapeFunctionGradient (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetQuadPolynomialShapeFunctionGradientHessian (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;

  // TRI specialized functions
  void ProjectTriNodalToPolynomialCoefficients (std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetTriPolynomialShapeFunction (std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetTriPolynomialShapeFunctionGradient (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetTriPolynomialShapeFunctionGradientHessian (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;

  // HEX specialized functions
  void ProjectHexNodalToPolynomialCoefficients (std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetHexPolynomialShapeFunction (std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetHexPolynomialShapeFunctionGradient (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetHexPolynomialShapeFunctionGradientHessian (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;

  // TET specialized functions
  void ProjectTetNodalToPolynomialCoefficients (std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetTetPolynomialShapeFunction (std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetTetPolynomialShapeFunctionGradient (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetTetPolynomialShapeFunctionGradientHessian (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;

  // WEDGE specialized functions
  void ProjectWedgeNodalToPolynomialCoefficients (std::vector < std::vector <double > > &aP, const std::vector < std::vector <double > > &aN, const unsigned &solutionType) ;
  void GetWedgePolynomialShapeFunction (std::vector < double >& phi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetWedgePolynomialShapeFunctionGradient (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, const std::vector < double >& xi, const unsigned & solType) ;
  void GetWedgePolynomialShapeFunctionGradientHessian (std::vector < double >& phi, std::vector < std::vector < double > >& gradPhi, std::vector < std::vector < std::vector < double > > >& hessPhi, const std::vector < double >& xi, const unsigned & solType) ;

  bool CheckIfPointIsInsideReferenceDomain (std::vector<double> &xi, const short unsigned &ielType, const double &eps = 0.);
  bool CheckIfPointIsInsideReferenceDomainHex (std::vector<double> &xi, const double &eps = 0.);
  bool CheckIfPointIsInsideReferenceDomainTet (std::vector<double> &xi, const double &eps = 0.);
  bool CheckIfPointIsInsideReferenceDomainWedge (std::vector<double> &xi, const double &eps = 0.);
  bool CheckIfPointIsInsideReferenceDomainQuad (std::vector<double> &xi, const double &eps = 0.);
  bool CheckIfPointIsInsideReferenceDomainTri (std::vector<double> &xi, const double &eps = 0.);
  bool CheckIfPointIsInsideReferenceDomainLine (std::vector<double> &xi, const double &eps = 0.);

  bool SPDCheck2D (const std::vector< std::vector <double> > &A);
  bool SPDCheck3D (const std::vector< std::vector <double> > &A);

  bool GetNewLocalCoordinates (std::vector <double> &xi, const std::vector< double > &x, const std::vector <double> &phi,
                               const std::vector < std::vector <double > > &gradPhi,
                               const std::vector < std::vector <double > > &a);

  bool GetNewLocalCoordinatesHess (std::vector <double> &xi, const std::vector< double > &x, const std::vector <double> &phi,
                                   const std::vector < std::vector <double > > &gradPhi, const std::vector < std::vector < std::vector <double> > > hessPhi,
                                   const std::vector < std::vector <double > > &a);


  void InverseMatrix (const std::vector< std::vector <double> > &A, std::vector< std::vector <double> > &invA);

  void GetConvexHullSphere (const std::vector< std::vector < double > > &xv, std::vector <double> &xc, double & r, const double tolerance = 1.0e-10);
  void GetBoundingBox (const std::vector< std::vector < double > > &xv, std::vector< std::vector < double > > &xe, const double tolerance = 1.0e-10);
  void GetInverseMapping (const unsigned &solType, short unsigned &ielType, const std::vector < std::vector < std::vector <double > > > &aP,
                          const std::vector <double > &xl, std::vector <double > &xi);
  bool GetInverseMapping (const unsigned &solType, short unsigned &ielType, const std::vector < std::vector < std::vector <double > > > &aP,
                          const std::vector <double > &xl, std::vector <double > &xi, const unsigned &MaxNumberOfIteration);
  void GetClosestPointInReferenceElement (const std::vector< std::vector < double > > &xv, const std::vector <double> &x,
                                          const short unsigned &ieltype, std::vector < double > &xi);

  void PrintLine (const std::string output_path, const std::string filename, const std::vector < std::vector< std::vector<double> > > &xn, const unsigned &step);


}
#endif

