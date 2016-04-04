/*=========================================================================

 Program: FEMUS
 Module: Mesh
 Authors: Eugenio Aulisa, Giacomo Capodaglio

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Marker.hpp"
#include "NumericVector.hpp"

const unsigned facePointNumber[6] = {0, 0, 0, 9, 7, 0};


const unsigned facePoints[6][9] = {
  { },
  { },
  { },
  { 0, 4, 1, 5, 2, 6, 3, 7, 0},
  { 0, 3, 1, 4, 2, 5, 0},
  { }
};

namespace femus {

  void Marker::GetElement() {

    unsigned dim = _mesh->GetDimension();
    std::vector< std::vector < double > > xv(dim);
    for(unsigned k = 0; k < dim; k++) {
      xv[k].reserve(9);
    }

    bool ElementHasBeenFound = false;
    for(int iel = _mesh->_elementOffset[_iproc]; iel < _mesh->_elementOffset[_iproc + 1]; iel++) {
      short unsigned ielType = _mesh->GetElementType(iel);
      for(unsigned k = 0; k < dim; k++) {
        xv[k].resize(facePointNumber[ielType]);
      }
      for(unsigned i = 0; i < facePointNumber[ielType]; i++) {
        unsigned iDof  = _mesh->GetSolutionDof(facePoints[ielType][i], iel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
          xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof) - _x[k];     // global extraction and local storage for the element coordinates
        }
      }
      double w = GetWindingNumber(xv, iel);
      if(w > 0) {
        _elem = iel;
        ElementHasBeenFound = true;
        break;
      }
      else if(w < 0) {
        std::cout << "Error negative Winding Number with counterclockwise oriented points" << std::endl;
        abort();
      }

    }
    if(ElementHasBeenFound)
      std::cout << "The marker belongs to element " << _elem << std::endl;
    else {
      std::cout << " The marker does not belong to this portion of the mesh" << std::endl;
    }
  }

  double Marker::GetWindingNumber(std::vector< std::vector < double > > &xv, const int &iel) {

    double length = 0.;
    for(unsigned i = 0; i < xv[0].size() - 1; i++) {
      length += sqrt((xv[0][i + 1] - xv[0][i]) * (xv[0][i + 1] - xv[0][i]) +
                     (xv[1][i + 1] - xv[1][i]) * (xv[1][i + 1] - xv[1][i]));
    }
    
    length /= xv[0].size();
    
    for(unsigned i = 0; i < xv[0].size(); i++) {
      xv[0][i] /= length;
      xv[1][i] /= length;
    }
    
    double epsilon  = 10e-10;

    double w = 0.;
    for(unsigned i = 0; i < xv[0].size() - 1; i++) {
      double Delta = -xv[0][i] * (xv[1][i + 1] - xv[1][i]) + xv[1][i] * (xv[0][i + 1] - xv[0][i]);
      if(iel == 60) {
        std::cout << "Delta for element" << iel << " is =" << Delta  << " , " << xv[0][i] << " , " << xv[1][i] << " , " << xv[0][i + 1] << " , " << xv[1][i + 1] << std::endl;
      }

      std::cout << "Delta=" << Delta << " " << epsilon << std::endl;

      if(fabs(Delta) > epsilon) {   // the edge does not pass for the origin
        std::cout << " xv[1][i]*xv[1][i+1] = " << xv[1][i]*xv[1][i + 1] << std::endl;
        if(fabs(xv[1][i]) < epsilon && xv[0][i] > 0) {  // the first vertex is on the positive x-axis
          std::cout << "the first vertex is on the positive x-axis" << std::endl;
          if(xv[1][i + 1] > 0) w += .5;
          else w -= .5;
        }
        else if(fabs(xv[1][i + 1]) < epsilon && xv[0][i + 1] > 0) { // the second vertex is on the positive x-axis
          std::cout << "the second vertex is on the positive x-axis" << std::endl;
          if(xv[1][i] < 0) w += .5;
          else w -= .5;
        }
        else if(xv[1][i]*xv[1][i + 1] < 0) { // the edge crosses the x-axis but doesn't pass through the origin
          double r = xv[0][i] - xv[1][i] * (xv[0][i + 1] - xv[0][i]) / (xv[1][i + 1] - xv[1][i]);
          std::cout << " r = " << r << std::endl;
          if(r > 0) {
            if(xv[1][i] < 0) w += 1.;
            else w -= 1;
          }
        }
      }
      else { // the line trought the edge passes for the origin
        std::cout << " xv[0][i]*xv[0][i+1] = " << xv[0][i]*xv[0][i + 1] << std::endl;
        if(fabs(xv[0][i]) < epsilon  && fabs(xv[1][i]) < epsilon) {  // vertex 1 is the origin
          w = 1; // set to 1 by default
          std::cout << "w set to 1 by default (vertex 1 is in the origin)" << std::endl;
	  break;
        }
        else if(fabs(xv[0][i + 1]) < epsilon && fabs(xv[1][i + 1]) < epsilon) { // vertex 2 is the origin
          w = 1; // set to 1 by default
          std::cout << "w set to 1 by default (vertex 2 is in the origin)" << std::endl;
	  break;
        }
        else if(xv[0][i] * xv[0][i + 1] < 0 || xv[1][i] * xv[1][i + 1] < 0) { //the edge crosses the origin
          w = 1; // set to 1 by default
          std::cout << "w set to 1 by default (the edge passes through the origin)" << std::endl;
	  break;
        }
      }
      std::cout << " w = " << w << " and iel = " << iel << std::endl;
    }
    return w;
  }

}
