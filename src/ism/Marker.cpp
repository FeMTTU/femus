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

    // let's look for the element of the mesh for which any of its vertices has least distance from the marker
    // let's initialize modulus and el (needed in this search)
    // modulus is the least distance from the origin of the vertices of iel and jel

    int iel = _mesh->_elementOffset[_iproc];
    int jel = iel + 1;

    std::vector< std::vector < double > > xiel(dim);    // will store the coordinates of the vertices of element iel
    std::vector< std::vector < double > > xjel(dim);    // will store the coordinates of the vertices of element jel
    std::vector< std::vector < double > > xkel(dim);    // will store the coordinates of the vertices of element kel
    for(unsigned k = 0; k < dim; k++) {
        xiel[k].reserve(9);
        xjel[k].reserve(9);
        xkel[k].reserve(9);
    }

    double modulusEl0 = 1e10; //modulus of element iel
    double modulusEl1 = 1e10; //modulus of element jel
    double modulus = 0; // least modulus among the two above
    int el = 0; // el will be the element from which we are going to start our search of the marker

    short unsigned ielType = _mesh->GetElementType(iel);
    short unsigned jelType = _mesh->GetElementType(jel);
    for(unsigned k = 0; k < dim; k++) {
        xiel[k].resize(facePointNumber[ielType]); //stores the coordinates of the vertices of element 0
        xjel[k].resize(facePointNumber[jelType]); //stores the coordinates of the vertices of element 1
    }

    // scale the coordinates of the vertices of element iel
    for(unsigned i = 0; i < facePointNumber[ielType]; i++) {
        unsigned iDof  = _mesh->GetSolutionDof(facePoints[ielType][i], iel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
            xiel[k][i] = (*_mesh->_topology->_Sol[k])(iDof) - _x[k];     // global extraction and local storage for the element coordinates
        }
    }

    //scale the coordinates of the vertices of element jel
    for(unsigned j = 0; j < facePointNumber[jelType]; j++) {
        unsigned jDof  = _mesh->GetSolutionDof(facePoints[jelType][j], jel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
            xjel[k][j] = (*_mesh->_topology->_Sol[k])(jDof) - _x[k];     // global extraction and local storage for the element coordinates
        }
    }

    for(unsigned i = 1; i < xiel[0].size() - 1; i++) {
        double modulusXv = sqrt(xiel[0][i]*xiel[0][i]+xiel[1][i]*xiel[1][i]);    // absolute value of vertex i of element iel
        std::cout << " modulusXv = " << modulusXv << " , " << " modulusEl0 = " << modulusEl0 << std::endl;
        modulusEl0 = (modulusXv < modulusEl0) ? modulusXv : modulusEl0;
    }

    for(unsigned j = 0; j < xjel[0].size() - 1; j++) {
        double modulusYv = sqrt(xjel[0][j]*xjel[0][j]+xjel[1][j]*xjel[1][j]);    // absolute value of vertex i of element jel = iel + 1
        std::cout << " modulusYv = " << modulusYv << " , " << " modulusEl1 = " << modulusEl1 << std::endl;
        modulusEl1 = (modulusYv < modulusEl1) ? modulusYv : modulusEl1;
    }
    std::cout << " modulusEl0 = " << modulusEl0 << " , " << " modulusEl1 = " << modulusEl1 << std::endl;

    el = (modulusEl0 <= modulusEl1) ? iel : jel;
    modulus = (modulusEl0 <= modulusEl1) ? modulusEl0 : modulusEl1;

    std::cout << " modulus = " << modulus << std::endl;

    // now let's find what element has the least distance from the origin (vertex wise)

    std::cout << " jel = " << jel << std::endl;
    std::cout << " el = " << el << std::endl;

    for(int kel = jel + 1; kel < _mesh->_elementOffset[_iproc + 1]; kel++) {
        double modulusKel = 1e10;
        short unsigned kelType = _mesh->GetElementType(kel);
        for(unsigned k = 0; k < dim; k++) {
            xkel[k].resize(facePointNumber[kelType]);
        }
        for(unsigned i = 0; i < facePointNumber[ielType]; i++) {
            unsigned kDof  = _mesh->GetSolutionDof(facePoints[kelType][i], kel, 2);    // global to global mapping between coordinates node and coordinate dof
            for(unsigned k = 0; k < dim; k++) {
                xkel[k][i] = (*_mesh->_topology->_Sol[k])(kDof) - _x[k];     // global extraction and local storage for the element coordinates
            }
        }
        for(unsigned i = 1; i < xkel[0].size() - 1; i++) {
            double modulusKel1 = sqrt(xkel[0][i]*xkel[0][i]+xkel[1][i]*xkel[1][i]);    // absolute value of vertex i of element kel
            modulusKel = (modulusKel1 < modulusKel) ? modulusKel1 : modulusKel;
        }
        std::cout << " modulus = " << modulus << " modulusKel = " << modulusKel << " el =  " << el << " kel =  " << kel <<  std::endl;
        el = (modulus <= modulusKel) ? el : kel;
        std::cout << " el = " << el << std::endl;
        modulus = (modulus <= modulusKel) ? modulus : modulusKel;
        std::cout << " modulus = " << modulus << " el =  " << el << " kel =  " << kel <<  std::endl;
        std::cout << " ----------------------------------------------------------------- " << std::endl;
    }

    std::cout << " el = " << el << std::endl;

    bool ElementHasBeenFound = false;

    if (modulus == 0) {
        std::cout << " The marker is one of the vertices of element " << el << std::endl;
        _elem = el;
        ElementHasBeenFound = true;
    }

    else {
        
        double w = GetWindingNumber(dim, el);
        if(w > 0) {
            _elem = el;
            ElementHasBeenFound = true;
        }
        else if(w < 0) {
            std::cout << " Error negative Winding Number with counterclockwise oriented points " << std::endl;
            abort();
        }
    }
    if(ElementHasBeenFound)
        std::cout << " The marker belongs to element " << _elem << std::endl;
    else {
        std::cout << " The marker does not belong to element " << el <<  std::endl;
    }
}

//     for(int iel = _mesh->_elementOffset[_iproc]; iel < _mesh->_elementOffset[_iproc + 1]; iel++) {
//       short unsigned ielType = _mesh->GetElementType(iel);
//       for(unsigned k = 0; k < dim; k++) {
//         xv[k].resize(facePointNumber[ielType]);
//       }
//       for(unsigned i = 0; i < facePointNumber[ielType]; i++) {
//         unsigned iDof  = _mesh->GetSolutionDof(facePoints[ielType][i], iel, 2);    // global to global mapping between coordinates node and coordinate dof
//         for(unsigned k = 0; k < dim; k++) {
//           xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof) - _x[k];     // global extraction and local storage for the element coordinates
//         }
//       }
//       double w = GetWindingNumber(xv, iel);
//       if(w > 0) {
//         _elem = iel;
//         ElementHasBeenFound = true;
//         break;
//       }
//       else if(w < 0) {
//         std::cout << "Error negative Winding Number with counterclockwise oriented points" << std::endl;
//         abort();
//       }
//
//     }
//     if(ElementHasBeenFound)
//       std::cout << "The marker belongs to element " << _elem << std::endl;
//     else {
//       std::cout << " The marker does not belong to this portion of the mesh" << std::endl;
//     }
//   }

double Marker::GetWindingNumber(const unsigned &dim, const int &iel) {
  
    std::vector< std::vector < double > > xv(dim);    // will store the coordinates of the vertices of element el
        
    for(unsigned k = 0; k < dim; k++) {
            xv[k].reserve(9);
        }
        short unsigned ielType = _mesh->GetElementType(iel);
        for(unsigned k = 0; k < dim; k++) {
            xv[k].resize(facePointNumber[ielType]);
        }
        for(unsigned i = 0; i < facePointNumber[ielType]; i++) {
            unsigned ielDof  = _mesh->GetSolutionDof(facePoints[ielType][i], iel, 2);    // global to global mapping between coordinates node and coordinate dof
            for(unsigned k = 0; k < dim; k++) {
                xv[k][i] = (*_mesh->_topology->_Sol[k])(ielDof) - _x[k];     // global extraction and local storage for the element coordinates
            }
        }  

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
