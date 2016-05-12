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

const unsigned facePoints3D[6][6][9] = {
  { {0, 8, 1, 17, 5, 12, 4, 16, 0},
    {1, 9, 2, 18, 6, 13, 5, 17, 1},
    {2, 10, 3, 19, 7, 14, 6, 18, 2},
    {3, 11, 0, 16, 4, 15, 7, 19, 3},
    {0, 11, 3, 10, 2, 9, 1, 8, 0},
    {4, 12, 5, 13, 6, 14, 7, 15, 4},
  },
  {{}}, {{}}, {{}}, {{}}, {{}}
};

const unsigned trianglesPerFace[6][6] = { //type - number of faces
  {8, 8, 8, 8, 8, 8},
  {6, 6, 6, 6},
  {8, 8, 8, 6, 6},
  {},
  {},
  {}
};

const unsigned faceTriangleNodes[6][6][8][4] = { // type - number of faces - number of triangles per faces - 3 vertices of the triangle (the first is counted twice)

  { {{0, 8, 20, 0}, {8, 1, 20, 8}, {1, 17, 20, 1}, {17, 5, 20, 17}, {5, 12, 20, 5}, {12, 4, 20, 12}, {4, 16, 20, 4 }, {16, 0, 20, 16}},
    {{1, 9, 21, 1}, {9, 2, 21, 9}, {2, 18, 21, 2}, {18, 6, 21, 18}, {6, 13, 21, 6}, {13, 5, 21, 13}, {5, 17, 21, 5}, {17, 1, 21, 17}},
    {{2, 10, 22, 2}, {10, 3, 22, 10}, {3, 19, 22, 3}, {19, 7, 22, 19}, {7, 14, 22, 7}, {14, 6, 22, 14}, {6, 18, 22, 6}, {18, 2, 22, 18}},
    {{3, 11, 23, 3}, {11, 0, 23, 11}, {0, 16, 23, 0}, {16, 4, 23, 16}, {4, 15, 23, 4}, {15, 7, 23, 15}, {7, 19, 23, 7}, {19, 3, 23, 19}},
    {{0, 24, 8, 0}, {8, 24, 1, 8}, {1, 24, 9, 1}, {9, 24, 2, 9}, {2, 24, 10, 2}, {10, 24, 3, 10}, {3, 24, 11, 3}, {11, 24, 0, 11}},
    {{4, 12, 25, 4}, {12, 5, 25, 12}, {5, 13, 25, 5}, {13, 6, 25, 13}, {6, 14, 25, 6}, {14, 7, 25, 14}, {7, 15, 25, 7}, {15, 4, 25, 15}}
  },
  { {{0, 6, 10, 0}, {6, 2, 10, 6}, {2, 5, 10, 2}, {5, 1, 10, 5}, {1, 4, 10, 1}, {4, 0, 10, 4}},
    {{0, 4, 11, 0}, {4, 1, 11, 4}, {1, 8, 11, 1}, {8, 3, 11, 8}, {3, 7, 11, 3}, {7, 0, 11, 7}},
    {{1, 5, 12, 1}, {5, 2, 12, 5}, {2, 9, 12, 2}, {9, 3, 12, 9}, {3, 8, 12, 3}, {8, 1, 12, 8}},
    {{2, 6, 13, 2}, {6, 0, 13, 6}, {0, 7, 13, 0}, {7, 3, 13, 7}, {3, 9, 13, 3}, {9, 2, 13, 9}}
  },
  { {{0, 6, 15, 0 }, {6, 1, 15, 6}, {1, 13, 15, 1}, {13, 4, 15, 13}, {4, 9, 15, 4}, {9, 3, 15, 9}, {3, 12, 15, 3}, {12, 0, 15, 12}},
    {{1, 7, 16, 1}, {7, 2, 16, 7}, {2, 14, 16, 2}, {14, 5, 16, 14}, {5, 10, 16, 5}, {10, 4, 16, 10}, {4, 13, 16, 4}, {13, 1, 16, 13}},
    {{2, 8, 17, 2}, {8, 0, 17, 8}, {0, 12, 17, 0}, {12, 3, 17, 12}, {3, 11, 17, 3}, {11, 5, 17, 11}, {5, 14, 17, 5}, {14, 2, 17, 14}},
    {{0, 8, 18, 0}, {8, 2, 18, 8}, {2, 7 , 18, 2}, {7, 1, 18, 7} , {1, 6 , 18, 1} , {6, 0, 18, 6}},
    {{3, 9, 19, 3}, {9, 4, 19, 9}, {4, 10, 19, 4}, {10, 5, 19, 10}, {5, 11, 19, 5}, {11, 3, 19, 11}}
  },
  {{{}}},
  {{{}}},
  {{{}}},
};

namespace femus {

  void Marker::GetElement(const bool &debug) {

    unsigned dim = _mesh->GetDimension();

    std::vector < unsigned > processorMarkerFlag(_nprocs, 3);

    //BEGIN SMART search
    // look to the closest element among a restricted list
    double modulus = 1.e10;
    unsigned iel = UINT_MAX;
    for(int jel = _mesh->_elementOffset[_iproc]; jel < _mesh->_elementOffset[_iproc + 1]; jel += 25) {

      unsigned interiorNode = _mesh->GetElementDofNumber(jel, 2) - 1;
      unsigned jDof  = _mesh->GetSolutionDof(interiorNode, jel, 2);    // global to global mapping between coordinates node and coordinate dof

      double distance2 = 0;
      for(unsigned k = 0; k < dim; k++) {
        double dk = (*_mesh->_topology->_Sol[k])(jDof) - _x[k];     // global extraction and local storage for the element coordinates
        distance2 += dk * dk;
      }

      double modulusKel = sqrt(distance2);

      if(modulusKel < modulus) {
        iel = jel;
        modulus = modulusKel;
      }
    }
    if(debug) {
      if(iel == UINT_MAX) {
        std::cout << "Warning the marker is located on unreasonable distance from the mesh >= 1.e10" << std::endl;
      }
      else {
        std::cout << "the smart search starts from element " << iel << std::endl;
      }
    }
    //END SMART search:


    bool elementHasBeenFound = false;
    bool pointIsOutsideThisProcess = false;
    bool pointIsOutsideTheDomain = false;

    std::vector< unsigned > nextElem(_nprocs, UINT_MAX);
    std::vector< unsigned > previousElem(_nprocs, UINT_MAX);
    previousElem[_iproc] = iel;
    unsigned nextProc = _iproc;

    while(!elementHasBeenFound) {

      //BEGIN next element search
      while(elementHasBeenFound + pointIsOutsideThisProcess + pointIsOutsideTheDomain == 0) {
        if(dim == 2) {
          nextElem[_iproc] = GetNextElement2D(dim, iel, previousElem[_iproc]);
        }
        else if(dim == 3) {
          nextElem[_iproc] = GetNextElement3D(dim, iel, previousElem[_iproc]);
        }
        previousElem[_iproc] = iel;

        if(nextElem[_iproc] == iel) {
          _elem = iel;
          elementHasBeenFound = true;
          processorMarkerFlag[_iproc] = 1;
        }
        else if(nextElem[_iproc] == UINT_MAX) {
          pointIsOutsideTheDomain = true;
          processorMarkerFlag[_iproc] = 0;
        }
        else {
          nextProc = _mesh->IsdomBisectionSearch(nextElem[_iproc], 3);
          if(nextProc != _iproc) {
            pointIsOutsideThisProcess = true;
            processorMarkerFlag[_iproc] = 2;
          }
          else {
            iel = nextElem[_iproc];
          }
        }
      }
      std::cout << std::flush;
      MPI_Barrier(PETSC_COMM_WORLD);

      if(debug) {
        if(elementHasBeenFound) {
          std::cout << " The marker belongs to element " << _elem << std::endl;
        }
        else if(pointIsOutsideTheDomain) {
          std::cout << " The marker does not belong to this domain" << std::endl;
        }
        else if(pointIsOutsideThisProcess) {
          std::cout << "proc " << _iproc << " believes the marker is in proc = " << nextProc << std::endl;
        }
      }
      //END next element search

      //BEGIN process exchange
      //send/receive if any process found the element
      for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
        if(jproc != _iproc) {
          if(processorMarkerFlag[_iproc] == 2 && jproc == nextProc) {
            unsigned three = 3;
            MPI_Send(&three, 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD);
          }
          else {
            MPI_Send(&processorMarkerFlag[_iproc], 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD);
          }
          MPI_Recv(&processorMarkerFlag[jproc], 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }
      // check if any process found the element
      unsigned sumFlag = 0;
      for(unsigned i = 0; i < _nprocs; i++) {
        sumFlag += processorMarkerFlag[i];
        if(processorMarkerFlag[i] == 1) {
          elementHasBeenFound = true;
          break;
        }
      }
      if(sumFlag == 0) {   // all the processes beleive that the marker is outside the domain
        std::cout << "Marker is outside the domain" << std::endl;
        _elem = UINT_MAX;
        if(debug) {
          for(unsigned j = 0 ; j < _nprocs; j++) {
            std::cout << " processorMarkerFlag[" << j << "] = " << processorMarkerFlag[j] <<  std::endl;
          }
        }
        break;
      }

      // _iproc sends its nextElem (which is in jproc) to jproc
      if(!elementHasBeenFound) {
        if(processorMarkerFlag[_iproc] == 2) {
          MPI_Send(&nextElem[_iproc], 1, MPI_UNSIGNED, nextProc, 1 , PETSC_COMM_WORLD);
          MPI_Send(&previousElem[_iproc], 1, MPI_UNSIGNED, nextProc, 2 , PETSC_COMM_WORLD);
        }
        for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
          if(processorMarkerFlag[jproc] == 3) {
            MPI_Recv(&nextElem[jproc], 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&previousElem[jproc], 1, MPI_UNSIGNED, jproc, 2 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          }
        }

        if(debug) {
          for(unsigned j = 0 ; j < _nprocs; j++) {
            std::cout << " processorMarkerFlag[" << j << "] = " << processorMarkerFlag[j] << "  " << nextElem[j] << std::endl;
          }
        }


        //BEGIN SMART search among all the elements received by _iproc and sent from jprocs
        double modulus = 1.e10;
        iel = _mesh->_elementOffset[_iproc + 1] ;
        for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
          if(processorMarkerFlag[jproc] == 3) {

            unsigned jel = nextElem[jproc];
            unsigned interiorNode = _mesh->GetElementDofNumber(jel, 2) - 1;

            unsigned jDof  = _mesh->GetSolutionDof(interiorNode, jel, 2);    // global to global mapping between coordinates node and coordinate dof
            double distance2 = 0;
            for(unsigned k = 0; k < dim; k++) {
              double dk = (*_mesh->_topology->_Sol[k])(jDof) - _x[k];     // global extraction and local storage for the element coordinates
              distance2 += dk * dk;
            }
            double modulusKel = sqrt(distance2);

            if(modulusKel < modulus) {
              iel = jel;
              previousElem[_iproc] = previousElem[jproc];
              modulus = modulusKel;

            }
          }
          //END SMART search
        }

        if(iel != _mesh->_elementOffset[_iproc + 1]) {
          std::cout << "start element= " << iel << std::endl;
          pointIsOutsideTheDomain = false;
          pointIsOutsideThisProcess = false;
          nextProc = _iproc;
          //previousElem = iel;
        }
        else {
          pointIsOutsideTheDomain = true;
          pointIsOutsideThisProcess = true;
          processorMarkerFlag[_iproc] = 0;
        }
      }

      //END process exchange

    }
  }


  unsigned Marker::GetNextElement2D(const unsigned &dim, const unsigned &currentElem, const unsigned &previousElem) {

    int nextElem;

    std::vector< std::vector < double > > xv(dim);

    for(unsigned k = 0; k < dim; k++) {
      xv[k].reserve(9);
    }
    short unsigned ielType = _mesh->GetElementType(currentElem);
    for(unsigned k = 0; k < dim; k++) {
      xv[k].resize(facePointNumber[ielType]);
    }
    for(unsigned i = 0; i < facePointNumber[ielType]; i++) {
      unsigned ielDof  = _mesh->GetSolutionDof(facePoints[ielType][i], currentElem, 2);    // global to global mapping between coordinates node and coordinate dof
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
      //std::cout << "Delta=" << Delta << " " << epsilon << std::endl;
      if(fabs(Delta) > epsilon) {   // the edge does not pass for the origin
        //std::cout << " xv[1][i]*xv[1][i+1] = " << xv[1][i]*xv[1][i + 1] << std::endl;
        if(fabs(xv[1][i]) < epsilon && xv[0][i] > 0) {  // the first vertex is on the positive x-axis
          //std::cout << "the first vertex is on the positive x-axis" << std::endl;
          if(xv[1][i + 1] > 0) w += .5;
          else w -= .5;
        }
        else if(fabs(xv[1][i + 1]) < epsilon && xv[0][i + 1] > 0) { // the second vertex is on the positive x-axis
          //std::cout << "the second vertex is on the positive x-axis" << std::endl;
          if(xv[1][i] < 0) w += .5;
          else w -= .5;
        }
        else if(xv[1][i]*xv[1][i + 1] < 0) { // the edge crosses the x-axis but doesn't pass through the origin
          double r = xv[0][i] - xv[1][i] * (xv[0][i + 1] - xv[0][i]) / (xv[1][i + 1] - xv[1][i]);
          //std::cout << " r = " << r << std::endl;
          if(r > 0) {
            if(xv[1][i] < 0) w += 1.;
            else w -= 1;
          }
        }
      }
      else { // the line trought the edge passes for the origin
        //std::cout << " xv[0][i]*xv[0][i+1] = " << xv[0][i]*xv[0][i + 1] << std::endl;
        if(fabs(xv[0][i]) < epsilon  && fabs(xv[1][i]) < epsilon) {  // vertex 1 is the origin
          w = 1; // set to 1 by default
          //std::cout << "w set to 1 by default (vertex 1 is in the origin)" << std::endl;
          break;
        }
        else if(fabs(xv[0][i + 1]) < epsilon && fabs(xv[1][i + 1]) < epsilon) { // vertex 2 is the origin
          w = 1; // set to 1 by default
          //std::cout << "w set to 1 by default (vertex 2 is in the origin)" << std::endl;
          break;
        }
        else if(xv[0][i] * xv[0][i + 1] < 0 || xv[1][i] * xv[1][i + 1] < 0) { //the edge crosses the origin
          w = 1; // set to 1 by default
          //std::cout << "w set to 1 by default (the edge passes through the origin)" << std::endl;
          break;
        }
      }
      //std::cout << " w = " << w << " and currentElem = " << currentElem << std::endl;
    }

    nextElem = currentElem;

    if(w == 0) {

      double distance = 1.e10;
      for(unsigned j = 1; j < xv[0].size() - 1; j += 2) {
        double distancej = 0.;
        for(unsigned k = 0; k < dim; k++) {
          distancej += xv[k][j] * xv[k][j];
        }
        distancej = sqrt(distancej);

        if(distancej < distance) {
          int jel = (_mesh->el->GetFaceElementIndex(currentElem, (j - 1) / 2) - 1);
          if(jel != previousElem) {
            nextElem = jel;
            distance = distancej;
          }
        }
      }

    }
    else if(w < 0) {
      //std::cout << " Error negative Winding Number with counterclockwise oriented points " << std::endl;
      abort();
    }

    std::cout << "the next element is " << nextElem << std::endl;


    return (nextElem >= 0) ? nextElem : UINT_MAX;
  }

  unsigned Marker::GetNextElement3D(const unsigned &dim, const unsigned &currentElem, const unsigned &previousElem) {

    int nextElem;
    //for on the element faces
    unsigned faceIntersectionCounter = 0;
    unsigned faceIntersectionCounterOld = 0; // if it is even, the marker is outside the element (0 is considered even)
    bool markerIsInElement = false;


    for(unsigned iface = 0; iface < _mesh->GetElementFaceNumber(currentElem); iface++) {

      // std::cout << " iface = " << iface << std::endl;

      short unsigned currentElementType = _mesh->GetElementType(currentElem);

      for(unsigned itri = 0; itri < trianglesPerFace[currentElementType][iface]; itri ++) {

        // std::cout << "faceIntersectionCounter  = " << faceIntersectionCounter  << " , " << " markerIsInElement = " << markerIsInElement <<  std::endl;

        faceIntersectionCounterOld = faceIntersectionCounter ;
        bool lineIntersection = false;

        std::vector< std::vector < double > > xv(dim); //stores the coordinates of the nodes of the triangle itri

        for(unsigned k = 0; k < dim; k++) {
          //xv[k].reserve(9);
          xv[k].reserve(4);
        }
//       short unsigned currentElementType = _mesh->GetElementType(currentElem);
//       short unsigned ifaceType = _mesh->GetElementFaceType(currentElem, iface);
        for(unsigned k = 0; k < dim; k++) {
          // xv[k].resize(facePointNumber[ifaceType]);
          xv[k].resize(4);
        }
//             for(unsigned i = 0; i < facePointNumber[ifaceType]; i++) {
//                 unsigned ifaceDof  = _mesh->GetSolutionDof(facePoints3D[currentElementType][iface][i], currentElem, 2);
//                 // std::cout << "ifaceDof = " << ifaceDof << std::endl;
//                 for(unsigned k = 0; k < dim; k++) {
//                     xv[k][i] = (*_mesh->_topology->_Sol[k])(ifaceDof) - _x[k];     // global extraction and local storage for the element coordinates
//                 }
//             }

        for(unsigned i = 0; i < 4; i++) {
          unsigned itriDof  = _mesh->GetSolutionDof(faceTriangleNodes[currentElementType][iface][itri][i], currentElem, 2);
          //std::cout << "itriDof = " << itriDof << std::endl;
          for(unsigned k = 0; k < dim; k++) {
            xv[k][i] = (*_mesh->_topology->_Sol[k])(itriDof) - _x[k];     // global extraction and local storage for the element coordinates
          }
        }


        double length = 0.;
        //  for(unsigned i = 0; i < xv[0].size() - 1; i++) {
        for(unsigned i = 0; i < 3; i++) {
          length += sqrt((xv[0][i + 1] - xv[0][i]) * (xv[0][i + 1] - xv[0][i]) +
                         (xv[1][i + 1] - xv[1][i]) * (xv[1][i + 1] - xv[1][i]) +
                         (xv[2][i + 1] - xv[2][i]) * (xv[2][i + 1] - xv[2][i]));
        }

//             length /= xv[0].size();
        length /= 4;
        // for(unsigned i = 0; i < xv[0].size(); i++) {
        for(unsigned i = 0; i < 4; i++) {
          xv[0][i] /= length;
          xv[1][i] /= length;
          xv[2][i] /= length;
        }

//             // now we have to check if the plane of the face intersects the positive x-axis
//             // let's find the plane passing through xv[][0], xv[][2] and xv[][4] (they will not be aligned)
//
//             // entries of the normal to the plane
//             double A = -(xv[1][4] - xv[1][0]) * (xv[2][2] - xv[2][0]) + (xv[2][4] - xv[2][0]) * (xv[1][2] - xv[1][0]);
//             double B = -(xv[2][4] - xv[2][0]) * (xv[0][2] - xv[0][0]) + (xv[0][4] - xv[0][0]) * (xv[2][2] - xv[2][0]);
//             double C = -(xv[0][4] - xv[0][0]) * (xv[1][2] - xv[1][0]) + (xv[1][4] - xv[1][0]) * (xv[0][2] - xv[0][0]);
//
//             double epsilon  = 10e-10;
//             double epsilon2  = epsilon * epsilon;
//             double By0 = B * xv[1][0];
//             double Cz0 = C * xv[2][0];

        // now we have to check if the plane of the triangle itri intersects the positive x-axis
        // let's find the plane passing through the vertices of the triangle itri

        // entries of the normal to the plane
        double A = -(xv[1][2] - xv[1][0]) * (xv[2][1] - xv[2][0]) + (xv[2][2] - xv[2][0]) * (xv[1][1] - xv[1][0]);
        double B = -(xv[2][2] - xv[2][0]) * (xv[0][1] - xv[0][0]) + (xv[0][2] - xv[0][0]) * (xv[2][1] - xv[2][0]);
        double C = -(xv[0][2] - xv[0][0]) * (xv[1][1] - xv[1][0]) + (xv[1][2] - xv[1][0]) * (xv[0][1] - xv[0][0]);

        double epsilon  = 10e-10;
        double epsilon2  = epsilon * epsilon;
        double By0 = B * xv[1][0];
        double Cz0 = C * xv[2][0];


        // std::cout << "A = " << A << " , " << "By0 = " << By0 << " , " << " Cz0 = " << Cz0 <<  std::endl;


        if(fabs(A) < epsilon && epsilon < fabs(By0 + Cz0)) { // A = 0 and By0 != -Cz0
          // std::cout << "The plane of face " << iface << "and the x-axis don't intersect " << std::endl;

        }
        else {

          std::vector < double > xTilde(dim, 0);
          double r = 0;

          if(fabs(A) < epsilon && fabs(By0 + Cz0) < epsilon) { // A = 0 and By0 = -Cz0
            //   std::cout << "The plane of face " << iface << "and the x-axis intersect on a line" << std::endl;
            // the marker and the face are already on the same plane so there is no need for further shifting
            lineIntersection = true ;
          }

          else if(epsilon < fabs(A)) {  // A != 0
            // std::cout << "The plane of face " << iface << "and the x-axis intersect at a point" << std::endl;

            r = (A * xv[0][0] + B * xv[1][0] + C * xv[2][0]) / A;
            xTilde[0] = r;
            if(fabs(r) < epsilon) { // r = 0
              lineIntersection = true; // the intersection point is the actual marker
            }
          }
          if(r > 0 || fabs(r) < epsilon) {

            //  for(unsigned i = 0; i < facePointNumber[ifaceType]; i++) {
            for(unsigned i = 0; i < 4; i++) {
              for(unsigned k = 0; k < dim; k++) {
                xv[k][i] = xv[k][i] - xTilde[k];     //transate again the reference frame so that the origin is xTilde
              }
            }

            unsigned scalarCount = 0;

            // for(unsigned i = 0; i < xv[0].size() - 1; i++) {
            for(unsigned i = 0; i < 3; i++) {
              //entries of the vector (xTilde - xi) X ( xi+1 -xi)
              double q0 = xv[1][i] * (xv[2][i] - xv[2][i + 1]) + xv[2][i] * (xv[1][i + 1] - xv[1][i]);
              double q1 = xv[2][i] * (xv[0][i] - xv[0][i + 1]) + xv[0][i] * (xv[2][i + 1] - xv[2][i]);
              double q2 = xv[0][i] * (xv[1][i] - xv[1][i + 1]) + xv[1][i] * (xv[0][i + 1] - xv[0][i]);

              //  std::cout << "q0 = " << q0 << " , " << "q1 = " << q1 << " , " << " q2 = " << q2 <<  std::endl;

              double  scalarProduct = q0 * A + q1 * B + q2 * C;

              //  std::cout << "scalarProduct = " << scalarProduct << std::endl;

              if(scalarProduct > 0) {  // this can be cancelled once the code is working ok
                //    std::cout << "the marker is outside face " << iface <<  std::endl;
                break;

              }
              else if(fabs(scalarProduct) < epsilon) {  //scalarProduct == 0
                std::cout << " the marker and the edge are aligned " << std::endl; //check if xTilde is actually on the edge.
                if(xv[0][i]*xv[0][i + 1] < 0 || xv[1][i]*xv[1][i + 1] < 0 || xv[2][i]*xv[2][i + 1] < 0) {
                  if(lineIntersection == true) {
                    // std::cout << " the marker belongs to an edge of face " << iface << std::endl;
                    markerIsInElement = true;
                    break;
                  }
                  else {
                    faceIntersectionCounter++;
                    // std::cout << " faceIntersectionCounter = " << faceIntersectionCounter << std::endl;
                    break;
                  }

                }
                else if((xv[0][i] * xv[0][i]  + xv[1][i] * xv[1][i] + xv[2][i] * xv[2][i]) < epsilon2 ||
                        (xv[0][i + 1]*xv[0][i + 1] + xv[1][i + 1]*xv[1][i + 1] + xv[2][i + 1]*xv[2][i + 1]) < epsilon2) {
                  if(lineIntersection == true) {
                    // std::cout << " one of the vertices is the marker" << std::endl;
                    markerIsInElement = true;
                    break;
                  }
                  else {
                    faceIntersectionCounter++;
                    //std::cout << " faceIntersectionCounter " << faceIntersectionCounter << std::endl;
                    break;
                  }
                }
              }
              else if(scalarProduct < 0) {
                //     std::cout << "increase scalarCount" <<std::endl;
                scalarCount++;
                //  std::cout << "scalarCount = " << scalarCount << std::endl;
              }
            }

            if(faceIntersectionCounterOld < faceIntersectionCounter) {
              break;
            }
//           if(scalarCount == acePointNumber[ifaceType] - 1 && lineIntersection == true) {
            if(scalarCount == 3 && lineIntersection == true) {
              markerIsInElement = true ;
              break;
            }
            //else if(scalarCount == facePointNumber[ifaceType] - 1 && lineIntersection == false) {
            else if(scalarCount == 3 && lineIntersection == false) {
              faceIntersectionCounter++;
              break;
            }
          }

        } //end of else
        if(markerIsInElement == true) {
          break;
        }

      } //end of the loop on tri

      if(markerIsInElement == true) {
        break;
      }
    }// end of the loop on iface

    //std::cout << "markerIsInElement = " << markerIsInElement << " and faceIntersectionCounter  = " << faceIntersectionCounter  <<  std::endl;

    if(markerIsInElement == true || faceIntersectionCounter % 2 != 0) {
      nextElem = currentElem;
    }
    else if(markerIsInElement == false && faceIntersectionCounter % 2 == 0) {
      // std::cout << " The marker doesn't belong to element " << currentElem << std::endl;
      double modulus = 1.e10;

      for(unsigned iface = 0; iface < _mesh->GetElementFaceNumber(currentElem); iface++) {

        // double xg[3] = {0., 0., 0.};
        for(int i = 0; i < 3; i ++) {
          unsigned j = _mesh->GetLocalFaceVertexIndex(currentElem, iface, i);
          unsigned iDof = _mesh->GetSolutionDof(j, currentElem, 2);

          //  for(unsigned k = 0; k < dim; k++) {
          //std::cout << (*_mesh->_topology->_Sol[k])(iDof)  << " ";
          // xg[k] += 1. / 3.*(*_mesh->_topology->_Sol[k])(iDof); ///COS 'E' QUESTO ??
          //  }
          //std::cout << std::endl;
        }
        // std::cout << "xg= " << xg[0] << " " << xg[1] << " " << xg[2] << std::endl;

        unsigned faceNodeNumber = _mesh->GetElementFaceDofNumber(currentElem, iface, 2);
        unsigned i = _mesh->GetLocalFaceVertexIndex(currentElem, iface, faceNodeNumber - 1);
        //std::cout << faceNodeNumber-1 <<" "<<i << " ";

        unsigned faceCentralDof = _mesh->GetSolutionDof(i, currentElem, 2);

        //std::cout << "faceCentralDof =" << faceCentralDof << std::endl;
        double distance2 = 0;
        for(unsigned k = 0; k < dim; k++) {
          // std::cout << (*_mesh->_topology->_Sol[k])(faceCentralDof)  << " ";
          double dk = (*_mesh->_topology->_Sol[k])(faceCentralDof) - _x[k];     // global extraction and local storage for the element coordinates
          distance2 += dk * dk;
        }
        //std::cout << std::endl;
        double ifaceModulus = sqrt(distance2);

        // std::cout << faceCentralDof << " " << ifaceModulus << std::endl;

        if(ifaceModulus < modulus) {
          int jel = (_mesh->el->GetFaceElementIndex(currentElem, iface) - 1);
          // std::cout << "jel = " << jel << "iface = " << iface <<  std::endl;
          if(jel != previousElem) {
            nextElem = jel;
            modulus = ifaceModulus;
          }
        }
      }
    }
    std::cout << "nextElem = " << nextElem << std::endl;

    return (nextElem >= 0) ? nextElem : UINT_MAX;
  }



  std::vector< double > Marker::InverseMappingQuad(const unsigned &currentElem, const unsigned &solutionType,
      std::vector< double > &x) {

    unsigned dim = 2;
    std::vector< std::vector < double > > xv(dim);
    std::vector < std::vector < double > > a(dim);

    unsigned nDofs = _mesh->GetElementDofNumber(currentElem, solutionType);
    short unsigned currentElementType = _mesh->GetElementType(currentElem);

    for(unsigned k = 0; k < dim; k++) {
      xv[k].resize(nDofs);
      a[k].resize(nDofs);
    }
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned iDof  = _mesh->GetSolutionDof(i, currentElem, 2);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
      }
    }


    if(solutionType == 0) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = 0.25 * (xv[k][0] + xv[k][1] + xv[k][2] + xv[k][3]);
        a[k][1] = 0.25 * (-xv[k][0] + xv[k][1] + xv[k][2] - xv[k][3]);
        a[k][2] = 0.25 * (-xv[k][0] - xv[k][1] + xv[k][2] + xv[k][3]);
        a[k][3] = 0.25 * (xv[k][0] - xv[k][1] + xv[k][2] - xv[k][3]);
      }
    }
    else if(solutionType == 1) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = -0.25 * (xv[k][0] + xv[k][1] + xv[k][2] + xv[k][3]) +
                  0.5 * (xv[k][4] + xv[k][5] + xv[k][6] + xv[k][7]);
        a[k][1] = 0.5 * (xv[k][5] - xv[k][7]);
        a[k][2] = 0.5 * (xv[k][6] - xv[k][4]);
        a[k][3] = 0.25 * (xv[k][0] - xv[k][1] + xv[k][2] - xv[k][3]);
        a[k][4] = 0.25 * (xv[k][0] + xv[k][1] + xv[k][2] + xv[k][3]) - 0.5 * (xv[k][4] + xv[k][6]);
        a[k][5] = 0.25 * (xv[k][0] + xv[k][1] + xv[k][2] + xv[k][3]) - 0.5 * (xv[k][5] + xv[k][7]);
        a[k][6] = 0.25 * (-xv[k][0] - xv[k][1] + xv[k][2] + xv[k][3]) + 0.5 * (xv[k][4] - xv[k][6]);
        a[k][7] = 0.25 * (-xv[k][0] + xv[k][1] + xv[k][2] - xv[k][3]) + 0.5 * (xv[k][7] - xv[k][5]);
      }
    }
    else if(solutionType == 2) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = xv[k][8];
        a[k][1] = 0.5 * (xv[k][5] - xv[k][7]);
        a[k][2] = 0.5 * (xv[k][6] - xv[k][4]);
        a[k][3] = 0.25 * (xv[k][0] - xv[k][1] + xv[k][2] - xv[k][3]);
        a[k][4] =  0.5 * (xv[k][5] + xv[k][7]) - xv[k][8];
        a[k][5] =  0.5 * (xv[k][4] + xv[k][6]) - xv[k][8];
        a[k][6] = 0.25 * (-xv[k][0] - xv[k][1] + xv[k][2] + xv[k][3]) +
                  0.5 * (xv[k][4] - xv[k][6]);
        a[k][7] = 0.25 * (-xv[k][0] + xv[k][1] + xv[k][2] - xv[k][3]) +
                  0.5 * (-xv[k][5] + xv[k][7]);
        a[k][8] = 0.25 * (xv[k][0] + xv[k][1] + xv[k][2] + xv[k][3]) -
                  0.5 * (xv[k][4] + xv[k][5] + xv[k][6] + xv[k][7]) + xv[k][8];
      }
    }

    std::vector < double > phi(nDofs);

    double xi = x[0];
    double eta = x[1];

    phi[0] = 1.;
    phi[1] = xi; // x
    phi[2] = eta;  // y
    phi[3] = xi * eta;  // x y
    if(solutionType > 0) {
      phi[4] = xi * xi; // x x
      phi[5] = eta * eta; // y y
      phi[6] = phi[4] * eta; // xx y
      phi[7] = phi[5] * xi; // x yy
    }
    if(solutionType > 1) {
      phi[8] = phi[3] * phi[3]; // xx yy
    }

    std::vector < double > xp(dim);
    for(int i = 0; i < nDofs; i++) {
      for(int k = 0; k < dim; k++) {
        xp[k] += a[k][i] * phi[i];
      }
    }
    return xp;

  }


  std::vector< double > Marker::InverseMappingTri(const unsigned &currentElem, const unsigned &solutionType,
      std::vector< double > &x) {

    unsigned dim = 2;
    std::vector< std::vector < double > > xv(dim);
    std::vector < std::vector < double > > a(dim);

    unsigned nDofs = _mesh->GetElementDofNumber(currentElem, solutionType);
    short unsigned currentElementType = _mesh->GetElementType(currentElem);

    for(unsigned k = 0; k < dim; k++) {
      xv[k].resize(nDofs);
      a[k].resize(nDofs);
    }
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned iDof  = _mesh->GetSolutionDof(i, currentElem, 2);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
      }
    }


    if(solutionType == 0) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = xv[k][0];
        a[k][1] = - xv[k][0] + xv[k][1];
        a[k][2] = - xv[k][0] + xv[k][2];
      }
    }
    else if(solutionType == 1) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = xv[k][0];
        a[k][1] = - 3 * xv[k][0] - xv[k][1] + 4 * xv[k][3];
        a[k][2] = - 3 * xv[k][0] - xv[k][2] + 4 * xv[k][5];
        a[k][3] = 4 * xv[k][0] - 4 * xv[k][3] + 4 * xv[k][4] - 4 * xv[k][5];
        a[k][4] = 2 * xv[k][0] + 2 * xv[k][1] - 4 * xv[k][3];
        a[k][5] = 2 * xv[k][0] + 2 * xv[k][2] - 4 * xv[k][5];
      }
    }
    else if(solutionType == 2) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = xv[k][0];
        a[k][1] = - 3 * xv[k][0] - xv[k][1] + 4 * xv[k][3];
        a[k][2] = - 3 * xv[k][0] - xv[k][2] + 4 * xv[k][5];
        a[k][3] = 7 * xv[k][0] + 3 * xv[k][1] + 3 * xv[k][2] - 16 * xv[k][3] -
                  8 * xv[k][4] - 16 * xv[k][5] + 27 * xv[k][6];
        a[k][4] = 2 * xv[k][0] + 2 * xv[k][1] - 4 * xv[k][3];
        a[k][5] = 2 * xv[k][0] + 2 * xv[k][2] - 4 * xv[k][5];
        a[k][6] = - 3 * xv[k][0] - 3 * xv[k][1] - 3 * xv[k][2] + 12 * xv[k][3] +
                  12 * xv[k][4] + 12 * xv[k][5] - 27 * xv[k][6];
      }
    }

    std::vector < double > phi(nDofs);

    double xi = x[0];
    double eta = x[1];

    phi[0] = 1.;
    phi[1] = xi; // x
    phi[2] = eta; // y
    if(solutionType > 0) {
      phi[3] = xi * eta; // x y
      phi[4] = xi * xi;  // x x
      phi[5] = eta * eta; // y y
    }
    if(solutionType > 1) {
      phi[6] = phi[4] * eta + phi[5] * xi; // xx y + x yy
    }

    std::vector < double > xp(dim);
    for(int i = 0; i < nDofs; i++) {
      for(int k = 0; k < dim; k++) {
        xp[k] += a[k][i] * phi[i];
      }
    }
    return xp;
  }


  std::vector< double > Marker::InverseMappingHex(const unsigned &currentElem, const unsigned &solutionType,
      std::vector< double > &x) {

    unsigned dim = 3;
    std::vector< std::vector < double > > xv(dim);
    std::vector < std::vector < double > > a(dim);

    unsigned nDofs = _mesh->GetElementDofNumber(currentElem, solutionType);
    short unsigned currentElementType = _mesh->GetElementType(currentElem);

    for(unsigned k = 0; k < dim; k++) {
      xv[k].resize(nDofs);
      a[k].resize(nDofs);
    }
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned iDof  = _mesh->GetSolutionDof(i, currentElem, 2);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
      }
    }


    if(solutionType == 0) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = 0.125 * (xv[k][0] + xv[k][1] + xv[k][2] + xv[k][3] + xv[k][4] +
                           xv[k][5] + xv[k][6] + xv[k][7]) ;
        a[k][1] = 0.125 * (- xv[k][0] + xv[k][1] + xv[k][2] - xv[k][3] - xv[k][4] +
                           xv[k][5] + xv[k][6] - xv[k][7]) ;
        a[k][2] = 0.125 * (- xv[k][0] - xv[k][1] + xv[k][2] + xv[k][3] - xv[k][4] -
                           xv[k][5] + xv[k][6] + xv[k][7]) ;
        a[k][3] = 0.125 * (- xv[k][0] - xv[k][1] - xv[k][2] - xv[k][3] + xv[k][4] +
                           xv[k][5] + xv[k][6] + xv[k][7]) ;
        a[k][4] = 0.125 * (xv[k][0] - xv[k][1] + xv[k][2] - xv[k][3] + xv[k][4] -
                           xv[k][5] + xv[k][6] - xv[k][7]) ;
        a[k][5] = 0.125 * (xv[k][0] - xv[k][1] - xv[k][2] + xv[k][3] - xv[k][4] +
                           xv[k][5] + xv[k][6] - xv[k][7]) ;
        a[k][6] = 0.125 * (xv[k][0] + xv[k][1] - xv[k][2] - xv[k][3] - xv[k][4] -
                           xv[k][5] + xv[k][6] + xv[k][7]) ;
        a[k][7] = 0.125 * (- xv[k][0] + xv[k][1] - xv[k][2] + xv[k][3] + xv[k][4] -
                           xv[k][5] + xv[k][6] - xv[k][7]) ;

      }
    }
    else if(solutionType == 1) { //questo non funziona quando si fanno i test

      for(int k = 0; k < dim; k++) {

        a[k][0] = 0.25 * (- xv[k][0] + xv[k][9] + xv[k][10] + xv[k][11] + xv[k][12] + xv[k][13] + xv[k][14] +
                          xv[k][15] + xv[k][16] + xv[k][17] + xv[k][18] - xv[k][1] + xv[k][19] - xv[k][2] -
                          xv[k][3] - xv[k][4] - xv[k][5] - xv[k][6] - xv[k][7] + xv[k][8]);
        a[k][1] = 0.125 * (xv[k][0] + 2 * xv[k][9] - 2 * xv[k][11] + 2 * xv[k][13] - 2 * xv[k][15] - 2 * xv[k][16] +
                           2 * xv[k][17] + 2 * xv[k][18] - xv[k][1] - 2 * xv[k][19] - xv[k][2] + xv[k][3] +
                           xv[k][4] - xv[k][5] - xv[k][6] + xv[k][7]);
        a[k][2] = 0.125 * (xv[k][0] + 2 * xv[k][10] - 2 * xv[k][12] + 2 * xv[k][14] - 2 * xv[k][16] - 2 * xv[k][17] +
                           2 * xv[k][18] + xv[k][1] + 2 * xv[k][19] - xv[k][2] - xv[k][3] + xv[k][4] + xv[k][5] - xv[k][6] - xv[k][7] - 2 * xv[k][8]);
        a[k][3] = 0.125 * (xv[k][0] - 2 * xv[k][9] - 2 * xv[k][10] - 2 * xv[k][11] + 2 * xv[k][12] + 2 * xv[k][13] + 2 * xv[k][14] +
                           2 * xv[k][15] + xv[k][1] + xv[k][2] + xv[k][3] - xv[k][4] - xv[k][5] - xv[k][6] - xv[k][7] - 2 * xv[k][10]);
        a[k][4] = 0.25 * (xv[k][16] - xv[k][17] + xv[k][18] - xv[k][19]);
        a[k][5] = 0.25 * (- xv[k][9] + xv[k][11] + xv[k][13] - xv[k][15]);
        a[k][6] = 0.25 * (- xv[k][10] - xv[k][12] + xv[k][14] + xv[k][8]);
        a[k][7] = 0.125 * (xv[k][0] - 2 * xv[k][10] - 2 * xv[k][12] - 2 * xv[k][14] + xv[k][1] + xv[k][2] + xv[k][3] + xv[k][4] +
                           xv[k][5] + xv[k][6] + xv[k][7] - 2 * xv[k][8]);
        a[k][8] = 0.125 * (xv[k][0] - 2 * xv[k][9] - 2 * xv[k][11] - 2 * xv[k][13] - 2 * xv[k][15] + xv[k][1] + xv[k][2] + xv[k][3] +
                           xv[k][4] + xv[k][5] + xv[k][6] + xv[k][7]);
        a[k][9] = 0.125 * (xv[k][0] - 2  * xv[k][16] - 2 * xv[k][17] - 2 *  xv[k][18] + xv[k][1] - 2 * xv[k][19] + xv[k][2] + xv[k][3] +
                           xv[k][4] + xv[k][5] + xv[k][6] + xv[k][7]);
        a[k][10] = 0.125 * (- xv[k][0] + xv[k][1] - xv[k][2] + xv[k][3] + xv[k][4] - xv[k][5] + xv[k][6] - xv[k][7]);
        a[k][11] = 0.125 * (- xv[k][0] - 2 * xv[k][10] + 2 *  xv[k][12] - 2 * xv[k][14] - xv[k][1] + xv[k][2] + xv[k][3] - xv[k][4] -
                            xv[k][5] + xv[k][6] + xv[k][7] + 2 * xv[k][8]);
        a[k][12] = 0.125 * (- xv[k][0] + 2 * xv[k][10] - 2 * xv[k][12] - 2 * xv[k][14] - xv[k][1] - xv[k][2] - xv[k][3] + xv[k][4] +
                            xv[k][5] + xv[k][6] + xv[k][7] + 2 * xv[k][8]);
        a[k][13] = 0.125 * (- xv[k][0] + 2 * xv[k][9] + 2 * xv[k][11] - 2 * xv[k][13] - 2 * xv[k][15] - xv[k][1] - xv[k][2] - xv[k][3] +
                            xv[k][4] + xv[k][5] + xv[k][6] + xv[k][7]);
        a[k][14] = 0.125 * (- xv[k][0] - 2 * xv[k][9] + 2 * xv[k][11] - 2 * xv[k][13] + 2 * xv[k][15] + xv[k][1] + xv[k][2] - xv[k][3] -
                            xv[k][4] + xv[k][5] + xv[k][6] - xv[k][7]);
        a[k][15] = 0.125 * (- xv[k][0] + 2 * xv[k][16] - 2 * xv[k][17] - 2 * xv[k][18] + xv[k][1] + 2 * xv[k][19] + xv[k][2] - xv[k][3] - xv[k][4] +
                            xv[k][5] + xv[k][6] - xv[k][7]);
        a[k][16] = 0.125 * (- xv[k][0] + 2 * xv[k][16] + 2 * xv[k][17] - 2 * xv[k][18] - xv[k][1] - 2 * xv[k][19] + xv[k][2] + xv[k][3] - xv[k][4] -
                            xv[k][5] + xv[k][6] + xv[k][7]);
        a[k][17] = 0.125 * (xv[k][0] + 2 * xv[k][10] + 2 * xv[k][12] - 2 * xv[k][14] + xv[k][1] - xv[k][2] - xv[k][3] - xv[k][4] - xv[k][5] + xv[k][6] +
                            xv[k][7] - 2 * xv[k][8]);
        a[k][18] = 0.125 * (xv[k][0] + 2 * xv[k][9] - 2 * xv[k][11] - 2 * xv[k][13] + 2 * xv[k][15] - xv[k][1] - xv[k][2] + xv[k][3] - xv[k][4] +
                            xv[k][5] + xv[k][6] - xv[k][7]);
        a[k][19] = 0.125 * (xv[k][0] - 2 * xv[k][16] + 2 * xv[k][17] - 2 * xv[k][18] - xv[k][1] + 2 * xv[k][19] + xv[k][2] - xv[k][3] + xv[k][4] -
                            xv[k][5] + xv[k][6] - xv[k][7]);
      }
    }
    else if(solutionType == 2) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = xv[k][26];
        a[k][1] = 0.5 * (xv[k][21] - xv[k][23]);
        a[k][2] = 0.5 * (xv[k][22] - xv[k][20]);
        a[k][3] = 0.5 * (xv[k][25] - xv[k][24]);
        a[k][4] = 0.25 * (xv[k][16] - xv[k][17] + xv[k][18] - xv[k][19]);
        a[k][5] = 0.25 * (- xv[k][9] + xv[k][11] + xv[k][13] - xv[k][15]);
        a[k][6] = 0.25 * (- xv[k][10] - xv[k][12] + xv[k][14] + xv[k][8]);
        a[k][7] = 0.5 * (xv[k][21] + xv[k][23] - 2 * xv[k][26]);
        a[k][8] = 0.5 * (xv[k][20] + xv[k][22] - 2 * xv[k][26]);
        a[k][9] = 0.5 * (xv[k][24] + xv[k][25] - 2 * xv[k][26]);
        a[k][10] = 0.125 * (- xv[k][0] + xv[k][1] - xv[k][2] + xv[k][3] + xv[k][4] - xv[k][5] + xv[k][6] - xv[k][7]);
        a[k][11] = 0.25 * (- xv[k][16] - xv[k][17] + xv[k][18] + xv[k][19] + 2 * xv[k][20] - 2 * xv[k][22]);
        a[k][12] = 0.25 * (- xv[k][9] - xv[k][11] + xv[k][13] + xv[k][15] + 2 * xv[k][24] - 2 * xv[k][25]);
        a[k][13] = 0.25 * (- xv[k][10] + xv[k][12] + xv[k][14] + 2 * xv[k][24] - 2 * xv[k][25] - xv[k][8]);
        a[k][14] = 0.25 * (- xv[k][16] + xv[k][17] + xv[k][18] - xv[k][19] - 2 * xv[k][21] + 2 * xv[k][23]);
        a[k][15] = 0.25 * (xv[k][9] - xv[k][11] + xv[k][13] - xv[k][15] - 2 * xv[k][21] + 2 * xv[k][23]);
        a[k][16] = 0.25 * (xv[k][10] - xv[k][12] + xv[k][14] + 2 * xv[k][20] - 2 * xv[k][22] - xv[k][8]);
        a[k][17] = 0.125 * (xv[k][0] + 2 * xv[k][10] + 2 * xv[k][12] - 2 * xv[k][14] + xv[k][1] - xv[k][2] - xv[k][3] - xv[k][4] - xv[k][5] + xv[k][6] + xv[k][7] - 2 * xv[k][8]);
        a[k][18] = 0.125 * (xv[k][0] + 2 * xv[k][9] - 2 * xv[k][11] - 2 * xv[k][13] + 2 * xv[k][15] - xv[k][1] - xv[k][2] + xv[k][3] - xv[k][4] + xv[k][5] + xv[k][6] - xv[k][7]);
        a[k][19] = 0.125 * (xv[k][0] - 2 * xv[k][16] + 2 * xv[k][17] - 2 * xv[k][18] - xv[k][1] + 2 * xv[k][19] + xv[k][2] - xv[k][3] + xv[k][4] - xv[k][5] + xv[k][6] - xv[k][7]);
        a[k][20] = 0.25 * (xv[k][16] + xv[k][17] + xv[k][18] + xv[k][19] - 2 * (xv[k][20] + xv[k][21] + xv[k][22] + xv[k][23])) + xv[k][26];
        a[k][21] = 0.25 * (xv[k][9] + xv[k][11] + xv[k][13] + xv[k][15] - 2 * (xv[k][21] + xv[k][23] + xv[k][24] + xv[k][25])) + xv[k][26];
        a[k][22] = 0.25 * (xv[k][10] + xv[k][12] + xv[k][14] - 2 * (xv[k][20] + xv[k][22] + xv[k][24] + xv[k][25] - 2 * xv[k][26]) + xv[k][8]);
        a[k][23] = 0.125 * (- xv[k][0] - 2 * xv[k][10] + 2 * xv[k][12] - 2 * xv[k][14] + 2 * xv[k][16] + 2 * xv[k][17] - 2 * xv[k][18] - xv[k][1] -
                            2 * xv[k][19] - 4 * xv[k][20] + 4 * xv[k][22] + xv[k][2] + xv[k][3] - xv[k][4] - xv[k][5] + xv[k][6] + xv[k][7] + 2 * xv[k][8]);
        a[k][24] = 0.125 * (- xv[k][0] + 2 * xv[k][9] + 2 * xv[k][10] + 2 * xv[k][11] - 2 * xv[k][12] - 2 * xv[k][13] - 2 * xv[k][14] -
                            2 * xv[k][15] - xv[k][1] - 4 * xv[k][24] + 4 * xv[k][25] - xv[k][2] - xv[k][3] + xv[k][4] + xv[k][5] + xv[k][6] + xv[k][7] + 2 * xv[k][8]);
        a[k][25] = 0.125 * (- xv[k][0] - 2 * xv[k][9] + 2 * xv[k][11] - 2 * xv[k][13] + 2 * xv[k][15] + 2 * xv[k][16] - 2 * xv[k][17] -
                            2 * xv[k][18] + xv[k][1] + 2 * xv[k][19] + 4 * xv[k][21] - 4 * xv[k][23] + xv[k][2] - xv[k][3] - xv[k][4] + xv[k][5] + xv[k][6] - xv[k][7]);
        a[k][26] = 0.125 * (xv[k][0] - 2 * xv[k][9] - 2 * xv[k][10] - 2 * xv[k][11] - 2 * xv[k][12] - 2 * xv[k][13] - 2 * xv[k][14] - 2 * xv[k][15] -
                            2 * xv[k][16] - 2 * xv[k][17] - 2 * xv[k][18] + xv[k][1] - 2 * xv[k][19] + 4 * xv[k][20] + 4 * xv[k][21] + 4 * xv[k][22] + 4 * xv[k][23] +
                            4 * xv[k][24] + 4 * xv[k][25] - 8 * xv[k][26] + xv[k][2] + xv[k][3] + xv[k][4] + xv[k][5] + xv[k][6] + xv[k][7] - 2 * xv[k][8]);
      }
    }

    std::vector < double > phi(nDofs);

    double xi = x[0];
    double eta = x[1];
    double zita = x[2];

    phi[0] = 1.;
    phi[1] = xi;
    phi[2] = eta;
    phi[3] = zita;
    phi[4] = phi[1] * phi[2]; // x y
    phi[5] = phi[1] * phi[3]; // x z
    phi[6] = phi[2] * phi[3]; // y z
    phi[7] = phi[4] * phi[3];  // x y z
    if(solutionType > 0) {
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
    }
    if(solutionType > 1) {
      phi[20] = phi[7] * phi[8];  // xx yy
      phi[21] = phi[7] * phi[9];  // xx zz
      phi[22] = phi[8] * phi[9];  // yy zz
      phi[23] = phi[11] * phi[9]; // xx y zz
      phi[24] = phi[7] * phi[13]; // xx yy z
      phi[25] = phi[14] * phi[9]; // x yy zz
      phi[26] = phi[20] * phi[9]; // xx yy zz
    }

    std::vector < double > xp(dim);
    for(int i = 0; i < nDofs; i++) {
      for(int k = 0; k < dim; k++) {
        xp[k] += a[k][i] * phi[i];
      }
    }
    return xp;
  }


  std::vector< double > Marker::InverseMappingTet(const unsigned &currentElem, const unsigned &solutionType,
      std::vector< double > &x) {
    unsigned dim = 3;
    std::vector< std::vector < double > > xv(dim);
    std::vector < std::vector < double > > a(dim);

    unsigned nDofs = _mesh->GetElementDofNumber(currentElem, solutionType);
    short unsigned currentElementType = _mesh->GetElementType(currentElem);

    for(unsigned k = 0; k < dim; k++) {
      xv[k].resize(nDofs);
      a[k].resize(nDofs);
    }
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned iDof  = _mesh->GetSolutionDof(i, currentElem, 2);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
      }
    }

    if(solutionType == 0) {
      for(int k = 0; k < dim; k++) {
// 	xv[k][0] = 0;
//         xv[k][1] = 1;
//         xv[k][2] = 2;
//         xv[k][3] = 3;
	
        a[k][0] = xv[k][0] ;
        a[k][1] = - xv[k][0] + xv[k][1] ;
        a[k][2] = - xv[k][0] + xv[k][2] ;
        a[k][3] = - xv[k][0] + xv[k][3] ;
	
// 	if(k == 0) {
//           for(int j = 3; j >= 0; j--) {
//             std::cout << j << " " << a[k][j] << std::endl;
//           }
//         }
      }
    }
    else if(solutionType == 1) {
      
      for(int k = 0; k < dim; k++) {
	
// 	xv[k][0] = 0;
//         xv[k][1] = 1;
//         xv[k][2] = 2;
//         xv[k][3] = 3;
//         xv[k][4] = 4;
//         xv[k][5] = 5;
//         xv[k][6] = 6;
//         xv[k][7] = 7;
//         xv[k][8] = 8;
//         xv[k][9] = 9;
	
        a[k][0] = xv[k][0];
        a[k][1] =  - 3 * xv[k][0] - xv[k][1] + 4 * xv[k][4];
        a[k][2] = - 3 * xv[k][0] - xv[k][2] + 4 * xv[k][6];
        a[k][3] =  - 3 * xv[k][0] - xv[k][3] + 4 * xv[k][7];
        a[k][4] = 4 * xv[k][0] - 4 * xv[k][4] + 4 * xv[k][5] - 4 * xv[k][6];
        a[k][5] = 4 * xv[k][0] - 4 * xv[k][4] - 4 * xv[k][7] + 4 * xv[k][8];
        a[k][6] = 4 * xv[k][0] + 4 * xv[k][9] - 4 * xv[k][6] - 4 * xv[k][7];
        a[k][7] =  2 * xv[k][0] + 2 * xv[k][1] - 4 * xv[k][4];
        a[k][8] = 2 * xv[k][0] + 2 * xv[k][2] - 4 * xv[k][6];
        a[k][9] = 2 * xv[k][0] + 2 * xv[k][3] - 4 * xv[k][7];
	
// 	if(k == 0) {
//           for(int j = 9; j >= 0; j--) {
//             std::cout << j << " " << a[k][j] << std::endl;
//           }
//         }
      }
    }
    else if(solutionType == 2) {
      for(int k = 0; k < dim; k++) {

//         xv[k][0] = 0;
//         xv[k][1] = 1;
//         xv[k][2] = 2;
//         xv[k][3] = 3;
//         xv[k][4] = 4;
//         xv[k][5] = 5;
//         xv[k][6] = 6;
//         xv[k][7] = 7;
//         xv[k][8] = 8;
//         xv[k][9] = 9;
//         xv[k][10] = 10;
//         xv[k][11] = 11;
//         xv[k][12] = 12;
//         xv[k][13] = 13;
//         xv[k][14] = 14;


        a[k][14] = 4 * (xv[k][0] - 8 * xv[k][9] + 27 * xv[k][10] + 27 * xv[k][11] + 27 * xv[k][12] + 27 * xv[k][13] - 64 * xv[k][14] +
                        xv[k][1] + xv[k][2] + xv[k][3] - 8 * (xv[k][4] + xv[k][5] + xv[k][6] + xv[k][7] + xv[k][8]));
        a[k][13] = -3 * (xv[k][0] - 4 * xv[k][9] + 9 * xv[k][13] + xv[k][2] + xv[k][3] - 4 * (xv[k][6] + xv[k][7]));
        a[k][12] = -3  * (xv[k][0] + 9 * xv[k][11] + xv[k][1] + xv[k][3] - 4 * (xv[k][4] + xv[k][7] + xv[k][8]));
        a[k][11] = -3 * (xv[k][0] + 9  * xv[k][10] + xv[k][1] + xv[k][2] - 4 * (xv[k][4] + xv[k][5] + xv[k][6]));
        a[k][10] = -13 * xv[k][0] + 32  * xv[k][9] - 135 * xv[k][10] - 135 * xv[k][11] - 81  * xv[k][12] - 135 * xv[k][13] +
		    256 * xv[k][14] - 7 * xv[k][1] - 7 * xv[k][2] - 7 * xv[k][3] + 8 * (7 * xv[k][4] + 4 * xv[k][5] + 7 * (xv[k][6] + xv[k][7]) + 4 * xv[k][8]);
        a[k][9] = 2 * (xv[k][0] + xv[k][3] - 2 * xv[k][7]);
        a[k][8] = 2 * (xv[k][0] + xv[k][2] - 2 * xv[k][6]);
        a[k][7] = 2 * (xv[k][0] + xv[k][1] - 2 * xv[k][4]);
        a[k][6] = 7 * xv[k][0] - 8 * xv[k][9] + 27 * xv[k][13] + 3 * xv[k][2] + 3 * xv[k][3] - 16 * (xv[k][6] + xv[k][7]);
        a[k][5] = 7 * xv[k][0] + 27 * xv[k][11] + 3 * xv[k][1] + 3 * xv[k][3] - 8 * (2 * (xv[k][4] + xv[k][7]) + xv[k][8]);
        a[k][4] = 7 * xv[k][0] + 27 * xv[k][10] + 3 * xv[k][1] + 3 * xv[k][2] - 8 * (2 * xv[k][4] + xv[k][5] + 2 * xv[k][6]);
        a[k][3] = -3 * xv[k][0] - xv[k][3] + 4 * xv[k][7];
        a[k][2] = -3 * xv[k][0] - xv[k][2] + 4 * xv[k][6];
        a[k][1] = -3 * xv[k][0] - xv[k][1] + 4 * xv[k][4];
        a[k][0] = xv[k][0];

//         if(k == 0) {
//           for(int j = 14; j >= 0; j--) {
//             std::cout << j << " " << a[k][j] << std::endl;
//           }
//         }
      }
    }

    std::vector < double > phi(nDofs);

    double xi = x[0];
    double eta = x[1];
    double zita = x[2];

    phi[0] = 1.;
    phi[1] = xi; // x
    phi[2] = eta;  // y
    phi[3] = zita; // z
    if(solutionType > 0) {
      phi[4] = phi[1] * phi[2];  // x y
      phi[5] = phi[1] * phi[3];  // x z
      phi[6] = phi[2] * phi[3];  // y z
      phi[7] = phi[1] * phi[1]; // x x
      phi[8] = phi[2] * phi[2]; // y y
      phi[9] = phi[3] * phi[3]; // z z
    }
    if(solutionType > 1) {
      phi[10] = phi[4] * phi[3]; // x y z
      phi[11] = phi[1] * phi[8] + phi[7] * phi[2] ; // x y y + x x y
      phi[12] = phi[1] * phi[9] + phi[7] * phi[3] ; // x z z + x x z;
      phi[13] = phi[2] * phi[9] + phi[8] * phi[3] ; // y z z + y y z;
      phi[14] = phi[4] * phi[9] + phi[4] * phi[6] + phi[7] * phi[6]; // x y z z + x y y z + x x y z
    }

    std::vector < double > xp(dim);
    for(int i = 0; i < nDofs; i++) {
      for(int k = 0; k < dim; k++) {
        xp[k] += a[k][i] * phi[i];
      }
    }
    return xp;
  }

  std::vector< double > Marker::InverseMappingWedge(const unsigned &currentElem, const unsigned &solutionType,
      std::vector< double > &x) {
    unsigned dim = 3;
    std::vector< std::vector < double > > xv(dim);
    std::vector < std::vector < double > > a(dim);

    unsigned nDofs = _mesh->GetElementDofNumber(currentElem, solutionType);
    short unsigned currentElementType = _mesh->GetElementType(currentElem);

    for(unsigned k = 0; k < dim; k++) {
      xv[k].resize(nDofs);
      a[k].resize(nDofs);
    }
    
     for(unsigned i = 0; i < nDofs; i++) {
        unsigned iDof  = _mesh->GetSolutionDof(i, currentElem, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
	  xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
	}
     }
    
    
    
//     for(unsigned i = 0; i < nDofs; i++) {
//       unsigned iDof  = _mesh->GetSolutionDof(i, currentElem, 2);    // global to global mapping between coordinates node and coordinate dof
//       for(unsigned k = 0; k < dim; k++) {
//         xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates local storage for the element coordinates
//       }
//     }

    if(solutionType == 0) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = 0.5 * (xv[k][0] + xv[k][3]);
        a[k][1] = 0.5 * (- xv[k][0] + xv[k][1] - xv[k][3] + xv[k][4]);
        a[k][2] = 0.5 * (- xv[k][0] + xv[k][2] - xv[k][3] + xv[k][5]);
        a[k][3] = 0.5 * (- xv[k][0] + xv[k][3]);
        a[k][4] = 0.5 * (xv[k][0] - xv[k][1] - xv[k][3] + xv[k][4]);
        a[k][5] = 0.5 * (xv[k][0] - xv[k][2] - xv[k][3] + xv[k][5]);
      }
    }
    else if(solutionType == 1) {
      for(int k = 0; k < dim; k++) {
        a[k][0] = xv[k][12];
        a[k][1] = - 3 * xv[k][12] - xv[k][13];
        a[k][2] = - 3 * xv[k][12] - xv[k][14];
        a[k][3] = 0.5 * (- xv[k][0] + xv[k][3]);
        a[k][4] = 4 * xv[k][12];
        a[k][5] = 0.5 * (3 * xv[k][0] + 4 * xv[k][9] + xv[k][1] - 3 * xv[k][3] - xv[k][4] - 4 * xv[k][6]);
        a[k][6] = 0.5 * (3 * xv[k][0] + 4 * xv[k][11] + xv[k][2] - 3 * xv[k][3] - xv[k][5] - 4 * xv[k][8]);
        a[k][7] = 2 * (xv[k][12] + xv[k][13]);
        a[k][8] = 2 * (xv[k][12] + xv[k][14]);
        a[k][9] = 0.5 * (xv[k][0] - 2 * xv[k][12] + xv[k][3]);
        a[k][10] = - 2 * (xv[k][0] + xv[k][9] - xv[k][10] + xv[k][11] - xv[k][3] - xv[k][6] + xv[k][7] - xv[k][8]);
        a[k][11] = - xv[k][0] - 2 * xv[k][9] - xv[k][1] + xv[k][3] + xv[k][4] + 2 * xv[k][6];
        a[k][12] = - xv[k][0] - 2 * xv[k][11] - xv[k][2] + xv[k][3] + xv[k][5] + 2 * xv[k][8];
        a[k][13] = 0.5 * (- 3 * xv[k][0] + 4 * xv[k][9] + 6 * xv[k][12] + 2 * xv[k][13] - xv[k][1] - 3 * xv[k][3] - xv[k][4] + 4 * xv[k][6]);
        a[k][14] = 0.5 * (- 3 * xv[k][0] + 4 * xv[k][11] + 6 * xv[k][12] + 2 * xv[k][14] - xv[k][2] - 3 * xv[k][3] - xv[k][5] + 4 * xv[k][8]);
        a[k][15] = xv[k][0] - 2 * xv[k][9] - 2 * xv[k][12] - 2 * xv[k][13] + xv[k][1] + xv[k][3] + xv[k][4] - 2 * xv[k][6];
        a[k][16] = xv[k][0] - 2 * xv[k][11] - 2 * xv[k][12] - 2 * xv[k][14] + xv[k][2] + xv[k][3] + xv[k][5] - 2 * xv[k][8];
        a[k][17] = 2 * (xv[k][0] - xv[k][9] + xv[k][10] - xv[k][11] - 2 * xv[k][12] + xv[k][3] - xv[k][6] + xv[k][7] - xv[k][8]);
      }
    }
    else if(solutionType == 2) {
      for(int k = 0; k < dim; k++) {
	 
 	/*for(int i=0 ; i <nDofs;  i++){
 	  xv[k][i] = i;
 	} */ 
	    
            a[k][0] = xv[k][12];
	    a[k][1] = - 3 * xv[k][12] - xv[k][13] + 4 * xv[k][15];
	    a[k][2] = - 3 * xv[k][12] - xv[k][14] + 4 * xv[k][17];
	    a[k][3] = 0.5 * (- xv[k][0] + xv[k][3]);
            a[k][4] = 7 * xv[k][12] + 3 * xv[k][13] + 3 * xv[k][14] - 8 * (2 * xv[k][15] + xv[k][16] + 2 * xv[k][17]) + 27 * xv[k][20];
	    a[k][5] = 0.5 * (3 * xv[k][0] + 4 * xv[k][9] + xv[k][1] - 3 * xv[k][3] - xv[k][4] - 4 * xv[k][6]);
	    a[k][6] = 0.5 * (3 * xv[k][0] + 4 * xv[k][11] + xv[k][2] - 3 * xv[k][3] - xv[k][5] - 4 * xv[k][8]);
	    a[k][7] = 2 * (xv[k][12] + xv[k][13] - 2 * xv[k][15]);
	    a[k][8] = 2 * (xv[k][12] + xv[k][14] - 2 * xv[k][17]);
	    a[k][9] = 0.5 * (xv[k][0] - 2 * xv[k][12] + xv[k][3]);
	    a[k][10] = 0.5 * (- 7 * xv[k][0] - 16 * xv[k][9] - 8 * xv[k][10] - 16 * xv[k][11] - 27 * xv[k][18] -
                              3 * xv[k][1] + 27 * xv[k][19] - 3 * xv[k][2] + 7 * xv[k][3] + 3 * xv[k][4] + 3 * xv[k][5] +
                              8 * (2 * xv[k][6] + xv[k][7] + 2 * xv[k][8])) ;
	    a[k][11] = - 3 * (xv[k][12] + xv[k][13] + xv[k][14] - 4 * (xv[k][15] + xv[k][16] + xv[k][17]) + 9 * xv[k][20]);
	    a[k][12] = - xv[k][0] - 2 * xv[k][9] - xv[k][1] + xv[k][3] + xv[k][4] + 2 * xv[k][6];
	    a[k][13] = - xv[k][0] - 2 * xv[k][11] - xv[k][2] + xv[k][3] + xv[k][5] + 2 * xv[k][8];
	    a[k][14] = 0.5 * (3 * xv[k][0] + 4 * xv[k][9] + 6 * xv[k][12] + 2 * xv[k][13] - 8 * xv[k][15] - xv[k][1] - 3 * xv[k][3] - xv[k][4] + 4 * xv[k][6]);
	    a[k][15] = 0.5 * (-3 * xv[k][0] + 4 * xv[k][11] + 6 * xv[k][12] + 2 * xv[k][14] - 8 * xv[k][17] - xv[k][2] - 3 * xv[k][3] - xv[k][5] + 4 * xv[k][8]);
	    a[k][16] = xv[k][0] - 2 * xv[k][9] - 2 * xv[k][12] - 2 * xv[k][13] + 4 * xv[k][15] + xv[k][1] + xv[k][3] + xv[k][4] - 2 * xv[k][6];
	    a[k][17] = xv[k][0] - 2 * xv[k][11] - 2 * xv[k][12] - 2 * xv[k][14] + 4 * xv[k][17] + xv[k][2] + xv[k][3] + xv[k][5] - 2 * xv[k][8];
	    a[k][18] = 1.5 * (xv[k][0] + 4 * xv[k][9] + 4 * xv[k][10] + 4 * xv[k][11] + 9 * xv[k][18] + xv[k][1] - 9 * xv[k][19] + xv[k][2] - xv[k][3] - xv[k][4] - xv[k][5] - 4 * (xv[k][6] + xv[k][7] + xv[k][8]));
	    a[k][19] = 0.5 * (7 * xv[k][0] - 16 * xv[k][9] - 8 * xv[k][10] - 16 * xv[k][11] - 14 * xv[k][12] - 6 * xv[k][13] - 6 * xv[k][14] + 32 * xv[k][15] + 
                              16 * xv[k][16] + 32 * xv[k][17] + 27 * xv[k][18] + 3 * xv[k][1] + 27 * xv[k][19] - 54 * xv[k][20] + 3 * xv[k][2] + 7 * xv[k][3] + 
                              3 * xv[k][4] + 3 * xv[k][5] - 8 * (2 * xv[k][6] + xv[k][7] + 2 * xv[k][8]));
	    a[k][20] = -1.5 * (xv[k][0] - 4 * xv[k][9] - 4 * xv[k][10] - 4 * xv[k][11] - 2 * xv[k][12] - 2 * xv[k][13] - 2 * xv[k][14] + 8 * xv[k][15] + 8 * xv[k][16] + 8 * xv[k][17] + 9 * xv[k][18] + xv[k][1] + 9 * xv[k][19] - 18 * xv[k][20] + 
	                       xv[k][2] + xv[k][3] + xv[k][4] + xv[k][5] - 4 * (xv[k][6] + xv[k][7] + xv[k][8]));
	            if(k == 0) {
//           for(int j = 20; j >= 0; j--) {
//             std::cout << j << " " << a[k][j] << std::endl;
//           }
        }
      }
    }

    std::vector < double > phi(nDofs);

    double xi = x[0];
    double eta = x[1];
    double zita = x[2];

    phi[0] = 1.;
    phi[1] = xi; // x
    phi[2] = eta;  // y
    phi[3] = zita; // z
//     phi[4] = phi[1] * phi[3];  // x z
//     phi[5] = phi[2] * phi[3];  // y z
    if(solutionType > 0) {
      phi[4] = xi * eta;  // x y
      phi[5] = xi * zita;  // x z
      phi[6] = eta * zita;  // y z
      phi[7] = xi * xi; // x x
      phi[8] = eta * eta; // y y
      phi[9] = zita * zita; // z z
      phi[10] = phi[4] * zita; // x y z
//       phi[11] = phi[8] * phi[3]; // x x z
//       phi[12] = phi[9] * phi[3]; // y y z
//       phi[13] = phi[9] * phi[3]; // y y z
//       phi[14] = phi[2] * phi[10]; // y z z;
//       phi[15] = phi[8] * phi[10]; // x x z z
//       phi[16] = phi[9] * phi[10]; // y y z z
//       phi[17] = phi[5] * phi[10]; // x y z z;
    }
    if(solutionType > 1) {
      
      std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "<<xi << " "<<eta<<" "<< zita<<std::endl;
      phi[11] = phi[7] * eta + xi * phi[8]; // xx y + x yy
      phi[12] = phi[7] * zita; // x x z
      phi[13] = phi[8] * zita; // y y z
      phi[14] = xi  * phi[9]; // x z z
      phi[15] = eta * phi[9]; // y z z
      phi[16] = phi[7] * phi[9]; // xx zz
      phi[17] = phi[8] * phi[9]; // yy zz
      phi[18] = phi[11] * zita; // x yy z + xx y z
      phi[19] = phi[10] * zita; // x y zz
      phi[20] = phi[18] * zita; // x yy zz + xx y zz
    }

    std::vector < double > xp(dim);
    for(int i = 0; i < nDofs; i++) {
      for(int k = 0; k < dim; k++) {
        xp[k] += a[k][i] * phi[i];
      }
    }
    return xp;
  }


  void Marker::InverseMappingTEST(std::vector < double > &x) {
    unsigned dim = _mesh->GetDimension();
    std::vector < std::vector < double > > xv(dim);
    for(int solType = 0; solType < 3; solType+= 2) {
      std::cout << "solType = " << solType << std::endl;
      int iel = _mesh->_elementOffset[_iproc + 1] - 1 ;

      unsigned nDofs = _mesh->GetElementDofNumber(iel, solType);
      short unsigned ielType = _mesh->GetElementType(iel);


      for(int j = 0; j < nDofs; j++) {
        for(unsigned k = 0; k < dim; k++) {
          x[k] = *(_mesh->_finiteElement[ielType][0]->GetBasis()->GetXcoarse(j) + k);
          xv[k].resize(nDofs);
        }
        for(unsigned i = 0; i < nDofs; i++) {
          unsigned iDof  = _mesh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
          }
        }

        // This is the test
        std::vector < double > xp;
        if(ielType == 0) {
          xp = InverseMappingHex(iel, solType, x);
        }
        else if(ielType == 1) {
          xp = InverseMappingTet(iel, solType, x);
        }
        else if(ielType == 2) {
          xp = InverseMappingWedge(iel, solType, x);
        }
        else if(ielType == 3) {
          xp = InverseMappingQuad(iel, solType, x);
        }
        else if(ielType == 4) {
          xp = InverseMappingTri(iel, solType, x);
        }

	std::cout << j << " ";
        for(int k = 0; k < dim; k++) {
	  double a = ( fabs ( xp[k] - j)< 1.0e-10 ) ? 0 : xp[k] - xv[k][j];
	  std::cout << a << " ";
	  //std::cout << xp[k] <<" ";
	  
        }
	std::cout << std::endl;
      }

    }
  }

}




