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
#include <math.h>

#include "PolynomialBases.hpp"

const unsigned facePointNumber[6] = {0, 0, 0, 9, 7, 0};


const unsigned facePoints[6][9] = {
  { },
  { },
  { },
  { 0, 4, 1, 5, 2, 6, 3, 7, 0},
  { 0, 3, 1, 4, 2, 5, 0},
  { }
};

// const unsigned facePoints3D[6][6][9] = {
//     {   {0, 8, 1, 17, 5, 12, 4, 16, 0},
//         {1, 9, 2, 18, 6, 13, 5, 17, 1},
//         {2, 10, 3, 19, 7, 14, 6, 18, 2},
//         {3, 11, 0, 16, 4, 15, 7, 19, 3},
//         {0, 11, 3, 10, 2, 9, 1, 8, 0},
//         {4, 12, 5, 13, 6, 14, 7, 15, 4},
//     },
//     {{}}, {{}}, {{}}, {{}}, {{}}
// };

const unsigned trianglesPerFace[6][6] = { //type - number of faces
  {8, 8, 8, 8, 8, 8},
  {6, 6, 6, 6},
  {8, 8, 8, 6, 6},
  {8},
  {6},
  {}
};

const unsigned faceTriangleNodes[6][6][8][4] = { // type - number of faces - number of triangles per faces - 3 vertices of the triangle (the first is counted twice)

  { {{0, 8, 20, 0}, {8, 1, 20, 8}, {1, 17, 20, 1}, {17, 5, 20, 17}, {5, 12, 20, 5}, {12, 4, 20, 12}, {4, 16, 20, 4 }, {16, 0, 20, 16}},
    {{1, 9, 21, 1}, {9, 2, 21, 9}, {2, 18, 21, 2}, {18, 6, 21, 18}, {6, 13, 21, 6}, {13, 5, 21, 13}, {5, 17, 21, 5}, {17, 1, 21, 17}},
    {{2, 10, 22, 2}, {10, 3, 22, 10}, {3, 19, 22, 3}, {19, 7, 22, 19}, {7, 14, 22, 7}, {14, 6, 22, 14}, {6, 18, 22, 6}, {18, 2, 22, 18}},
    {{3, 11, 23, 3}, {11, 0, 23, 11}, {0, 16, 23, 0}, {16, 4, 23, 16}, {4, 15, 23, 4}, {15, 7, 23, 15}, {7, 19, 23, 7}, {19, 3, 23, 19}},
    {{0, 11, 24, 0}, {11, 3, 24, 11}, {3, 10, 24, 3}, {10, 2, 24, 10}, {2, 9, 24, 2}, {9, 1, 24, 9}, {1, 8, 24, 1}, {8, 0, 24, 8}},
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
  {{{0, 4, 8, 0}, {4, 1, 8, 4}, {1, 5, 8, 1}, {5, 2, 8, 5}, {2, 6, 8, 2}, {6, 3, 8, 6}, {3, 7, 8, 3}, {7, 0, 8, 7}}},
  {{{0, 3, 6, 0}, {3, 1, 6, 3}, {1, 4, 6, 1}, {4, 2, 6, 4}, {2, 5, 6, 2}, {5, 0, 6, 5}}},
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
          nextElem[_iproc] = GetNextElement2D(iel, previousElem[_iproc]);
        }
        else if(dim == 3) {
          nextElem[_iproc] = GetNextElement3D(iel, previousElem[_iproc]);
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






  unsigned Marker::GetNextElement2D(const unsigned &currentElem, const unsigned &previousElem) {


    unsigned dim = _mesh->GetDimension();
    int nextElem ;
    bool markerIsInElement = false;
    bool nextElementFound = false;
    short unsigned currentElementType = _mesh->GetElementType(currentElem);
    double epsilon  = 10.e-10;
    double epsilon2  = epsilon * epsilon;
    double t;

    std::vector<double> xc(dim, 0); //stores the coordinates of the face node of currentElem
    unsigned faceNodeLocalIndex;

    if(currentElementType == 3) faceNodeLocalIndex = 8;
    else if(currentElementType == 4) faceNodeLocalIndex = 6;

    unsigned faceNodeDof = _mesh->GetSolutionDof(faceNodeLocalIndex, currentElem, 2);
    //std::cout << "faceNodeDof = " << faceNodeDof << std::endl;

    for(unsigned k = 0; k < dim; k++) {
      xc[k] = (*_mesh->_topology->_Sol[k])(faceNodeDof) - _x[k];    // coordinates are translated so that the marker is the new origin
    }

    if(xc[0]*xc[0] < epsilon2 && xc[1]*xc[1] < epsilon2) {
      std::cout << "the marker is the central face node" << std::endl;
      markerIsInElement = true; //the marker is xc
    }

    else {

      std::vector<double> xcc(dim, 0); // will store the coordinates of the center scaled

      std::vector<double> r(dim, 0);   //coordinates of the intersection point between the line of the edges and the line that connects the marker and the face node
      std::vector< std::vector < double > > xv(dim);   //stores the coordinates of the vertices and midpoints of the element, the first and the last are the same

      for(unsigned k = 0; k < dim; k++) {
        xv[k].reserve(9);
      }

      for(unsigned k = 0; k < dim; k++) {
        xv[k].resize(facePointNumber[currentElementType]);
      }

      for(unsigned i = 0; i < facePointNumber[currentElementType]; i++) {
        unsigned inodeDof  = _mesh->GetSolutionDof(facePoints[currentElementType][i], currentElem, 2);
        std::cout << "inodeDof = " << inodeDof << std::endl;

        for(unsigned k = 0; k < dim; k++) {
          xv[k][i] = (*_mesh->_topology->_Sol[k])(inodeDof) - _x[k];
        }
      }

      // rescaling coordinates to properly handle different scales of meshes
      double length = 0.;
      double sum = 0.;

      for(unsigned i = 0; i < facePointNumber[currentElementType] - 1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          sum += (xv[k][i + 1] - xv[k][i]) * (xv[k][i + 1] - xv[k][i]);
        }

        length += sqrt(sum);
      }

      length /= facePointNumber[currentElementType];
      std::cout << "length= " << length << std::endl;

      for(unsigned k = 0; k < dim; k++) {
        xcc[k] = xc[k] / length;

        for(unsigned i = 0; i < facePointNumber[currentElementType]; i++) {
          xv[k][i] /= length;
        }
      }


      for(unsigned i = 0 ; i < facePointNumber[currentElementType] - 1; i++) {

        // let's find the plane passing through the points xv[][i], xv[][i+1] and xv[][2] = xv[][i] but with z = length .
        double A = (xv[1][i + 1] - xv[1][i]);
        double B = -(xv[0][i + 1] - xv[0][i]);

        //std::cout << "A= " << A << " , " <<"B= " << B <<std::endl;

        double tBottom = (A * xcc[0] + B * xcc[1]) ;
        double tTop = A * xv[0][i] + B * xv[1][i];
        std::cout << "tBottom = " << tBottom << " , " << "A= " << A << " , " <<  "B= " << B << " , " << "xv[1][" << i << "] =" << xv[1][i] << " , " <<  "tTop = " <<   tTop << std::endl;

        if(fabs(tBottom) < epsilon && tTop != 0) {
          // std::cout << "The plane of edge " << i << "does not intersect the line" <<std::endl;
        }

        else {
          //now let's find the coordinates of the intersection point r
          t = tTop / tBottom ;
          std::cout << "t = " << t << std::endl;

          for(unsigned k = 0; k < dim; k++) {
            r[k] = t * xcc[k];
            //std::cout << "r[" << k << "] = " << r[k] <<std::endl;
          }

          if(t < 1) {   //if not, it means the point r is far away from the marker, and we don't want to go in that direction

            std::vector< std::vector < double > > xvr(dim);

            for(unsigned k = 0; k < dim; k++) {
              xvr[k].reserve(9);
            }

            for(unsigned k = 0; k < dim; k++) {
              xvr[k].resize(facePointNumber[currentElementType]);
            }

            //now we have to determine if r is inside edge i
            for(unsigned j = 0; j < facePointNumber[currentElementType]; j++) {
              for(unsigned k = 0; k < dim; k++) {
                xvr[k][j] = xv[k][j] - r[k];     //transate again the reference frame so that the origin is r
              }
            }


            if((xvr[0][i] * xvr[0][i]  + xvr[1][i] * xvr[1][i]) < epsilon2 ||
                (xvr[0][i + 1]*xvr[0][i + 1] + xvr[1][i + 1]*xvr[1][i + 1]) < epsilon2) {
              std::cout << "intersection on a vertex of the edge" << std::endl;

              if(fabs(t) < epsilon || t < 0) {   //this means the marker is on one of the edges

                if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is one of the nodes" << std::endl;

                if(t < 0) std::cout << "setting markerIsInElement = true because r is one of the nodes" << std::endl;

                markerIsInElement = true;
                break;
              }
              else {
                unsigned nodeIndex;

                if(i % 2 == 0 && i != facePointNumber[currentElementType]) nodeIndex = i / 2 ;
                else if(i == facePointNumber[currentElementType]) nodeIndex = (i - 2) / 2 ;
                else if(i % 2 != 0) nodeIndex = (i - 1) / 2 ;

                nextElem = (_mesh->el->GetFaceElementIndex(currentElem, nodeIndex) - 1);
                nextElementFound = true;
                break;
              }

            }


            else if(xvr[0][i]*xvr[0][i + 1] < 0 || xvr[1][i]*xvr[1][i + 1] < 0) {
              std::cout << "intersection on an edge" << std::endl;

              if(fabs(t) < epsilon || t < 0) {   //this means the marker is on one of the edges

                if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges " << std::endl;

                if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges " << std::endl;

                markerIsInElement = true;
                break;
              }
              else {
                unsigned nodeIndex;

                if(i % 2 == 0 && i != facePointNumber[currentElementType]) nodeIndex = i / 2 ;
                else if(i == facePointNumber[currentElementType]) nodeIndex = (i - 2) / 2 ;
                else if(i % 2 != 0) nodeIndex = (i - 1) / 2 ;

                nextElem = (_mesh->el->GetFaceElementIndex(currentElem, nodeIndex) - 1);
                nextElementFound = true;
                break;
              }
            }
          } // closes the " if t < 1 "
        } // closes the else
      } //closes the for on the nodes
    }// closes the else before the for

    if(markerIsInElement == true) {
      nextElem = currentElem;
      std::cout << "The marker belongs to element " << currentElem << std::endl;
    }

    if(nextElementFound == true) {
      std::cout << "The marker does not belong to element " << currentElem << std::endl;
    }


    std::cout << "markerIsInElement = " << markerIsInElement << " , " << "nextElementFound= " << nextElementFound << ", " << "nextElem = " << nextElem << std::endl;

    return (nextElem >= 0) ? nextElem : UINT_MAX;

  }





  unsigned Marker::GetNextElement3D(const unsigned &currentElem, const unsigned &previousElem) {


    unsigned dim = _mesh->GetDimension();
    int nextElem ;
    bool markerIsInElement = false;
    bool nextElementFound = false;
    short unsigned currentElementType = _mesh->GetElementType(currentElem);
    double epsilon  = 10.e-10;
    double epsilon2  = epsilon * epsilon;
    double t;

    std::vector<double> xc(dim, 0); //stores the coordinates of the central node of currentElem
    unsigned centralNodeLocalIndex;

    if(currentElementType == 0) centralNodeLocalIndex = 26;
    else if(currentElementType == 1) centralNodeLocalIndex = 14;
    else if(currentElementType == 2) centralNodeLocalIndex = 20;

    unsigned centralNodeDof = _mesh->GetSolutionDof(centralNodeLocalIndex, currentElem, 2);

    for(unsigned k = 0; k < dim; k++) {
      xc[k] = (*_mesh->_topology->_Sol[k])(centralNodeDof) - _x[k];    // coordinates are translated so that the marker is the new origin
      //std::cout << " xc[" << k << "]= " <<xc[k] <<std::endl;
    }

    if(xc[0]*xc[0] < epsilon2 && xc[1]*xc[1] < epsilon2 && xc[2]*xc[2] < epsilon2) {
      std::cout << "the marker is the central element node" << std::endl;
      markerIsInElement = true; //the marker is xc
    }

    else {


      for(unsigned iface = 0; iface < _mesh->GetElementFaceNumber(currentElem); iface++) {

        std::cout << "iface = " << iface << std::endl;

        for(unsigned itri = 0; itri < trianglesPerFace[currentElementType][iface]; itri ++) {

          std::cout << "itri = " << itri << std::endl;

          std::vector<double> xcc(dim, 0); // will store the coordinates of the center scaled

          unsigned scalarCount = 0;
          std::vector<double> r(dim, 0);   //coordinates of the intersection point between the plane of itri and the line through the element center point and the marker
          std::vector< std::vector < double > > xv(dim);   //stores the coordinates of the nodes of the triangle itri

          // fill in the coordinates of the vertices of itri
          for(unsigned k = 0; k < dim; k++) {
            xv[k].reserve(4);
          }
          for(unsigned k = 0; k < dim; k++) {
            xv[k].resize(4);
          }

          for(unsigned i = 0; i < 4; i++) {
            unsigned itriDof  = _mesh->GetSolutionDof(faceTriangleNodes[currentElementType][iface][itri][i], currentElem, 2);
            // std::cout << "itriDof = " << itriDof << std::endl;
            for(unsigned k = 0; k < dim; k++) {
              xv[k][i] = (*_mesh->_topology->_Sol[k])(itriDof) - _x[k];     // coordinates are translated so that the marker is the new origin
            }
          }

          // rescaling coordinates to properly handle different scales of meshes
          double length = 0.;
          double sum = 0.;
          for(unsigned i = 0; i < 3; i++) {
            for(unsigned k = 0; k < dim; k++) {
              sum += (xv[k][i + 1] - xv[k][i]) * (xv[k][i + 1] - xv[k][i]);
            }
            length += sqrt(sum);
          }

          length /= 4;

          for(unsigned k = 0; k < dim; k++) {
            xcc[k] = xc[k] / length;
            //std::cout << " xcc[" << k << "]= " <<xcc[k] <<std::endl;
            for(unsigned i = 0; i < 4; i++) {
              xv[k][i] /= length;
            }
          }

          // let's find the plane passing through the vertices of the triangle itri
          double A = -(xv[1][2] - xv[1][0]) * (xv[2][1] - xv[2][0]) + (xv[2][2] - xv[2][0]) * (xv[1][1] - xv[1][0]);
          double B = -(xv[2][2] - xv[2][0]) * (xv[0][1] - xv[0][0]) + (xv[0][2] - xv[0][0]) * (xv[2][1] - xv[2][0]);
          double C = -(xv[0][2] - xv[0][0]) * (xv[1][1] - xv[1][0]) + (xv[1][2] - xv[1][0]) * (xv[0][1] - xv[0][0]);

          //std::cout << "A= " << A << " , " <<"B= " << B << " , " << "C = " << C << " , " <<std::endl;

          double tBottom = (A * xcc[0] + B * xcc[1] + C * xcc[2]);
          double tTop = A * xv[0][0] + B * xv[1][0] + C * xv[2][0];

          //std::cout << " tTop = " << tTop <<std::endl;
          //std::cout << " tBottom = " << tBottom <<std::endl;

          if(fabs(tBottom) < epsilon && tTop != 0) {
            // std::cout << "The plane of face" << itri << "does not intersect the line" <<std::endl;
            break; // must exit the loop on itri
          }

          else {
            //now let's find the coordinates of the intersection point r
            t = tTop / tBottom ;
            std::cout << "t = " << t << std::endl;

            for(unsigned k = 0; k < dim; k++) {
              r[k] = t * xcc[k];
              // std::cout << "r[" << k << "] = " << r[k] <<std::endl;
            }

            if(t < 1) {   //if not, it means the point r is far away from the marker, and we don't want to go in that direction

              //now we have to determine if r is inside itri
              for(unsigned i = 0; i < 4; i++) {
                for(unsigned k = 0; k < dim; k++) {
                  xv[k][i] = xv[k][i] - r[k];     //transate again the reference frame so that the origin is r
                }
              }

              for(unsigned i = 0; i < 3; i++) {
                double q0 = xv[1][i] * (xv[2][i] - xv[2][i + 1]) + xv[2][i] * (xv[1][i + 1] - xv[1][i]);
                double q1 = xv[2][i] * (xv[0][i] - xv[0][i + 1]) + xv[0][i] * (xv[2][i + 1] - xv[2][i]);
                double q2 = xv[0][i] * (xv[1][i] - xv[1][i + 1]) + xv[1][i] * (xv[0][i + 1] - xv[0][i]);

                // std::cout << "q0 = " << q0 << " , " << "q1 = " << q1 << " , " << " q2 = " << q2 <<  std::endl;

                double  scalarProduct = q0 * A + q1 * B + q2 * C;

                std::cout << "fabs(scalarProduct) = " << fabs(scalarProduct) << std::endl;

                if(scalarProduct > epsilon) {
                  std::cout << "r is outside triangle " << itri <<  std::endl;
                  break;

                }
                else if(fabs(scalarProduct) < epsilon) {  //scalarProduct == 0

                  if((xv[0][i] * xv[0][i]  + xv[1][i] * xv[1][i] + xv[2][i] * xv[2][i]) < epsilon2 ||
                      (xv[0][i + 1]*xv[0][i + 1] + xv[1][i + 1]*xv[1][i + 1] + xv[2][i + 1]*xv[2][i + 1]) < epsilon2) {
                    std::cout << "intersection on a vertex of itri" << std::endl;
                    if(fabs(t) < epsilon || t < 0) {   //this means the marker is on one of the faces

                      if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is one vertex of triangle " << itri << std::endl;
                      if(t < 0) std::cout << "setting markerIsInElement = true because r is one vertex of triangle " << itri << std::endl;

                      markerIsInElement = true;
                      break;
                    }
                    else {
                      std::cout << "r is in triangle " << itri << std::endl;
                      nextElem = (_mesh->el->GetFaceElementIndex(currentElem, iface) - 1);
                      nextElementFound = true;
                      break;
                    }

                  }


                  else if(xv[0][i]*xv[0][i + 1] < 0 || xv[1][i]*xv[1][i + 1] < 0 || xv[2][i]*xv[2][i + 1] < 0) {
                    std::cout << "intersection on an edge of itri" << std::endl;
                    if(fabs(t) < epsilon || t < 0) {   //this means the marker is on one of the faces

                      if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges of triangle " << itri << std::endl;
                      if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges of triangle " << itri << std::endl;

                      markerIsInElement = true;
                      break;
                    }
                    else {
                      std::cout << "r is in triangle " << itri << std::endl;
                      nextElem = (_mesh->el->GetFaceElementIndex(currentElem, iface) - 1);
                      nextElementFound = true;
                      break;
                    }
                  }
                }
                else if(scalarProduct < 0) {
                  std::cout << " scalarProduct = " << scalarProduct << std::endl;
                  scalarCount++;
                }
              } // closes the for loop
            } // closes " if t < 1 "
          } // closes the "else" on tBottom = 0


          if(scalarCount == 3) {
            if(fabs(t) < epsilon || t < 0) {   //this means the marker is on one of the faces

              if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges of triangle " << itri << std::endl;
              if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges of triangle " << itri << std::endl;

              markerIsInElement = true;
              break;
            }
            else {
              std::cout << "r is in triangle " << itri << std::endl;
              nextElem = (_mesh->el->GetFaceElementIndex(currentElem, iface) - 1);
              nextElementFound = true;
              break;
            }
          }

          if(nextElementFound == true) {
            break;
          }

          if(markerIsInElement == true) {
            break;
          }
        } //end for on itri

        if(nextElementFound == true) {
          break;
        }

        if(markerIsInElement == true) {
          break;
        }
      } //end for on iface
    } //end of else before for on iface

    if(markerIsInElement == true) {
      nextElem = currentElem;
      std::cout << "The marker belongs to element " << currentElem << std::endl;
    }

    if(nextElementFound == true) {
      std::cout << "The marker does not belong to element " << currentElem << std::endl;
    }


    std::cout << "markerIsInElement = " << markerIsInElement << " , " << "nextElementFound= " << nextElementFound << ", " << "nextElem = " << nextElem << std::endl;

    return (nextElem >= 0) ? nextElem : UINT_MAX;

  }

//   bool GetNewLocalCoordinates(std::vector <double> &xi, const std::vector< double > &x, const std::vector <double> &phi,
//                               const std::vector < std::vector <double > > &gradPhi, const std::vector < std::vector < std::vector <double> > > hessPhi,
//                               const std::vector < std::vector <double > > &a, const unsigned &dim, const unsigned &nDofs) {
// 
//     bool convergence = false;
//     std::vector < double > xp(dim, 0.);
//     std::vector < std::vector < double > > gradXp(dim);
//     std::vector < std::vector < std::vector < double > > > hessXp(dim);
// 
//     for(int k = 0; k < dim; k++) {
//       gradXp[k].assign(dim, 0.);
//       hessXp[k].resize(dim);
// 
//       for(int i1 = 0; i1 < dim; i1++) {
//         hessXp[k][i1].assign(dim, 0.);
//       }
//     }
// 
//     for(int k = 0; k < dim; k++) {
//       for(int i = 0; i < nDofs; i++) {
//         xp[k] += a[k][i] * phi[i];
// 
//         for(int i1 = 0; i1 < dim; i1++) {
//           gradXp[k][i1] += a[k][i] * gradPhi[i][i1];
// 
//           for(int i2 = 0; i2 < dim; i2++) {
//             hessXp[k][i1][i2] += a[k][i] * hessPhi[i][i1][i2];
//           }
//         }
//       }
//     }
// 
//     std::vector < double > gradF(dim, 0.);
//     std::vector < std::vector < double > >  hessF(dim);
// 
//     for(int i1 = 0; i1 < dim; i1++) {
//       hessF[i1].assign(dim, 0.);
//     }
// 
//     for(int k = 0; k < dim; k++) {
//       for(int i1 = 0; i1 < dim; i1++) {
//         gradF[i1] += -2. * (x[k] - xp[k]) * gradXp[k][i1];
// 
//         for(int i2 = 0; i2 < dim; i2++) {
//           hessF[i1][i2] += -2. * (x[k] - xp[k]) * hessXp[k][i1][i2] + 2. * gradXp[k][i1] * gradXp[k][i2];
//         }
//       }
//     }
// 
//     std::vector < std::vector < double > >  hessFm1(dim);
// 
//     for(int i1 = 0; i1 < dim; i1++) {
//       hessFm1[i1].resize(dim);
//     }
// 
//     if(dim == 2) {
//       double det = hessF[0][0] * hessF[1][1] - hessF[0][1] * hessF[1][0];
//       hessFm1[0][0] = hessF[1][1] / det;
//       hessFm1[0][1] = -hessF[0][1] / det;
//       hessFm1[1][0] = -hessF[1][0] / det;
//       hessFm1[1][1] = hessF[0][0] / det;
//     }
//     else if(dim == 3) {
//       double det = (hessF[0][0] * hessF[1][1] * hessF[2][2] + hessF[0][1] * hessF[1][2] * hessF[2][0] + hessF[0][2] * hessF[1][0] * hessF[2][1])
//                    - (hessF[2][0] * hessF[1][1] * hessF[0][2] + hessF[2][1] * hessF[1][2] * hessF[0][0] + hessF[2][2] * hessF[1][0] * hessF[0][1]) ;
// 
//       hessFm1[0][0] = (hessF[1][1] * hessF[2][2] - hessF[2][1] * hessF[1][2]) / det ;
//       hessFm1[0][1] = (hessF[0][2] * hessF[2][1] - hessF[2][2] * hessF[0][1]) / det ;
//       hessFm1[0][2] = (hessF[0][1] * hessF[1][2] - hessF[1][1] * hessF[0][2]) / det ;
//       hessFm1[1][0] = (hessF[1][2] * hessF[2][0] - hessF[2][2] * hessF[1][0]) / det ;
//       hessFm1[1][1] = (hessF[0][0] * hessF[2][2] - hessF[2][0] * hessF[0][2]) / det ;
//       hessFm1[1][2] = (hessF[0][2] * hessF[1][0] - hessF[0][0] * hessF[1][2]) / det ;
//       hessFm1[2][0] = (hessF[1][0] * hessF[2][1] - hessF[2][0] * hessF[1][1]) / det ;
//       hessFm1[2][1] = (hessF[0][1] * hessF[2][0] - hessF[2][1] * hessF[0][0]) / det ;
//       hessFm1[2][2] = (hessF[0][0] * hessF[1][1] - hessF[1][0] * hessF[0][1]) / det ;
//     }
// 
//     double delta2 = 0.;
// 
//     for(int i1 = 0; i1 < dim; i1++) {
//       double deltak = 0.;
// 
//       for(int i2 = 0; i2 < dim; i2++) {
//         deltak -= hessFm1[i1][i2] * gradF[i2];
//       }
// 
//       xi[i1] += deltak;
//       delta2 += deltak * deltak;
//     }
// 
// //     for(int k = 0; k < dim; k++) {
// //       std::cout << "xT[" << k << "]= " << xp[k] <<  " ";
// //     }
// //     std::cout << std::endl;
// 
//     if(delta2 < 1.0e-6) {
//       convergence = true;
//     }
// 
//     return convergence;
//   }



  bool GetNewLocalCoordinates(std::vector <double> &xi, const std::vector< double > &x, const std::vector <double> &phi,
                              const std::vector < std::vector <double > > &gradPhi, const std::vector < std::vector < std::vector <double> > > hessPhi,
                              const std::vector < std::vector <double > > &a, const unsigned &dim, const unsigned &nDofs) {

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
	  if(i1 == 1){
	    std::cout << i << " " << "gradphi[" << i << "][" << i1 << "]=" << gradPhi[i][i1] <<std::endl;
	  }
        }
      }
      if(k == 1){
	std::cout << "F[" << k << "]= " << F[k] << " " << " J[" << k << "][0]= " << J[k][0] <<" " << "J[" << k << "][1] = " << J[k][1] <<" " << "J[" << k <<"][2] =" << J[k][2] << std::endl; 
      }
      
      F[k] -= x[k];
    }
    
    

    std::vector < std::vector < double > >  Jm1(dim);

    for(int i1 = 0; i1 < dim; i1++) {
      Jm1[i1].resize(dim);
    }

    if(dim == 2) {
      double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
      Jm1[0][0] =  J[1][1] / det;
      Jm1[0][1] = -J[0][1] / det;
      Jm1[1][0] = -J[1][0] / det;
      Jm1[1][1] =  J[0][0] / det;
    }
    else if(dim == 3) {
      double det = (J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] + J[0][2] * J[1][0] * J[2][1])
                   - (J[2][0] * J[1][1] * J[0][2] + J[2][1] * J[1][2] * J[0][0] + J[2][2] * J[1][0] * J[0][1]) ;

      Jm1[0][0] = (J[1][1] * J[2][2] - J[2][1] * J[1][2]) / det ;
      Jm1[0][1] = (J[0][2] * J[2][1] - J[2][2] * J[0][1]) / det ;
      Jm1[0][2] = (J[0][1] * J[1][2] - J[1][1] * J[0][2]) / det ;
      Jm1[1][0] = (J[1][2] * J[2][0] - J[2][2] * J[1][0]) / det ;
      Jm1[1][1] = (J[0][0] * J[2][2] - J[2][0] * J[0][2]) / det ;
      Jm1[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det ;
      Jm1[2][0] = (J[1][0] * J[2][1] - J[2][0] * J[1][1]) / det ;
      Jm1[2][1] = (J[0][1] * J[2][0] - J[2][1] * J[0][0]) / det ;
      Jm1[2][2] = (J[0][0] * J[1][1] - J[1][0] * J[0][1]) / det ;
    }

    double delta2 = 0.;

    for(int i1 = 0; i1 < dim; i1++) {
      double deltak = 0.;

      for(int i2 = 0; i2 < dim; i2++) {
        deltak -= Jm1[i1][i2] * F[i2];
      }

      xi[i1] += deltak;
      delta2 += deltak * deltak;
    }

//     for(int k = 0; k < dim; k++) {
//       std::cout << "xT[" << k << "]= " << xp[k] <<  " ";
//     }
//     std::cout << std::endl;

    if(delta2 < 1.0e-6) {
      convergence = true;
    }

    return convergence;
  }

  void Marker::InverseMapping(const unsigned &iel, const unsigned &solType,
                              const std::vector< double > &x, std::vector< double > &xi) {

    
    
    unsigned dim =  _mesh->GetDimension();
    unsigned nDofs = _mesh->GetElementDofNumber(iel, solType);
    short unsigned ielType = _mesh->GetElementType(iel);

    std::cout << solType <<" "<< nDofs <<std::endl;
    
    //BEGIN extraction nodal coordinate values
    std::vector< std::vector < double > > xv(dim);

    for(unsigned k = 0; k < dim; k++) {
      xv[k].resize(nDofs);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned iDof  = _mesh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
      }
    }
    //END extraction


    //BEGIN projection nodal to polynomial coefficients
    std::vector < std::vector < double > > a;
    ProjectNodalToPolynomialCoefficients(a, xv, ielType, solType);
    
    for(int k=0;k<dim;k++){
      for(int i=0;i<nDofs;i++){
	std::cout << a[k][i] <<std::endl;
      }
      std::cout << std::endl;
    }
    
    //exit(0);
    
    //END projection


    //BEGIN inverse mapping search
    //std::vector <double> xi(dim, 0.);

    bool convergence = false;
    for(unsigned k = 0; k < dim; k++) {
      std::cout << "xi[" << k << "] = " << xi[k] << " ";
    }
    std::cout << std::endl;

    while(!convergence) {

      std::vector < double > phi;
      std::vector < std::vector < double > > gradPhi;
      std::vector < std::vector < std::vector < double > > > hessPhi;

      GetPolynomialShapeFunctionGradientHessian(phi, gradPhi, hessPhi, xi, ielType, solType);

      convergence = GetNewLocalCoordinates(xi, x, phi, gradPhi, hessPhi, a, dim, nDofs);
      for(unsigned k = 0; k < dim; k++) {
        std::cout << "xi[" << k << "] = " << xi[k] << " ";
      }
      std::cout << std::endl;
    }
    //END inverse mapping search

    return;


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


    // initial guess
    double xi = 1. / 3.;  //baricenter coordinates
    double eta = 1. / 3.;

    bool convergence = false;
    std::cout << "xi = " << xi << " , " << "eta = " << eta << std::endl;

    while(!convergence) {

      std::vector < double > phi(nDofs);
      std::vector < std::vector < double > > gradPhi(nDofs);
      std::vector < std::vector < std::vector < double > > > hessPhi(nDofs);

      for(int i = 0; i < nDofs; i++) {
        gradPhi[i].resize(dim);
        hessPhi[i].resize(dim);

        for(int i1 = 0; i1 < dim; i1++) {
          hessPhi[i][i1].resize(dim);
        }
      }

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

      //phi_x
      gradPhi[0][0] = 0.; // 0
      gradPhi[1][0] = 1.; // 1
      gradPhi[2][0] = 0.; // 0

      if(solutionType > 0) {
        gradPhi[3][0] = eta; // y
        gradPhi[4][0] = 2.*xi ;  // 2 x
        gradPhi[5][0] = 0.; // 0
      }

      if(solutionType > 1) {
        gradPhi[6][0] = gradPhi[4][0] * eta + phi[5]; // 2 x y + y y
      }

      //phi_y
      gradPhi[0][1] = 0.; // 0
      gradPhi[1][1] = 0.;  // 0
      gradPhi[2][1] = 1.;  // 1

      if(solutionType > 0) {
        gradPhi[3][1] = xi; // x
        gradPhi[4][1] = 0.;  // 0
        gradPhi[5][1] = 2.*eta; // 2*y
      }

      if(solutionType > 1) {
        gradPhi[6][1] = phi[4] + gradPhi[5][1] * xi; // xx  + 2 x y
      }

      //phi_xx
      hessPhi[0][0][0] = 0.; // 0
      hessPhi[1][0][0] = 0.; // 0
      hessPhi[2][0][0] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[3][0][0] = 0;  // 0
        hessPhi[4][0][0] = 2.; // 2
        hessPhi[5][0][0] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[6][0][0] = 2. * eta; // 2 y
      }

      //phi_xy
      hessPhi[0][1][0] = hessPhi[0][0][1] = 0.; // 0
      hessPhi[1][1][0] = hessPhi[1][0][1] = 0.; // 0
      hessPhi[2][1][0] = hessPhi[2][0][1] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[3][1][0] = hessPhi[3][0][1] = 1.; // 1
        hessPhi[4][1][0] = hessPhi[4][0][1] = 0.; // 0
        hessPhi[5][1][0] = hessPhi[5][0][1] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[6][1][0] = hessPhi[6][0][1] = 2. * xi + 2. * eta; // 2. x + 2 * y
      }

      //phi_yy
      hessPhi[0][1][1] = 0.; // 0
      hessPhi[1][1][1] = 0.; // 0
      hessPhi[2][1][1] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[3][1][1] = 0.; // 0
        hessPhi[4][1][1] = 0.; // 0
        hessPhi[5][1][1] = 2.; // 2
      }

      if(solutionType > 1) {
        hessPhi[6][1][1] = 2. * xi; // 2 * x
      }


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

      std::vector < std::vector < double > >  hessFm1(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessFm1[i1].resize(dim);
      }

      double det = hessF[0][0] * hessF[1][1] - hessF[0][1] * hessF[1][0];

      hessFm1[0][0] = hessF[1][1] / det;
      hessFm1[0][1] = -hessF[0][1] / det;
      hessFm1[1][0] = -hessF[1][0] / det;
      hessFm1[1][1] = hessF[0][0] / det;

      double dxi = hessFm1[0][0] * gradF[0] + hessFm1[0][1] * gradF[1];
      double deta = hessFm1[1][0] * gradF[0] + hessFm1[1][1] * gradF[1];

      xi  -= dxi;
      eta -= deta;

      if(dxi * dxi + deta * deta < 1.0e-6) {
        convergence = true;

        for(int k = 0; k < dim; k++) {
          std::cout << xp[k] << " ";
        }

        std::cout << std::endl;
      }

      //std::cout << "xi = " << xi <<" , "<< "eta = " << eta << std::endl;
    }

    std::vector <double> xcord(2, 0);
    xcord[0] = xi;
    xcord[1] = eta;
    return xcord;
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
    else if(solutionType == 1) {

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
                           2 * xv[k][15] + xv[k][1] + xv[k][2] + xv[k][3] - xv[k][4] - xv[k][5] - xv[k][6] - xv[k][7] - 2 * xv[k][8]);
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

    //initial guess
    double xi = 0.;
    double eta = 0.;
    double zita = 0.;

    bool convergence = false;
    std::cout << "xi = " << xi << " , " << "eta = " << eta << " , " << "zita = " << zita << std::endl;

    while(!convergence) {

      std::vector < double > phi(nDofs);
      std::vector < std::vector < double > > gradPhi(nDofs);
      std::vector < std::vector < std::vector < double > > > hessPhi(nDofs);

      for(int i = 0; i < nDofs; i++) {
        gradPhi[i].resize(dim);
        hessPhi[i].resize(dim);

        for(int i1 = 0; i1 < dim; i1++) {
          hessPhi[i][i1].resize(dim);
        }
      }

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


      //phi_x
      gradPhi[0][0] = 0.; // 0
      gradPhi[1][0] = 1.; // 1
      gradPhi[2][0] = 0.; // 0
      gradPhi[3][0] = 0.; // 0
      gradPhi[4][0] = eta ; // y
      gradPhi[5][0] = zita ; // z
      gradPhi[6][0] = 0.; // 0
      gradPhi[7][0] = phi[6];  //  y z

      if(solutionType > 0) {
        gradPhi[7][0] = 2 * xi;   // 2 x
        gradPhi[8][0] = 0.;   // 0
        gradPhi[9][0] = 0.;   // 0
        gradPhi[10][0] = phi[6];  // y z
        gradPhi[11][0] = 2 * phi[4];  // 2 x y
        gradPhi[12][0] = 2 * phi[5];  // 2 x z
        gradPhi[13][0] = 0.;  // 0
        gradPhi[14][0] = phi[8];  // yy
        gradPhi[15][0] = phi[9];  // zz
        gradPhi[16][0] = 0.;  // 0
        gradPhi[17][0] = 2 * phi[10]; // 2 x y z
        gradPhi[18][0] = phi[13]; //  yy z
        gradPhi[19][0] = phi[16]; // y zz
      }

      if(solutionType > 1) {
        gradPhi[20][0] = 2 * phi[14];  // 2 x yy
        gradPhi[21][0] = 2 * phi[15];  // 2 x zz
        gradPhi[22][0] = 0.;  // 0
        gradPhi[23][0] = 2 * phi[19]; // 2 x y zz
        gradPhi[24][0] = 2 * phi[18]; // 2 x yy z
        gradPhi[25][0] = phi[22]; //  yy zz
        gradPhi[26][0] = 2 * phi[25]; // 2 x yy zz
      }


      //phi_y
      gradPhi[0][1] = 0.; // 0
      gradPhi[1][1] = 0.;  // 0
      gradPhi[2][1] = 1.; // 1
      gradPhi[3][1] = 0.; // 0
      gradPhi[4][1] = xi; // x
      gradPhi[5][1] = 0.; // 0
      gradPhi[6][1] = zita; // z
      gradPhi[7][1] = phi[5];  // x z

      if(solutionType > 0) {
        gradPhi[7][1] = 0.;   // 0
        gradPhi[8][1] = 2 * eta;   // 2 y
        gradPhi[9][1] = 0.;   // 0
        gradPhi[10][1] = phi[5];  // x z
        gradPhi[11][1] = phi[7];  // xx
        gradPhi[12][1] = 0.;  // 0
        gradPhi[13][1] = 2 * phi[6];  // 2 y z
        gradPhi[14][1] = 2 * phi[4];  // 2 x y
        gradPhi[15][1] = 0.;  // 0
        gradPhi[16][1] = phi[9];  // zz
        gradPhi[17][1] = phi[12]; // xx  z
        gradPhi[18][1] = 2 * phi[10]; // 2 x y z
        gradPhi[19][1] = phi[15]; // x zz
      }

      if(solutionType > 1) {
        gradPhi[20][1] = 2 * phi[11];  // 2 xx y
        gradPhi[21][1] = 0.;  // 0
        gradPhi[22][1] = 2 * phi[16];  // 2 y zz
        gradPhi[23][1] = phi[21]; // xx zz
        gradPhi[24][1] = 2 * phi[17]; // 2 xx y z
        gradPhi[25][1] = 2 * phi[19]; // 2 x y zz
        gradPhi[26][1] = 2 * phi[23]; // 2 xx y zz
      }


      //phi_z
      gradPhi[0][2] = 0.; // 0
      gradPhi[1][2] = 0.; // 0
      gradPhi[2][2] = 0.; // 0
      gradPhi[3][2] = 1.; // 1
      gradPhi[4][2] = 0.; // 0
      gradPhi[5][2] = xi; // x
      gradPhi[6][2] = eta; // y
      gradPhi[7][2] = phi[4];  // x y

      if(solutionType > 0) {
        gradPhi[7][2] = 0.;   // 0
        gradPhi[8][2] = 0.;   // 0
        gradPhi[9][2] = 2 * zita;   // 2 z
        gradPhi[10][2] = phi[4];  // x y
        gradPhi[11][2] = 0.;  // 0
        gradPhi[12][2] = phi[7];  // xx
        gradPhi[13][2] = phi[8];  // yy
        gradPhi[14][2] = 0.;  // 0
        gradPhi[15][2] = 2 * phi[5];  // 2 x z
        gradPhi[16][2] = 2 * phi[6];  // 2 y z
        gradPhi[17][2] = phi[11]; // xx y
        gradPhi[18][2] = phi[14]; // x yy
        gradPhi[19][2] = 2 * phi[10]; // 2 x y z
      }

      if(solutionType > 1) {
        gradPhi[20][2] = 0.;  // 0
        gradPhi[21][2] = 2 * phi[12];  // 2 xx z
        gradPhi[22][2] = 2 * phi[13];  // 2 yy z
        gradPhi[23][2] = 2 * phi[17]; // 2 xx y z
        gradPhi[24][2] = phi[20]; // xx yy
        gradPhi[25][2] = 2 * phi[18]; // 2 x yy z
        gradPhi[26][2] = 2 * phi[24]; // 2 xx yy z
      }


      //phi_xx
      hessPhi[0][0][0] = 0.; // 0
      hessPhi[1][0][0] = 0.; // 0
      hessPhi[2][0][0] = 0.; // 0
      hessPhi[3][0][0] = 0.; // 0
      hessPhi[4][0][0] = 0.; // 0
      hessPhi[5][0][0] = 0.; // 0
      hessPhi[6][0][0] = 0.; // 0
      hessPhi[7][0][0] = 0.;  //  0

      if(solutionType > 0) {
        hessPhi[7][0][0] = 2.;   // 2
        hessPhi[8][0][0] = 0.;   // 0
        hessPhi[9][0][0] = 0.;   // 0
        hessPhi[10][0][0] = 0.;  // 0
        hessPhi[11][0][0] = 2 * eta;  // 2 y
        hessPhi[12][0][0] = 2 * zita;  // 2 z
        hessPhi[13][0][0] = 0.;  // 0
        hessPhi[14][0][0] = 0.;  // 0
        hessPhi[15][0][0] = 0.;  // 0
        hessPhi[16][0][0] = 0.;  // 0
        hessPhi[17][0][0] = 2 * phi[6]; // 2 y z
        hessPhi[18][0][0] = 0.; //  0
        hessPhi[19][0][0] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[20][0][0] = 2 * phi[8];  // 2 yy
        hessPhi[21][0][0] = 2 * phi[9];  // 2 zz
        hessPhi[22][0][0] = 0.;  // 0
        hessPhi[23][0][0] = 2 * phi[16]; // 2  y zz
        hessPhi[24][0][0] = 2 * phi[13]; // 2  yy z
        hessPhi[25][0][0] = 0.; //  0
        hessPhi[26][0][0] = 2 * phi[22]; // 2  yy zz
      }


      //phi_xy
      hessPhi[0][1][0] = hessPhi[0][0][1] = 0.; // 0
      hessPhi[1][1][0] = hessPhi[1][0][1] = 0.; // 0
      hessPhi[2][1][0] = hessPhi[2][0][1] = 0.; // 0
      hessPhi[3][1][0] = hessPhi[3][0][1] = 0.; // 0
      hessPhi[4][1][0] = hessPhi[4][0][1] = 1.; // 1
      hessPhi[5][1][0] = hessPhi[5][0][1] = 0.; // 0
      hessPhi[6][1][0] = hessPhi[6][0][1] = 0.; // 0
      hessPhi[7][1][0] = hessPhi[7][0][1] = zita;  // z

      if(solutionType > 0) {
        hessPhi[7][1][0] = hessPhi[7][0][1] = 0.;   // 0
        hessPhi[8][1][0] = hessPhi[8][0][1] = 0.;   // 0
        hessPhi[9][1][0] = hessPhi[9][0][1] = 0.;   // 0
        hessPhi[10][1][0] = hessPhi[10][0][1] = zita;  // z
        hessPhi[11][1][0] = hessPhi[11][0][1] = 2 * xi;  // 2 x
        hessPhi[12][1][0] = hessPhi[12][0][1] = 0.;  // 0.
        hessPhi[13][1][0] = hessPhi[13][0][1] = 0.;  // 0
        hessPhi[14][1][0] = hessPhi[14][0][1] = 2 * eta;  // 2 y
        hessPhi[15][1][0] = hessPhi[15][0][1] = 0.;  // 0
        hessPhi[16][1][0] = hessPhi[16][0][1] = 0.;  // 0
        hessPhi[17][1][0] = hessPhi[17][0][1] = 2 * phi[5]; // 2 x z
        hessPhi[18][1][0] = hessPhi[18][0][1] = 2 * phi[6]; //  2 y z
        hessPhi[19][1][0] = hessPhi[19][0][1] = phi[9]; // zz
      }

      if(solutionType > 1) {
        hessPhi[20][1][0] = hessPhi[20][0][1] = 4 * phi[4];  // 4 x y
        hessPhi[21][1][0] = hessPhi[21][0][1] = 0.;  // 0
        hessPhi[22][1][0] = hessPhi[22][0][1] = 0.;  // 0
        hessPhi[23][1][0] = hessPhi[23][0][1] = 2 * phi[15]; // 2 x zz
        hessPhi[24][1][0] = hessPhi[24][0][1] = 4 * phi[10]; // 4 x y z
        hessPhi[25][1][0] = hessPhi[25][0][1] = 2 * phi[16]; //  2 y zz
        hessPhi[26][1][0] = hessPhi[26][0][1] = 4 * phi[19]; // 4 x y zz
      }


      //phi_xz
      hessPhi[0][2][0] = hessPhi[0][0][2] = 0.; // 0
      hessPhi[1][2][0] = hessPhi[1][0][2] = 0.; // 0
      hessPhi[2][2][0] = hessPhi[2][0][2] = 0.; // 0
      hessPhi[3][2][0] = hessPhi[3][0][2] = 0.; // 0
      hessPhi[4][2][0] = hessPhi[4][0][2] = 0.; // 0
      hessPhi[5][2][0] = hessPhi[5][0][2] = 1.; // 1
      hessPhi[6][2][0] = hessPhi[6][0][2] = 0.; // 0
      hessPhi[7][2][0] = hessPhi[7][0][2] = eta;  //  y

      if(solutionType > 0) {
        hessPhi[7][2][0] = hessPhi[7][0][2] = 0.;   // 0
        hessPhi[8][2][0] = hessPhi[8][0][2] = 0.;   // 0
        hessPhi[9][2][0] = hessPhi[9][0][2] = 0.;   // 0
        hessPhi[10][2][0] = hessPhi[10][0][2] = eta;  // y
        hessPhi[11][2][0] = hessPhi[11][0][2] = 0.;  // 0
        hessPhi[12][2][0] = hessPhi[12][0][2] = 2 * xi;  // 2 x
        hessPhi[13][2][0] = hessPhi[13][0][2] = 0.;  // 0
        hessPhi[14][2][0] = hessPhi[14][0][2] = 0.;  // 0
        hessPhi[15][2][0] = hessPhi[15][0][2] = 2 * zita;  // 2 z
        hessPhi[16][2][0] = hessPhi[16][0][2] = 0.;  // 0
        hessPhi[17][2][0] = hessPhi[17][0][2] = 2 * phi[4]; // 2 x y
        hessPhi[18][2][0] = hessPhi[18][0][2] = phi[8]; //  yy
        hessPhi[19][2][0] = hessPhi[19][0][2] = 2 * phi[6]; // 2 y z
      }

      if(solutionType > 1) {
        hessPhi[20][2][0] = hessPhi[20][0][2] = 0.;  // 0
        hessPhi[21][2][0] = hessPhi[21][0][2] = 4 * phi[5];  // 4 x z
        hessPhi[22][2][0] = hessPhi[22][0][2] = 0.;  // 0
        hessPhi[23][2][0] = hessPhi[23][0][2] = 4 * phi[10]; // 4 x y z
        hessPhi[24][2][0] = hessPhi[24][0][2] = 2 * phi[14]; // 2 x yy
        hessPhi[25][2][0] = hessPhi[25][0][2] = 2 * phi[13]; //  2 yy z
        hessPhi[26][2][0] = hessPhi[26][0][2] = 4 * phi[18]; // 4 x yy z
      }


      //phi_yy
      hessPhi[0][1][1] = 0.; // 0
      hessPhi[1][1][1] = 0.;  // 0
      hessPhi[2][1][1] = 0.; // 0
      hessPhi[3][1][1] = 0.; // 0
      hessPhi[4][1][1] = 0.; // 0
      hessPhi[5][1][1] = 0.; // 0
      hessPhi[6][1][1] = 0.; // 0
      hessPhi[7][1][1] = 0.;  // 0

      if(solutionType > 0) {
        hessPhi[7][1][1] = 0.;   // 0
        hessPhi[8][1][1] = 2.;   // 2
        hessPhi[9][1][1] = 0.;   // 0
        hessPhi[10][1][1] = 0.;  // 0
        hessPhi[11][1][1] = 0.;  // 0
        hessPhi[12][1][1] = 0.;  // 0
        hessPhi[13][1][1] = 2 * zita;  // 2 z
        hessPhi[14][1][1] = 2 * xi;  // 2 x
        hessPhi[15][1][1] = 0.;  // 0
        hessPhi[16][1][1] = 0.;  // 0
        hessPhi[17][1][1] = 0.; // 0
        hessPhi[18][1][1] = 2 * phi[5]; // 2 x z
        hessPhi[19][1][1] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[20][1][1] = 2 * phi[7];  // 2 xx
        hessPhi[21][1][1] = 0.;  // 0
        hessPhi[22][1][1] = 2 * phi[9];  // 2 zz
        hessPhi[23][1][1] = 0.; // 0
        hessPhi[24][1][1] = 2 * phi[12]; // 2 xx z
        hessPhi[25][1][1] = 2 * phi[15]; // 2 x zz
        hessPhi[26][1][1] = 2 * phi[21]; // 2 xx zz
      }


      //phi_yz
      hessPhi[0][2][1] = hessPhi[0][1][2] = 0.; // 0
      hessPhi[1][2][1] = hessPhi[1][1][2] = 0.;  // 0
      hessPhi[2][2][1] = hessPhi[2][1][2] = 0.; // 0
      hessPhi[3][2][1] = hessPhi[3][1][2] = 0.; // 0
      hessPhi[4][2][1] = hessPhi[4][1][2] = 0.; // 0
      hessPhi[5][2][1] = hessPhi[5][1][2] = 0.; // 0
      hessPhi[6][2][1] = hessPhi[6][1][2] = 1.; // 1
      hessPhi[7][2][1] = hessPhi[7][1][2] = xi ;  // x

      if(solutionType > 0) {
        hessPhi[7][2][1] = hessPhi[7][1][2] = 0.;   // 0
        hessPhi[8][2][1] = hessPhi[8][1][2] = 0.;   // 0
        hessPhi[9][2][1] = hessPhi[9][1][2] = 0.;   // 0
        hessPhi[10][2][1] = hessPhi[10][1][2] = xi ;  // x
        hessPhi[11][2][1] = hessPhi[11][1][2] = 0.;  // 0
        hessPhi[12][2][1] = hessPhi[12][1][2] = 0.;  // 0
        hessPhi[13][2][1] = hessPhi[13][1][2] = 2 * eta ;  // 2 y
        hessPhi[14][2][1] = hessPhi[14][1][2] = 0.;  // 0
        hessPhi[15][2][1] = hessPhi[15][1][2] = 0.;  // 0
        hessPhi[16][2][1] = hessPhi[16][1][2] = 2 * zita ;  // 2 z
        hessPhi[17][2][1] = hessPhi[17][1][2] = phi[7]; // xx
        hessPhi[18][2][1] = hessPhi[18][1][2] = 2 * phi[4]; // 2 x y
        hessPhi[19][2][1] = hessPhi[19][1][2] = 2 * phi[5]; // 2 x z
      }

      if(solutionType > 1) {
        hessPhi[20][2][1] = hessPhi[20][1][2] = 0.;  // 0
        hessPhi[21][2][1] = hessPhi[21][1][2] = 0.;  // 0
        hessPhi[22][2][1] = hessPhi[22][1][2] = 4 * phi[6];  // 4 y z
        hessPhi[23][2][1] = hessPhi[23][1][2] = 2 * phi[12]; // 2 xx z
        hessPhi[24][2][1] = hessPhi[24][1][2] = 2 * phi[11]; // 2 xx y
        hessPhi[25][2][1] = hessPhi[25][1][2] = 4 * phi[10]; // 4 x y z
        hessPhi[26][2][1] = hessPhi[26][1][2] = 4 * phi[17]; // 4 xx y z
      }


      //phi_zz
      hessPhi[0][2][2] = 0.; // 0
      hessPhi[1][2][2] = 0.; // 0
      hessPhi[2][2][2] = 0.; // 0
      hessPhi[3][2][2] = 0.; // 0
      hessPhi[4][2][2] = 0.; // 0
      hessPhi[5][2][2] = 0.; // 0
      hessPhi[6][2][2] = 0.; // 0
      hessPhi[7][2][2] = 0.;  // 0

      if(solutionType > 0) {
        hessPhi[7][2][2] = 0.;   // 0
        hessPhi[8][2][2] = 0.;   // 0
        hessPhi[9][2][2] = 2.;   // 2
        hessPhi[10][2][2] = 0.;  // 0
        hessPhi[11][2][2] = 0.;  // 0
        hessPhi[12][2][2] = 0.;  // 0
        hessPhi[13][2][2] = 0.;  // 0
        hessPhi[14][2][2] = 0.;  // 0
        hessPhi[15][2][2] = 2 * xi;  // 2 x
        hessPhi[16][2][2] = 2 * eta;  // 2 y
        hessPhi[17][2][2] = 0.; // 0
        hessPhi[18][2][2] = 0.; // 0
        hessPhi[19][2][2] = 2 * phi[4]; // 2 x y
      }

      if(solutionType > 1) {
        hessPhi[20][2][2] = 0.;  // 0
        hessPhi[21][2][2] = 2 * phi[7];  // 2 xx
        hessPhi[22][2][2] = 2 * phi[8];  // 2 yy
        hessPhi[23][2][2] = 2 * phi[11]; // 2 xx y
        hessPhi[24][2][2] = 0.; // 0
        hessPhi[25][2][2] = 2 * phi[14]; // 2 x yy
        hessPhi[26][2][2] = 2 * phi[20]; // 2 xx yy
      }


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

      std::vector < std::vector < double > >  hessFm1(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessFm1[i1].resize(dim);
      }

      double det = (hessF[0][0] * hessF[1][1] * hessF[2][2] + hessF[0][1] * hessF[1][2] * hessF[2][0] + hessF[0][2] * hessF[1][0] * hessF[2][1])
                   - (hessF[2][0] * hessF[1][1] * hessF[0][2] + hessF[2][1] * hessF[1][2] * hessF[0][0] + hessF[2][2] * hessF[1][0] * hessF[0][1]) ;

      hessFm1[0][0] = (hessF[1][1] * hessF[2][2] - hessF[2][1] * hessF[1][2]) / det ;
      hessFm1[0][1] = (hessF[0][2] * hessF[2][1] - hessF[2][2] * hessF[0][1]) / det ;
      hessFm1[0][2] = (hessF[0][1] * hessF[1][2] - hessF[1][1] * hessF[0][2]) / det ;
      hessFm1[1][0] = (hessF[1][2] * hessF[2][0] - hessF[2][2] * hessF[1][0]) / det ;
      hessFm1[1][1] = (hessF[0][0] * hessF[2][2] - hessF[2][0] * hessF[0][2]) / det ;
      hessFm1[1][2] = (hessF[0][2] * hessF[1][0] - hessF[0][0] * hessF[1][2]) / det ;
      hessFm1[2][0] = (hessF[1][0] * hessF[2][1] - hessF[2][0] * hessF[1][1]) / det ;
      hessFm1[2][1] = (hessF[0][1] * hessF[2][0] - hessF[2][1] * hessF[0][0]) / det ;
      hessFm1[2][2] = (hessF[0][0] * hessF[1][1] - hessF[1][0] * hessF[0][1]) / det ;

      double dxi = hessFm1[0][0] * gradF[0] + hessFm1[0][1] * gradF[1] + hessFm1[0][2] * gradF[2];
      double deta = hessFm1[1][0] * gradF[0] + hessFm1[1][1] * gradF[1] + hessFm1[1][2] * gradF[2];
      double dzita = hessFm1[2][0] * gradF[0] + hessFm1[2][1] * gradF[1] + hessFm1[2][2] * gradF[2];

      xi  -= dxi;
      eta -= deta;
      zita -= dzita;

      if(dxi * dxi + deta * deta + dzita * dzita < 1.0e-6) {
        convergence = true;

        for(int k = 0; k < dim; k++) {
          std::cout << xp[k] << " ";
        }

        std::cout << std::endl;
      }

      //std::cout << "xi = " << xi <<" , "<< "eta = " << eta << " , " << "zita= " << zita <<  std::endl;
    }

    std::vector <double> xcord(3, 0);
    xcord[0] = xi;
    xcord[1] = eta;
    xcord[2] = zita;

    return xcord;

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

        a[k][0] = xv[k][0] ;
        a[k][1] = - xv[k][0] + xv[k][1] ;
        a[k][2] = - xv[k][0] + xv[k][2] ;
        a[k][3] = - xv[k][0] + xv[k][3] ;
      }
    }
    else if(solutionType == 1) {

      for(int k = 0; k < dim; k++) {

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


    //initial guess
    double xi = 1 / 4 ;
    double eta = 1 / 4 ;
    double zita = 1 / 4 ;

    bool convergence = false;
    std::cout << "xi = " << xi << " , " << "eta = " << eta << " , " << "zita = " << zita << std::endl;

    while(!convergence) {


      std::vector < double > phi(nDofs);
      std::vector < std::vector < double > > gradPhi(nDofs);
      std::vector < std::vector < std::vector < double > > > hessPhi(nDofs);

      for(int i = 0; i < nDofs; i++) {
        gradPhi[i].resize(dim);
        hessPhi[i].resize(dim);

        for(int i1 = 0; i1 < dim; i1++) {
          hessPhi[i][i1].resize(dim);
        }
      }

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
        phi[12] = phi[1] * phi[9] + phi[7] * phi[3] ; // x z z + x x z
        phi[13] = phi[2] * phi[9] + phi[8] * phi[3] ; // y z z + y y z
        phi[14] = phi[4] * phi[9] + phi[4] * phi[6] + phi[7] * phi[6]; // x y z z + x y y z + x x y z
      }


      //phi_x
      gradPhi[0][0] = 0.; // 0
      gradPhi[1][0] = 1.; // 1
      gradPhi[2][0] = 0.;  // 0
      gradPhi[3][0] = 0.; // 0

      if(solutionType > 0) {
        gradPhi[4][0] = eta ;  //  y
        gradPhi[5][0] = zita ;  //  z
        gradPhi[6][0] = 0.;  // 0
        gradPhi[7][0] = 2 * xi ; // 2 x
        gradPhi[8][0] = 0.; // 0
        gradPhi[9][0] = 0.; // 0
      }

      if(solutionType > 1) {
        gradPhi[10][0] = phi[6]; //  y z
        gradPhi[11][0] = phi[8] + 2 * phi[4] ; //  y y + 2 x y
        gradPhi[12][0] = phi[9] + 2 * phi[5] ; // z z + 2 x z
        gradPhi[13][0] = 0.; // 0;
        gradPhi[14][0] = eta * (gradPhi[12][0] + phi[6]) ; //  y z z +  y y z + 2 x y z
      }


      //phi_y
      gradPhi[0][1] = 0.; // 0
      gradPhi[1][1] = 0.; // 0
      gradPhi[2][1] = 1.;  // 1
      gradPhi[3][1] = 0.; // 0

      if(solutionType > 0) {
        gradPhi[4][1] = xi ;  // x
        gradPhi[5][1] = 0.;  // 0
        gradPhi[6][1] = zita ;  //  z
        gradPhi[7][1] = 0.; // 0
        gradPhi[8][1] = 2 * eta; // 2 y
        gradPhi[9][1] = 0.; // 0
      }

      if(solutionType > 1) {
        gradPhi[10][1] = phi[5]; // x z
        gradPhi[11][1] = 2 * phi[4] + phi[7] ; //  2 x y + x x
        gradPhi[12][1] = 0.; // 0
        gradPhi[13][1] = phi[9] + 2 * phi[6]; //  z z + 2 y z
        gradPhi[14][1] = xi * (gradPhi[13][1] + phi[5]); // x z z + 2 x y z + x x z
      }


      //phi_z
      gradPhi[0][2] = 0.; // 0
      gradPhi[1][2] = 0.; // 0
      gradPhi[2][2] = 0.;  // 0
      gradPhi[3][2] = 1.; // 1

      if(solutionType > 0) {
        gradPhi[4][2] = 0.;  // 0
        gradPhi[5][2] = xi;  // x
        gradPhi[6][2] = eta;  // y
        gradPhi[7][2] = 0.; // 0
        gradPhi[8][2] = 0.; // 0
        gradPhi[9][2] = 2 * zita; // 2 z
      }

      if(solutionType > 1) {
        gradPhi[10][2] = phi[4]; // x y
        gradPhi[11][2] = 0.; // 0
        gradPhi[12][2] = 2 * phi[5] + phi[7]; // 2 x z  + x x
        gradPhi[13][2] = 2 * phi[6] + phi[8]; // 2 y z  + y y
        gradPhi[14][2] = xi * (gradPhi[13][2] + phi[4]); //  2 x y z  + x y y  + x x y
      }


      //phi_xx
      hessPhi[0][0][0] = 0.; // 0
      hessPhi[1][0][0] = 0.; // 0
      hessPhi[2][0][0] = 0.;  // 0
      hessPhi[3][0][0] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[4][0][0] = 0.;  //  0
        hessPhi[5][0][0] = 0.;  //  0
        hessPhi[6][0][0] = 0.;  // 0
        hessPhi[7][0][0] = 2.; // 2
        hessPhi[8][0][0] = 0.; // 0
        hessPhi[9][0][0] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[10][0][0] = 0.; //  0
        hessPhi[11][0][0] = 2 * eta; //  2 y
        hessPhi[12][0][0] = 2 * zita; // 2 z
        hessPhi[13][0][0] = 0.; // 0;
        hessPhi[14][0][0] = 2 * phi[6]; //  2 y z
      }


      //phi_xy
      hessPhi[0][1][0] = hessPhi[0][0][1] = 0.; // 0
      hessPhi[1][1][0] = hessPhi[1][0][1] = 0.; // 0
      hessPhi[2][1][0] = hessPhi[2][0][1] = 0.;  // 0
      hessPhi[3][1][0] = hessPhi[3][0][1] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[4][1][0] = hessPhi[4][0][1] = 1.;  //  1
        hessPhi[5][1][0] = hessPhi[5][0][1] = 0.;  //  0
        hessPhi[6][1][0] = hessPhi[6][0][1] = 0.;  // 0
        hessPhi[7][1][0] = hessPhi[7][0][1] = 0.; // 0
        hessPhi[8][1][0] = hessPhi[8][0][1] = 0.; // 0
        hessPhi[9][1][0] = hessPhi[9][0][1] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[10][1][0] = hessPhi[10][0][1] = zita ; //  z
        hessPhi[11][1][0] = hessPhi[11][0][1] = 2 * eta + 2 * xi ; //  2 y + 2 x
        hessPhi[12][1][0] = hessPhi[12][0][1] = 0.; // 0
        hessPhi[13][1][0] = hessPhi[13][0][1] = 0.; // 0;
        hessPhi[14][1][0] = hessPhi[14][0][1] = gradPhi[13][1] + 2 * phi[5] ; //  z z +  2 y z + 2 x z
      }


      //phi_xz
      hessPhi[0][2][0] = hessPhi[0][0][2] = 0.; // 0
      hessPhi[1][2][0] = hessPhi[1][0][2] = 0.; // 0
      hessPhi[2][2][0] = hessPhi[2][0][2] = 0.;  // 0
      hessPhi[3][2][0] = hessPhi[3][0][2] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[4][2][0] = hessPhi[4][0][2] = 0.;  //  0
        hessPhi[5][2][0] = hessPhi[5][0][2] = 1.;  //  1
        hessPhi[6][2][0] = hessPhi[6][0][2] = 0.;  // 0
        hessPhi[7][2][0] = hessPhi[7][0][2] = 0.; // 0
        hessPhi[8][2][0] = hessPhi[8][0][2] = 0.; // 0
        hessPhi[9][2][0] = hessPhi[9][0][2] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[10][2][0] = hessPhi[10][0][2] = eta ; //  y
        hessPhi[11][2][0] = hessPhi[11][0][2] = 0.; //  0
        hessPhi[12][2][0] = hessPhi[12][0][2] = 2 * (xi + zita) ; // 2 z + 2 x
        hessPhi[13][2][0] = hessPhi[13][0][2] = 0.; // 0;
        hessPhi[14][2][0] = hessPhi[14][0][2] = gradPhi[13][2] +  2 * phi[4] ; //  2 y z +  y y  + 2 x y
      }


      //phi_yy
      hessPhi[0][1][1] = 0.; // 0
      hessPhi[1][1][1] = 0.; // 0
      hessPhi[2][1][1] = 0.; // 0
      hessPhi[3][1][1] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[4][1][1] = 0.;  // 0
        hessPhi[5][1][1] = 0.;  // 0
        hessPhi[6][1][1] = 0.;  //  0
        hessPhi[7][1][1] = 0.; // 0
        hessPhi[8][1][1] = 2.; // 2
        hessPhi[9][1][1] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[10][1][1] = 0.; // 0
        hessPhi[11][1][1] = 2 * xi ; //  2 x
        hessPhi[12][1][1] = 0.; // 0
        hessPhi[13][1][1] = 2 * zita; //  2 z
        hessPhi[14][1][1] = 2 * phi[5]; // 2 x z
      }


      //phi_yz
      hessPhi[0][2][1] = hessPhi[0][1][2] = 0.; // 0
      hessPhi[1][2][1] = hessPhi[1][1][2] = 0.; // 0
      hessPhi[2][2][1] = hessPhi[2][1][2] = 0.;  // 0
      hessPhi[3][2][1] = hessPhi[3][1][2] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[4][2][1] = hessPhi[4][1][2] = 0.;  // 0
        hessPhi[5][2][1] = hessPhi[5][1][2] = 0.;  // 0
        hessPhi[6][2][1] = hessPhi[6][1][2] = 1.;  //  1
        hessPhi[7][2][1] = hessPhi[7][1][2] = 0.; // 0
        hessPhi[8][2][1] = hessPhi[8][1][2] = 0.; // 0
        hessPhi[9][2][1] = hessPhi[9][1][2] = 0.; // 0
      }

      if(solutionType > 1) {
        hessPhi[10][2][1] = hessPhi[10][1][2] = xi ; // x
        hessPhi[11][2][1] = hessPhi[11][1][2] = 0.; //  0
        hessPhi[12][2][1] = hessPhi[12][1][2] = 0.; // 0
        hessPhi[13][2][1] = hessPhi[13][1][2] = 2 * zita + 2 * eta; //  2 z + 2 y
        hessPhi[14][2][1] = hessPhi[14][1][2] = 2 * phi[5] + gradPhi[11][1] ; // 2 x z  + 2 x y  + x x
      }


      //phi_zz
      hessPhi[0][2][2] = 0.; // 0
      hessPhi[1][2][2] = 0.; // 0
      hessPhi[2][2][2] = 0.; // 0
      hessPhi[3][2][2] = 0.; // 0

      if(solutionType > 0) {
        hessPhi[4][2][2] = 0.;  // 0
        hessPhi[5][2][2] = 0.;  // 0
        hessPhi[6][2][2] = 0.;  // 0
        hessPhi[7][2][2] = 0.; // 0
        hessPhi[8][2][2] = 0.; // 0
        hessPhi[9][2][2] = 2.; // 2
      }

      if(solutionType > 1) {
        hessPhi[10][2][2] = 0.; // 0
        hessPhi[11][2][2] = 0.; // 0
        hessPhi[12][2][2] = 2 * xi ; // 2 x
        hessPhi[13][2][2] = 2 * eta ; // 2 y
        hessPhi[14][2][2] = 2 * phi[4]; //  2 x y
      }

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

      std::vector < std::vector < double > >  hessFm1(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessFm1[i1].resize(dim);
      }

      double det = (hessF[0][0] * hessF[1][1] * hessF[2][2] + hessF[0][1] * hessF[1][2] * hessF[2][0] + hessF[0][2] * hessF[1][0] * hessF[2][1])
                   - (hessF[2][0] * hessF[1][1] * hessF[0][2] + hessF[2][1] * hessF[1][2] * hessF[0][0] + hessF[2][2] * hessF[1][0] * hessF[0][1]) ;

      hessFm1[0][0] = (hessF[1][1] * hessF[2][2] - hessF[2][1] * hessF[1][2]) / det ;
      hessFm1[0][1] = (hessF[0][2] * hessF[2][1] - hessF[2][2] * hessF[0][1]) / det ;
      hessFm1[0][2] = (hessF[0][1] * hessF[1][2] - hessF[1][1] * hessF[0][2]) / det ;
      hessFm1[1][0] = (hessF[1][2] * hessF[2][0] - hessF[2][2] * hessF[1][0]) / det ;
      hessFm1[1][1] = (hessF[0][0] * hessF[2][2] - hessF[2][0] * hessF[0][2]) / det ;
      hessFm1[1][2] = (hessF[0][2] * hessF[1][0] - hessF[0][0] * hessF[1][2]) / det ;
      hessFm1[2][0] = (hessF[1][0] * hessF[2][1] - hessF[2][0] * hessF[1][1]) / det ;
      hessFm1[2][1] = (hessF[0][1] * hessF[2][0] - hessF[2][1] * hessF[0][0]) / det ;
      hessFm1[2][2] = (hessF[0][0] * hessF[1][1] - hessF[1][0] * hessF[0][1]) / det ;

      double dxi = hessFm1[0][0] * gradF[0] + hessFm1[0][1] * gradF[1] + hessFm1[0][2] * gradF[2];
      double deta = hessFm1[1][0] * gradF[0] + hessFm1[1][1] * gradF[1] + hessFm1[1][2] * gradF[2];
      double dzita = hessFm1[2][0] * gradF[0] + hessFm1[2][1] * gradF[1] + hessFm1[2][2] * gradF[2];

      xi  -= dxi;
      eta -= deta;
      zita -= dzita;

      if(dxi * dxi + deta * deta + dzita * dzita < 1.0e-6) {
        convergence = true;

        for(int k = 0; k < dim; k++) {
          std::cout << xp[k] << " ";
        }

        std::cout << std::endl;
      }

      //std::cout << "xi = " << xi <<" , "<< "eta = " << eta << " , " << "zita= " << zita <<  std::endl;
    }

    std::vector <double> xcord(3, 0);
    xcord[0] = xi;
    xcord[1] = eta;
    xcord[2] = zita;

    return xcord;

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

//             for(int i=0 ; i <nDofs;  i++) {
//                 xv[k][i] = i;
//             }

        a[k][0] = xv[k][12];
        a[k][1] = - xv[k][0] + 2 * xv[k][9] - xv[k][12] + xv[k][13] - xv[k][1] - xv[k][3] - xv[k][4] + 2 * xv[k][6];
        a[k][2] = - xv[k][0] + 2 * xv[k][11] - xv[k][12] + xv[k][14] - xv[k][2] - xv[k][3] - xv[k][5] + 2 * xv[k][8];
        a[k][3] = 0.5 * (- xv[k][0] + xv[k][3]);
        a[k][4] = 2 * (xv[k][0] - xv[k][9] + xv[k][10] - xv[k][11] + xv[k][3] - xv[k][6] + xv[k][7] - xv[k][8]);
        a[k][5] = 0.5 * (3 * xv[k][0] + 4 * xv[k][9] + xv[k][1] - 3 * xv[k][3] - xv[k][4] - 4 * xv[k][6]);
        a[k][6] = 0.5 * (3 * xv[k][0] + 4 * xv[k][11] + xv[k][2] - 3 * xv[k][3] - xv[k][5] - 4 * xv[k][8]);
        a[k][7] = xv[k][0] - 2 * xv[k][9] + xv[k][1] + xv[k][3] + xv[k][4] - 2 * xv[k][6];
        a[k][8] = xv[k][0] - 2 * xv[k][11] + xv[k][2] + xv[k][3] + xv[k][5] - 2 * xv[k][8];
        a[k][9] = 0.5 * (xv[k][0] - 2 * xv[k][12] + xv[k][3]);
        a[k][10] = -2 * (xv[k][0] + xv[k][9] - xv[k][10] + xv[k][11] - xv[k][3] - xv[k][6] + xv[k][7] - xv[k][8]);
        a[k][11] = - xv[k][0] - 2 * xv[k][9] - xv[k][1] + xv[k][3] + xv[k][4] + 2 * xv[k][6];
        a[k][12] = - xv[k][0] - 2 * xv[k][11] - xv[k][2] + xv[k][3] + xv[k][5] + 2 * xv[k][8];
        a[k][13] = 0.5 * (- xv[k][0] + 2 * xv[k][12] - 2 * xv[k][13] + xv[k][1] - xv[k][3] + xv[k][4]);
        a[k][14] = 0.5 * (- xv[k][0] + 2 * xv[k][12] - 2 * xv[k][14] + xv[k][2] - xv[k][3] + xv[k][5]);
//             if(k == 0) {
//                 for(int j = 14; j >= 0; j--) {
//                     std::cout << j << " " << a[k][j] << std::endl;
//                 }
//             }
      }
    }
    else if(solutionType == 2) {  // I think the coefficients get too small
      for(int k = 0; k < dim; k++) {

//             for(int i=0 ; i <nDofs;  i++){
//               xv[k][i] = i;
//             }

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
        a[k][14] = 0.5 * (-3 * xv[k][0] + 4 * xv[k][9] + 6 * xv[k][12] + 2 * xv[k][13] - 8 * xv[k][15] - xv[k][1] - 3 * xv[k][3] - xv[k][4] + 4 * xv[k][6]);
        a[k][15] = 0.5 * (-3 * xv[k][0] + 4 * xv[k][11] + 6 * xv[k][12] + 2 * xv[k][14] - 8 * xv[k][17] - xv[k][2] - 3 * xv[k][3] - xv[k][5] + 4 * xv[k][8]);
        a[k][16] = xv[k][0] - 2 * xv[k][9] - 2 * xv[k][12] - 2 * xv[k][13] + 4 * xv[k][15] + xv[k][1] + xv[k][3] + xv[k][4] - 2 * xv[k][6];
        a[k][17] = xv[k][0] - 2 * xv[k][11] - 2 * xv[k][12] - 2 * xv[k][14] + 4 * xv[k][17] + xv[k][2] + xv[k][3] + xv[k][5] - 2 * xv[k][8];
        a[k][18] = 1.5 * (xv[k][0] + 4 * xv[k][9] + 4 * xv[k][10] + 4 * xv[k][11] + 9 * xv[k][18] + xv[k][1] - 9 * xv[k][19] + xv[k][2] - xv[k][3] - xv[k][4] - xv[k][5] - 4 * (xv[k][6] + xv[k][7] + xv[k][8]));
        a[k][19] = 0.5 * (7 * xv[k][0] - 16 * xv[k][9] - 8 * xv[k][10] - 16 * xv[k][11] - 14 * xv[k][12] - 6 * xv[k][13] - 6 * xv[k][14] + 32 * xv[k][15] +
                          16 * xv[k][16] + 32 * xv[k][17] + 27 * xv[k][18] + 3 * xv[k][1] + 27 * xv[k][19] - 54 * xv[k][20] + 3 * xv[k][2] + 7 * xv[k][3] +
                          3 * xv[k][4] + 3 * xv[k][5] - 8 * (2 * xv[k][6] + xv[k][7] + 2 * xv[k][8]));
        a[k][20] = -1.5 * (xv[k][0] - 4 * xv[k][9] - 4 * xv[k][10] - 4 * xv[k][11] - 2 * xv[k][12] - 2 * xv[k][13] - 2 * xv[k][14] + 8 * xv[k][15] + 8 * xv[k][16] + 8 * xv[k][17] + 9 * xv[k][18] + xv[k][1] + 9 * xv[k][19] - 18 * xv[k][20] +
                           xv[k][2] + xv[k][3] + xv[k][4] + xv[k][5] - 4 * (xv[k][6] + xv[k][7] + xv[k][8]));
// 	            if(k == 0) {
//           for(int j = 20; j >= 0; j--) {
//             std::cout << j << " " << a[k][j] << std::endl;
//           }
//         }

// 	    	            if(k == 0) {
//           for(int j = 0; j < 20; j++) {
//             std::cout << "j = " << j << " a[0][" << j << "] =" << a[k][j] << std::endl;
//           }
//         }


      }
    }


    //initial guess
    double xi = 1 / 3;
    double eta = 1 / 3;
    double zita = 0.;

    bool convergence = false;
    std::cout << "xi = " << xi << " , " << "eta = " << eta << " , " << "zita = " << zita << std::endl;

    while(!convergence) {

      std::vector < double > phi(nDofs);
      std::vector < std::vector < double > > gradPhi(nDofs);
      std::vector < std::vector < std::vector < double > > > hessPhi(nDofs);

      for(int i = 0; i < nDofs; i++) {
        gradPhi[i].resize(dim);
        hessPhi[i].resize(dim);

        for(int i1 = 0; i1 < dim; i1++) {
          hessPhi[i][i1].resize(dim);
        }
      }


      phi[0] = 1.;
      phi[1] = xi; // x
      phi[2] = eta;  // y
      phi[3] = zita; // z
      phi[4] = phi[1] * phi[3];  // x z
      phi[5] = phi[2] * phi[3];  // y z

      if(solutionType > 0) {
        phi[4] = xi * eta;  // x y
        phi[5] = xi * zita;  // x z
        phi[6] = eta * zita;  // y z
        phi[7] = xi * xi; // x x
        phi[8] = eta * eta; // y y
        phi[9] = zita * zita; // z z
        phi[10] = phi[4] * zita; // x y z
        phi[11] = phi[7] * zita; // x x z
        phi[12] = phi[8] * zita; // y y z
        phi[13] = xi * phi[9]; // x z z
        phi[14] = eta * phi[9]; // y z z;
      }

      if(solutionType > 1) {

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


      //phi_x
      gradPhi[0][0] = 0.; // 0
      gradPhi[1][0] = 1.; // 1
      gradPhi[2][0] = 0.;  // 0
      gradPhi[3][0] = 0.; // 0
      gradPhi[4][0] = zita ;  //  z
      gradPhi[5][0] = 0.;  // 0

      if(solutionType > 0) {
        gradPhi[4][0] = eta ;  //  y
        gradPhi[5][0] = zita ;  //  z
        gradPhi[6][0] = 0.;  // 0
        gradPhi[7][0] = 2 * xi ; // 2 x
        gradPhi[8][0] = 0.; // 0
        gradPhi[9][0] = 0.; // 0
        gradPhi[10][0] = phi[6]; //  y z
        gradPhi[11][0] = 2 * phi[5]; // 2 x z
        gradPhi[12][0] = 0.; // 0
        gradPhi[13][0] = phi[9]; //  z z
        gradPhi[14][0] = 0.; // 0 ;
      }

      if(solutionType > 1) {

        gradPhi[11][0] = 2 * phi[4] + phi[8]; // 2 x y +  yy
        gradPhi[12][0] = 2 * phi[5]; // 2 x z
        gradPhi[13][0] = 0.; // 0
        gradPhi[14][0] = phi[9]; //  z z
        gradPhi[15][0] = 0.; // 0
        gradPhi[16][0] = 2 * phi[14]; // 2 x zz
        gradPhi[17][0] = 0.; // 0
        gradPhi[18][0] = phi[13] + 2 * phi[10]; //  yy z + 2 x y z
        gradPhi[19][0] = phi[15]; // y zz
        gradPhi[20][0] = gradPhi[18][0] * zita; //  yy zz + 2 x y zz
      }


      //phi_y
      gradPhi[0][1] = 0.; // 0
      gradPhi[1][1] = 0.; // 0
      gradPhi[2][1] = 1.;  // 1
      gradPhi[3][1] = 0.; // 0
      gradPhi[4][1] = 0.;  // 0
      gradPhi[5][1] = zita ;  // z

      if(solutionType > 0) {
        gradPhi[4][1] = xi ;  //  x
        gradPhi[5][1] = 0.;  // 0
        gradPhi[6][1] = zita ;  // z
        gradPhi[7][1] = 0.; // 0
        gradPhi[8][1] = 2 * eta; // 2 y
        gradPhi[9][1] = 0.; // 0
        gradPhi[10][1] = phi[5]; // x z
        gradPhi[11][1] = 0.; // 0
        gradPhi[12][1] = 2 * phi[6]; // 2 y z
        gradPhi[13][1] = 0.; // 0
        gradPhi[14][1] = phi[9]; //  z z;
      }

      if(solutionType > 1) {

        gradPhi[11][1] = phi[7] + 2 * phi[4]; // xx + 2 x y
        gradPhi[12][1] = 0.; // 0
        gradPhi[13][1] = 2 * phi[6]; // 2 y z
        gradPhi[14][1] = 0.; // 0
        gradPhi[15][1] = phi[9]; //  z z
        gradPhi[16][1] = 0.; // 0
        gradPhi[17][1] = 2 * phi[15]; // 2 y zz
        gradPhi[18][1] = 2 * phi[10] + phi[12]; // 2 x y z + xx z
        gradPhi[19][1] = phi[14]; // x zz
        gradPhi[20][1] = gradPhi[18][1] * zita; // 2 x y zz + xx zz
      }


      //phi_z
      gradPhi[0][2] = 0.; // 0
      gradPhi[1][2] = 0.; // 0
      gradPhi[2][2] = 0.;  // 0
      gradPhi[3][2] = 1.; // 1
      gradPhi[4][2] = xi ;  // x
      gradPhi[5][2] = eta ;  // y

      if(solutionType > 0) {
        gradPhi[4][2] = 0.;  // 0
        gradPhi[5][2] = xi ;  // x
        gradPhi[6][2] = eta ;  // y
        gradPhi[7][2] = 0.; // 0
        gradPhi[8][2] = 0.; // 0
        gradPhi[9][2] = 2 * zita; // 2 z
        gradPhi[10][2] = phi[4]; // x y
        gradPhi[11][2] = phi[7]; // x x
        gradPhi[12][2] = phi[8]; // y y
        gradPhi[13][2] = 2 * phi[5]; // 2 x z
        gradPhi[14][2] = 2 * phi[6]; // 2 y z;
      }

      if(solutionType > 1) {

        gradPhi[11][2] = 0.; // 0
        gradPhi[12][2] = phi[7]; // x x
        gradPhi[13][2] = phi[8]; // y y
        gradPhi[14][2] = 2  * phi[5]; // 2 x z
        gradPhi[15][2] = 2 * phi[6]; // 2 y z
        gradPhi[16][2] = 2 * phi[12]; // 2 xx z
        gradPhi[17][2] = 2 * phi[13]; // 2 yy z
        gradPhi[18][2] = eta * (phi[4] + phi[7]); // x yy + xx y
        gradPhi[19][2] = 2 * phi[10]; // 2 x y z
        gradPhi[20][2] = 2 * zita * gradPhi[18][2]; // 2 x yy z + 2 xx y z
      }


      //phi_xx
      hessPhi[0][0][0] = 0.; // 0
      hessPhi[1][0][0] = 0.; // 0
      hessPhi[2][0][0] = 0.;  // 0
      hessPhi[3][0][0] = 0.; // 0
      hessPhi[4][0][0] = 0.;  //  0
      hessPhi[5][0][0] = 0.;  // 0

      if(solutionType > 0) {
        hessPhi[4][0][0] = 0.;  //  0
        hessPhi[5][0][0] = 0.;  //  0
        hessPhi[6][0][0] = 0.;  // 0
        hessPhi[7][0][0] = 2.; // 2
        hessPhi[8][0][0] = 0.; // 0
        hessPhi[9][0][0] = 0.; // 0
        hessPhi[10][0][0] = 0.; //  0
        hessPhi[11][0][0] = 2 * zita; // 2 z
        hessPhi[12][0][0] = 0.; // 0
        hessPhi[13][0][0] = 0.; //  0
        hessPhi[14][0][0] = 0.; // 0
      }

      if(solutionType > 1) {

        hessPhi[11][0][0] = 2 * eta; // 2 y
        hessPhi[12][0][0] = 2 * zita; // 2 z
        hessPhi[13][0][0] = 0.; // 0
        hessPhi[14][0][0] = 0.; //  0
        hessPhi[15][0][0] = 0.; // 0
        hessPhi[16][0][0] = 2 * phi[9]; // 2 zz
        hessPhi[17][0][0] = 0.; // 0
        hessPhi[18][0][0] = 2 * phi[6]; // 2 y z
        hessPhi[19][0][0] = 0.; // 0
        hessPhi[20][0][0] = 2 * phi[15]; //  2 y zz
      }


      //phi_xy
      hessPhi[0][1][0] = hessPhi[0][0][1] = 0.; // 0
      hessPhi[1][1][0] = hessPhi[1][0][1] = 0.; // 0
      hessPhi[2][1][0] = hessPhi[2][0][1] = 0.;  // 0
      hessPhi[3][1][0] = hessPhi[3][0][1] = 0.; // 0
      hessPhi[4][1][0] = hessPhi[4][0][1] = 0.;  //  0
      hessPhi[5][1][0] = hessPhi[5][0][1] = 0.;  // 0

      if(solutionType > 0) {
        hessPhi[4][1][0] = hessPhi[4][0][1] = 1.;  //  1
        hessPhi[5][1][0] = hessPhi[5][0][1] = 0.;  //  0
        hessPhi[6][1][0] = hessPhi[6][0][1] = 0.;  // 0
        hessPhi[7][1][0] = hessPhi[7][0][1] = 0.; // 0
        hessPhi[8][1][0] = hessPhi[8][0][1] = 0.; // 0
        hessPhi[9][1][0] = hessPhi[9][0][1] = 0.; // 0
        hessPhi[10][1][0] = hessPhi[10][0][1] = zita ; //  z
        hessPhi[11][1][0] = hessPhi[11][0][1] = 0.; // 0
        hessPhi[12][1][0] = hessPhi[12][0][1] = 0.; // 0
        hessPhi[13][1][0] = hessPhi[13][0][1] = 0.; //  0
        hessPhi[14][1][0] = hessPhi[14][0][1] = 0.; // 0
      }

      if(solutionType > 1) {

        hessPhi[11][1][0] = hessPhi[11][0][1] = 2 * xi + 2 * eta; // 2 x + 2 y
        hessPhi[12][1][0] = hessPhi[12][0][1] = 0.; // 0
        hessPhi[13][1][0] = hessPhi[13][0][1] = 0.; // 0
        hessPhi[14][1][0] = hessPhi[14][0][1] = 0.; //  0
        hessPhi[15][1][0] = hessPhi[15][0][1] = 0.; // 0
        hessPhi[16][1][0] = hessPhi[16][0][1] = 0.; // 0
        hessPhi[17][1][0] = hessPhi[17][0][1] = 0.; // 0
        hessPhi[18][1][0] = hessPhi[18][0][1] = 2 * (phi[6] + phi[5]); // 2 y z + 2 x z
        hessPhi[19][1][0] = hessPhi[19][0][1] = phi[9]; // zz
        hessPhi[20][1][0] = hessPhi[20][0][1] = hessPhi[18][0][1] * zita; // 2 y zz + 2 x zz
      }


      //phi_xz
      hessPhi[0][0][2] = 0.; // 0
      hessPhi[1][0][2] = 0.; // 0
      hessPhi[2][0][2] = 0.;  // 0
      hessPhi[3][0][2] = 0.; // 0
      hessPhi[4][0][2] = 1.;  //  1
      hessPhi[5][0][2] = 0.;  // 0

      if(solutionType > 0) {
        hessPhi[4][0][2] = 0.;  //  0
        hessPhi[5][0][2] = 1.;  //  1
        hessPhi[6][0][2] = 0.;  // 0
        hessPhi[7][0][2] = 0.; // 0
        hessPhi[8][0][2] = 0.; // 0
        hessPhi[9][0][2] = 0.; // 0
        hessPhi[10][0][2] = eta ; //  y
        hessPhi[11][0][2] = 2 * xi; // 2 x
        hessPhi[12][0][2] = 0.; // 0
        hessPhi[13][0][2] = 2 * zita; //  2 z
        hessPhi[14][0][2] = 0.; // 0
      }

      if(solutionType > 1) {

        hessPhi[11][0][2] = 0.; // 0
        hessPhi[12][0][2] = 2 * xi; // 2 x
        hessPhi[13][0][2] = 0.; // 0
        hessPhi[14][0][2] = 2 * zita; //  2 z
        hessPhi[15][0][2] = 0.; // 0
        hessPhi[16][0][2] = 4 * phi[5]; // 4 x z
        hessPhi[17][0][2] = 0.; // 0
        hessPhi[18][0][2] = phi[8] + 2 * phi[4]; //  yy  + 2 x y
        hessPhi[19][0][2] = 2 * phi[6]; // 2 y z
        hessPhi[20][0][2] = 2 * zita * hessPhi[18][0][3]; //  2 yy z + 4 x y z
      }


      //phi_yy
      hessPhi[0][1][1] = 0.; // 0
      hessPhi[1][1][1] = 0.; // 0
      hessPhi[2][1][1] = 0.;  // 0
      hessPhi[3][1][1] = 0.; // 0
      hessPhi[4][1][1] = 0.;  // 0
      hessPhi[5][1][1] = 0.;  // 0

      if(solutionType > 0) {
        hessPhi[4][1][1] = 0.;  //  0
        hessPhi[5][1][1] = 0.;  // 0
        hessPhi[6][1][1] = 0.;  // 0
        hessPhi[7][1][1] = 0.; // 0
        hessPhi[8][1][1] = 2.; // 2
        hessPhi[9][1][1] = 0.; // 0
        hessPhi[10][1][1] = 0.; // 0
        hessPhi[11][1][1] = 0.; // 0
        hessPhi[12][1][1] = 2 * zita; // 2 z
        hessPhi[13][1][1] = 0.; // 0
        hessPhi[14][1][1] = 0.; //  0;
      }

      if(solutionType > 1) {

        hessPhi[11][1][1] = 2 * xi ; // 2 x
        hessPhi[12][1][1] = 0.; // 0
        hessPhi[13][1][1] = 2 * zita; // 2 z
        hessPhi[14][1][1] = 0.; // 0
        hessPhi[15][1][1] = 0.; //  0
        hessPhi[16][1][1] = 0.; // 0
        hessPhi[17][1][1] = 2 * phi[9]; // 2 zz
        hessPhi[18][1][1] = 2 * phi[5]; // 2 x z
        hessPhi[19][1][1] = 0.; // 0
        hessPhi[20][1][1] = 2 * phi[14]; // 2 x zz
      }

      //phi_yz
      hessPhi[0][2][1] = hessPhi[0][1][2] = 0.; // 0
      hessPhi[1][2][1] = hessPhi[1][1][2] = 0.; // 0
      hessPhi[2][2][1] = hessPhi[2][1][2] = 0.;  // 0
      hessPhi[3][2][1] = hessPhi[3][1][2] = 0.; // 0
      hessPhi[4][2][1] = hessPhi[4][1][2] = 0.;  // 0
      hessPhi[5][2][1] = hessPhi[5][1][2] = 1.;  // 1

      if(solutionType > 0) {
        hessPhi[4][2][1] = hessPhi[4][1][2] = 0.;  //  0
        hessPhi[5][2][1] = hessPhi[5][1][2] = 0.;  // 0
        hessPhi[6][2][1] = hessPhi[6][1][2] = 1.;  // 1
        hessPhi[7][2][1] = hessPhi[7][1][2] = 0.; // 0
        hessPhi[8][2][1] = hessPhi[8][1][2] = 0.; // 0
        hessPhi[9][2][1] = hessPhi[9][1][2] = 0.; // 0
        hessPhi[10][2][1] = hessPhi[10][1][2] = xi ; // x
        hessPhi[11][2][1] = hessPhi[11][1][2] = 0.; // 0
        hessPhi[12][2][1] = hessPhi[12][1][2] = 2 * eta; // 2 y
        hessPhi[13][2][1] = hessPhi[13][1][2] = 0.; // 0
        hessPhi[14][2][1] = hessPhi[14][1][2] = 2 * zita; //  2 z;
      }

      if(solutionType > 1) {

        hessPhi[11][2][1] = hessPhi[11][1][2] = 0.; // 0
        hessPhi[12][2][1] = hessPhi[12][1][2] = 0.; // 0
        hessPhi[13][2][1] = hessPhi[13][1][2] = 2 * eta; // 2 y
        hessPhi[14][2][1] = hessPhi[14][1][2] = 0.; // 0
        hessPhi[15][2][1] = hessPhi[15][1][2] = 2 * zita; //  2 z
        hessPhi[16][2][1] = hessPhi[16][1][2] = 0.; // 0
        hessPhi[17][2][1] = hessPhi[17][1][2] = 4 * phi[6]; // 4 y z
        hessPhi[18][2][1] = hessPhi[18][1][2] = 2 * phi[4] + phi[7]; // 2 x y  + xx
        hessPhi[19][2][1] = hessPhi[19][1][2] = 2 * phi[5]; // 2 x z
        hessPhi[20][2][1] = hessPhi[20][1][2] = 2 * zita * hessPhi[18][1][2]; // 4 x y z + 2 xx z
      }


      //phi_zz
      hessPhi[0][2][2] = 0.; // 0
      hessPhi[1][2][2] = 0.; // 0
      hessPhi[2][2][2] = 0.;  // 0
      hessPhi[3][2][2] = 0.; // 0
      hessPhi[4][2][2] = 0.;  // 0
      hessPhi[5][2][2] = 0.;  // 0

      if(solutionType > 0) {
        hessPhi[4][2][2] = 0.;  // 0
        hessPhi[5][2][2] = 0.;  // 0
        hessPhi[6][2][2] = 0.;  // 0
        hessPhi[7][2][2] = 0.; // 0
        hessPhi[8][2][2] = 0.; // 0
        hessPhi[9][2][2] = 2.; // 2
        hessPhi[10][2][2] = 0.; // 0
        hessPhi[11][2][2] = 0.; // 0
        hessPhi[12][2][2] = 0.; // 0
        hessPhi[13][2][2] = 2 * xi; // 2 x
        hessPhi[14][2][2] = 2 * eta; // 2 y ;
      }

      if(solutionType > 1) {

        hessPhi[11][2][2] = 0.; // 0
        hessPhi[12][2][2] = 0.; // 0
        hessPhi[13][2][2] = 0.; // 0
        hessPhi[14][2][2] = 2  * xi; // 2 x
        hessPhi[15][2][2] = 2 * eta; // 2 y
        hessPhi[16][2][2] = 2 * phi[7]; // 2 xx
        hessPhi[17][2][2] = 2 * phi[8]; // 2 yy
        hessPhi[18][2][2] = 0.; // 0
        hessPhi[19][2][2] = 2 * phi[4]; // 2 x y
        hessPhi[20][2][2] = eta * (hessPhi[19][2][2] + hessPhi[16][2][2]) ; // 2 x yy  + 2 xx y
      }

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

      std::vector < std::vector < double > >  hessFm1(dim);

      for(int i1 = 0; i1 < dim; i1++) {
        hessFm1[i1].resize(dim);
      }

      double det = (hessF[0][0] * hessF[1][1] * hessF[2][2] + hessF[0][1] * hessF[1][2] * hessF[2][0] + hessF[0][2] * hessF[1][0] * hessF[2][1])
                   - (hessF[2][0] * hessF[1][1] * hessF[0][2] + hessF[2][1] * hessF[1][2] * hessF[0][0] + hessF[2][2] * hessF[1][0] * hessF[0][1]) ;

      hessFm1[0][0] = (hessF[1][1] * hessF[2][2] - hessF[2][1] * hessF[1][2]) / det ;
      hessFm1[0][1] = (hessF[0][2] * hessF[2][1] - hessF[2][2] * hessF[0][1]) / det ;
      hessFm1[0][2] = (hessF[0][1] * hessF[1][2] - hessF[1][1] * hessF[0][2]) / det ;
      hessFm1[1][0] = (hessF[1][2] * hessF[2][0] - hessF[2][2] * hessF[1][0]) / det ;
      hessFm1[1][1] = (hessF[0][0] * hessF[2][2] - hessF[2][0] * hessF[0][2]) / det ;
      hessFm1[1][2] = (hessF[0][2] * hessF[1][0] - hessF[0][0] * hessF[1][2]) / det ;
      hessFm1[2][0] = (hessF[1][0] * hessF[2][1] - hessF[2][0] * hessF[1][1]) / det ;
      hessFm1[2][1] = (hessF[0][1] * hessF[2][0] - hessF[2][1] * hessF[0][0]) / det ;
      hessFm1[2][2] = (hessF[0][0] * hessF[1][1] - hessF[1][0] * hessF[0][1]) / det ;

      double dxi = hessFm1[0][0] * gradF[0] + hessFm1[0][1] * gradF[1] + hessFm1[0][2] * gradF[2];
      double deta = hessFm1[1][0] * gradF[0] + hessFm1[1][1] * gradF[1] + hessFm1[1][2] * gradF[2];
      double dzita = hessFm1[2][0] * gradF[0] + hessFm1[2][1] * gradF[1] + hessFm1[2][2] * gradF[2];

      xi  -= dxi;
      eta -= deta;
      zita -= dzita;

      if(dxi * dxi + deta * deta + dzita * dzita < 1.0e-6) {
        convergence = true;

        for(int k = 0; k < dim; k++) {
          std::cout << xp[k] << " ";
        }

        std::cout << std::endl;
      }

      //std::cout << "xi = " << xi <<" , "<< "eta = " << eta << " , " << "zita= " << zita <<  std::endl;
    }

    std::vector <double> xcord(3, 0);
    xcord[0] = xi;
    xcord[1] = eta;
    xcord[2] = zita;

    return xcord;

  }

  const double InitialGuess[6][3] = {
    {0., 0., 0.},
    {0.25, 0.25, 0.25},
    {1. / 3., 1. / 3., 0},
    {0., 0.},
    {1. / 3., 1. / 3.},
    {0.}
  };
  
  
  bool SPDCheck2D(std::vector< std::vector <double> > &A){
    bool SPD = true;
    
    if(A[0][1] != A[1][0]){
      std::cout << "ERROR : matrix A is not symmetric" <<std::endl;
      exit(0);
    }
    
    else if(A[0][0] < 0 || fabs(A[0][0]) < 1.0e-8){
      SPD = false;
    }
    else{
      double lambda = A[1][1] - ((A[0][1] * A[0][1]) / A[0][0]);
      if(lambda < 0 || fabs(lambda) < 1.0e-8){
	SPD = false;
      }
    }
    
    return SPD;
    
  }
  

    bool SPDCheck3D(std::vector< std::vector <double> > &A){
    
      bool SPD = true;
      bool notSymm = false;
    
    if(A[0][2] = A[2][0]){
       notSymm = true;
    }
    
    std::vector < std::vector <double> > B;
    B[0][0] = A[0][0];
    B[0][1] = A[0][1];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1];
    
    if(SPDCheck2D(B) == false || notSymm == true){
      std::cout << "ERROR : matrix A is not symmetric" <<std::endl;
      exit(0);
    }
    
    double l00 = sqrt(A[0][0]);
    double l11 = sqrt( A[1][1] - ( (A[0][1] * A[0][1]) / A[0][0] ) );
    double l01 = A[0][1] / l00;
    
    double detL =  l00 * l11;
    std::vector < std::vector < double > > Lm1 ;
    
    Lm1[0][0] = (1/detL) * l11;
    Lm1[1][0] = -(1/detL) * l01;
    Lm1[0][1] = 0. ;
    Lm1[1][1] = (1/detL) * l00 ;
    
    std::vector < double > K;
    
    K[0] = Lm1[0][0] * A[0][2] + Lm1[0][1] * A[1][2];
    K[1] = Lm1[1][0] * A[0][2] + Lm1[1][1] * A[1][2];
    
    double KK = sqrt(K[0]*K[0] + K[1]*K[1]) ;
    
    if(A[2][2] - KK < 0 || fabs(A[2][2] - KK) < 1.0e-8){
      SPD = false;
    }
    
    return SPD;
    
  }
  

  void Marker::InverseMappingTEST(std::vector < double > &x) {
    unsigned dim = _mesh->GetDimension();

    for(int solType = 1; solType < 2; solType++) {
      std::cout << "\n\n--------------------------------------------------" << std::endl;
      std::cout << "solType = " << solType << std::endl;

      //for(int iel = _mesh->_elementOffset[_iproc]; iel < _mesh->_elementOffset[_iproc + 1]; iel++) {
      for(int iel = 394; iel <395; iel++) {
        std::cout << "iel = " << iel << std::endl;
        std::cout << "--------------------------------------------------\n" << std::endl;

        unsigned nDofs = _mesh->GetElementDofNumber(iel, 2);
        short unsigned ielType = _mesh->GetElementType(iel);

        std::vector < std::vector < double > > xv(nDofs);

        for(unsigned i = 0; i < nDofs; i++) {
          xv[i].resize(dim);
          unsigned iDof  = _mesh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof

          for(unsigned k = 0; k < dim; k++) {
            xv[i][k] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
          }
          //std::cout <<"x"<< i+1 <<"="<< xv[i][2]<<"; " << std::endl;
          
        }

        //exit(0);
        for(int j = 0; j < nDofs; j++) {
          std::vector < double > xiT(dim);

          for(unsigned k = 0; k < dim; k++) {
            xiT[k] = *(_mesh->_finiteElement[ielType][0]->GetBasis()->GetXcoarse(j) + k);
          }

          // This is the test
          std::vector < double > xi(dim);
          for(int k = 0; k < dim; k++) {
            xi[k] = InitialGuess[ielType][k];
          }

          for(int k = 0; k < dim; k++) {
            std::cout << "xv[" << k << "]= " << xv[j][k] <<  " ";
          }

          std::cout << std::endl;


          //InverseMapping(iel, 0, xv[j], xi);
          if(solType > 0) {
            std::cout << std::endl;
            InverseMapping(iel, solType, xv[j], xi);
          }

          for(int k = 0; k < dim; k++) {
            std::cout << "xiT[" << k << "]= " << xiT[k] <<  " xi[" << k << "]= " << xi[k];
            std::cout << " error: " << xiT[k] - xi[k] << std::endl;
            if(fabs(xiT[k] - xi[k]) > 1.0e-3) {
              std::cout << "Inverse map test failed " << std::endl;
              abort();
            }
          }
          std::cout << "--------------------------------------------------\n" << std::endl;
        }
      }
    }
  }


//this function returns the position of the marker at time T given the position at time T0 = 0, given the function f and the stepsize h
  std::vector<double> Marker::GetPosition(std::vector<double> (*f)(std::vector<double>), int n, double T) {

    unsigned dim = _mesh->GetDimension();

    unsigned RKOrder = 4;  //order of the RK integration scheme

    std::vector< std::vector<double> > x(RKOrder);
    std::vector< std::vector<double> > K(RKOrder);
    std::vector< double > y(dim, 0);

    // determine the step size
    double h = T / n;

    for(unsigned i = 0; i < RKOrder; i++) {
      x[i].reserve(dim + 1); // x = (t, x, y, z)
      K[i].reserve(dim);
    }
    for(unsigned i = 0; i < RKOrder; i++) {
      x[i].resize(dim + 1);
      K[i].resize(dim);
    }

    //initialize time
    x[0][0] = 0;

    // initialize the position
    for(unsigned i = 1; i < dim + 1; i++) {
      x[0][i] = _x[i - 1] ;
      std::cout << "x[0][" << i << "]= " << x[0][i] << std::endl;
    }


    double step = 0;
    while(step < n) {

      std::cout << "------------------------------------- t = " << x[0][0] << "----------------------------------" << std::endl;
      std::cout << "------------------------------------- h = " << h << "------------------------------------" << std::endl;

      for(unsigned i = 0; i < dim; i++) {
        std::cout << "f(x)[ " << i << "]= " << (*f)(x[0])[i] << std::endl ;
        K[0][i] = h * (*f)(x[0])[i] ;
        std::cout << "K[0][[ " << i << "]= " << K[0][i] << std::endl ;
      }

      //compute x[1]
      x[1][0] = x[0][0] + (0.5 * h);
      for(unsigned j = 1; j < dim + 1; j++) {
        x[1][j] = x[0][j] + 0.5 * K[0][j - 1];
      }

      //compute K[1]
      for(unsigned i = 0; i < dim; i++) {
        K[1][i] = h * (*f)(x[1])[i] ;
        std::cout << "K[1][[ " << i << "]= " << K[1][i] << std::endl ;
      }

      //compute x[2]
      x[2][0] = x[0][0] + (0.5 * h);
      for(unsigned j = 1; j < dim + 1; j++) {
        x[2][j] = x[0][j] + 0.5 * K[1][j - 1];
      }

      //compute K[2]
      for(unsigned i = 0; i < dim; i++) {
        K[2][i] = h * (*f)(x[2])[i] ;
        std::cout << "K[2][[ " << i << "]= " << K[2][i] << std::endl ;
      }

      //compute x[3]
      x[3][0] = x[0][0] + h;
      for(unsigned j = 1; j < dim + 1; j++) {
        x[3][j] = x[0][j] + K[2][j - 1];
      }

      //compute K[3]
      for(unsigned i = 0; i < dim; i++) {
        K[3][i] = h * (*f)(x[3])[i] ;
        std::cout << "K[3][[ " << i << "]= " << K[3][i] << std::endl ;
      }

      // RK stepping
      for(unsigned j = 1; j < dim + 1; j++) {
        x[0][j] += (1. / 6) * (K[0][j - 1] + 2. * K[1][j - 1] + 2. * K[2][j - 1] + K[3][j - 1]);
        std::cout << "x[0][" << j << "]=" << x[0][j] << std::endl;
      }
      //update t
      x[0][0] += h ;

      //update the step
      step++;
    }

    for(unsigned j = 1; j < dim + 1; j++) {
      y[j - 1] = x[0][j];
    }

    return y;

  }


///////////////////////////////////////////// Here there are the OLD functions used for the 2D and 3D inclusion test ///////////////////////////////////////////////



  unsigned Marker::GetNextElement2DOLD(const unsigned &dim, const unsigned &currentElem, const unsigned &previousElem) {

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
        else if(fabs(xv[1][i + 1]) < epsilon && xv[0][i + 1] > 0) {  // the second vertex is on the positive x-axis
          //std::cout << "the second vertex is on the positive x-axis" << std::endl;
          if(xv[1][i] < 0) w += .5;
          else w -= .5;
        }
        else if(xv[1][i]*xv[1][i + 1] < 0) {  // the edge crosses the x-axis but doesn't pass through the origin
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
        else if(fabs(xv[0][i + 1]) < epsilon && fabs(xv[1][i + 1]) < epsilon) {  // vertex 2 is the origin
          w = 1; // set to 1 by default
          //std::cout << "w set to 1 by default (vertex 2 is in the origin)" << std::endl;
          break;
        }
        else if(xv[0][i] * xv[0][i + 1] < 0 || xv[1][i] * xv[1][i + 1] < 0) {  //the edge crosses the origin
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



  unsigned Marker::GetNextElement3DOLD(const unsigned &dim, const unsigned &currentElem, const unsigned &previousElem) {

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


        if(fabs(A) < epsilon && epsilon < fabs(By0 + Cz0)) {  // A = 0 and By0 != -Cz0
          // std::cout << "The plane of face " << iface << "and the x-axis don't intersect " << std::endl;

        }
        else {

          std::vector < double > xTilde(dim, 0);
          double r = 0;

          if(fabs(A) < epsilon && fabs(By0 + Cz0) < epsilon) {  // A = 0 and By0 = -Cz0
            //   std::cout << "The plane of face " << iface << "and the x-axis intersect on a line" << std::endl;
            // the marker and the face are already on the same plane so there is no need for further shifting
            lineIntersection = true ;
          }

          else if(epsilon < fabs(A)) {  // A != 0
            // std::cout << "The plane of face " << iface << "and the x-axis intersect at a point" << std::endl;

            r = (A * xv[0][0] + B * xv[1][0] + C * xv[2][0]) / A;
            xTilde[0] = r;

            if(fabs(r) < epsilon) {  // r = 0
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
//             for(int i = 0; i < 3; i ++) {
//                 unsigned j = _mesh->GetLocalFaceVertexIndex(currentElem, iface, i);
//                 unsigned iDof = _mesh->GetSolutionDof(j, currentElem, 2);

        //  for(unsigned k = 0; k < dim; k++) {
        //std::cout << (*_mesh->_topology->_Sol[k])(iDof)  << " ";
        // xg[k] += 1. / 3.*(*_mesh->_topology->_Sol[k])(iDof); ///COS 'E' QUESTO ??
        //  }
        //std::cout << std::endl;
        //}
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




}





