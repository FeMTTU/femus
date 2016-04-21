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

namespace femus {

  void Marker::GetElement(const bool &debug) {

    unsigned dim = _mesh->GetDimension();

    std::vector < unsigned > processorMarkerFlag(_nprocs, 3);

    //BEGIN SMART search
    // look to the closest element among a restricted list
    double modulus = 1.e10;
    int iel = _mesh->_elementOffset[_iproc + 1];
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
    if(debug){
      if(iel == _mesh->_elementOffset[_iproc + 1]) {
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
   
    std::vector< int > nextElem(_nprocs, -1);
    int previousElem = iel;
    unsigned nextProc = _iproc;
    
    while(!elementHasBeenFound) {
      
      //BEGIN next element search 
      while(elementHasBeenFound + pointIsOutsideThisProcess + pointIsOutsideTheDomain == 0) {
        if(dim == 2) {
          nextElem[_iproc] = GetNextElement2D(dim, iel, previousElem);
        }
        else if(dim == 3) {
          nextElem[_iproc] = GetNextElement3D(dim, iel, previousElem);
        }
        previousElem = iel;

        if(nextElem[_iproc] == iel) {
          _elem = iel;
          elementHasBeenFound = true;
          processorMarkerFlag[_iproc] = 1;
        }
        else if(nextElem[_iproc] < 0) {
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

      if( debug ){
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
      if( sumFlag == 0 ) { // all the processes beleive that the marker is outside the domain
        std::cout << "Marker is outside the domain" << std::endl;
        _elem = UINT_MAX;
	if(debug){
	  for(unsigned j = 0 ; j < _nprocs; j++) {
	    std::cout << " processorMarkerFlag[" << j << "] = " << processorMarkerFlag[j] <<  std::endl;
	  }	
	}
        break;
      }

      // _iproc sends its nextElem (which is in jproc) to jproc
      if(!elementHasBeenFound) {
        if(processorMarkerFlag[_iproc] == 2) {
          MPI_Send(&nextElem[_iproc], 1, MPI_INT, nextProc, 1 , PETSC_COMM_WORLD);
        }
        for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
          if(processorMarkerFlag[jproc] == 3) {
            MPI_Recv(&nextElem[jproc], 1, MPI_INT, jproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
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
          previousElem = iel;
        }
        else {
	  pointIsOutsideTheDomain = true;
	  pointIsOutsideThisProcess = true;
          processorMarkerFlag[_iproc] = 0;
        }
      }
      if(debug){
	for(unsigned j = 0 ; j < _nprocs; j++) {
	  std::cout << " processorMarkerFlag[" << j << "] = " << processorMarkerFlag[j] << "  " << nextElem[j] << std::endl;
	}
      }
      //END process exchange  
  
    }
  }


  int Marker::GetNextElement2D(const unsigned &dim, const int &currentElem, const int &previousElem) {

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
//       if(currentElem == 125) {
//         std::cout << "Delta for element" << currentElem << " is =" << Delta  << " , " << xv[0][i] << " , " << xv[1][i] << " , " << xv[0][i + 1] << " , " << xv[1][i + 1] << std::endl;
//       }

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
      std::cout << " Error negative Winding Number with counterclockwise oriented points " << std::endl;
      abort();
    }

    std::cout << "the next element is " << nextElem << std::endl;


    return nextElem;
  }

  int Marker::GetNextElement3D(const unsigned &dim, const int &currentElem, const int &previousElem) {

    int nextElem;
    //for on the element faces
    int faceIntersectionCounter = 0; // if it is even, the marker is outside the element (0 is considered even)
    bool markerIsInElement = false;


    for(unsigned iface = 0; iface < _mesh->GetElementFaceNumber(currentElem); iface++) {
      std::cout << "faceIntersectionCounter  = " << faceIntersectionCounter  << " , " << " markerIsInElement = " << markerIsInElement <<  std::endl;

      bool lineIntersection = false;

      std::vector< std::vector < double > > xv(dim); //stores the coordinates of the nodes of the face iface

      for(unsigned k = 0; k < dim; k++) {
        xv[k].reserve(9);
      }
      short unsigned currentElementType = _mesh->GetElementType(currentElem);
      short unsigned ifaceType = _mesh->GetElementFaceType(currentElem, iface);
      for(unsigned k = 0; k < dim; k++) {
        xv[k].resize(facePointNumber[ifaceType]);
      }
      for(unsigned i = 0; i < facePointNumber[ifaceType]; i++) {
        unsigned ifaceDof  = _mesh->GetSolutionDof(facePoints3D[currentElementType][iface][i], currentElem, 2);
        // std::cout << "ifaceDof = " << ifaceDof << std::endl;
        for(unsigned k = 0; k < dim; k++) {
          xv[k][i] = (*_mesh->_topology->_Sol[k])(ifaceDof) - _x[k];     // global extraction and local storage for the element coordinates
        }
      }


      double length = 0.;
      for(unsigned i = 0; i < xv[0].size() - 1; i++) {
        length += sqrt((xv[0][i + 1] - xv[0][i]) * (xv[0][i + 1] - xv[0][i]) +
                       (xv[1][i + 1] - xv[1][i]) * (xv[1][i + 1] - xv[1][i]) +
                       (xv[2][i + 1] - xv[2][i]) * (xv[2][i + 1] - xv[2][i]));
      }

      length /= xv[0].size();

      for(unsigned i = 0; i < xv[0].size(); i++) {
        xv[0][i] /= length;
        xv[1][i] /= length;
        xv[2][i] /= length;
      }

      // now we have to check if the plane of the face intersects the positive x-axis
      // let's find the plane passing through xv[][0], xv[][2] and xv[][4] (they will not be aligned)

      // entries of the normal to the plane
      double A = -(xv[1][4] - xv[1][0]) * (xv[2][2] - xv[2][0]) + (xv[2][4] - xv[2][0]) * (xv[1][2] - xv[1][0]);
      double B = -(xv[2][4] - xv[2][0]) * (xv[0][2] - xv[0][0]) + (xv[0][4] - xv[0][0]) * (xv[2][2] - xv[2][0]);
      double C = -(xv[0][4] - xv[0][0]) * (xv[1][2] - xv[1][0]) + (xv[1][4] - xv[1][0]) * (xv[0][2] - xv[0][0]);

      double epsilon  = 10e-10;
      double epsilon2  = epsilon * epsilon;
      double By0 = B * xv[1][0];
      double Cz0 = C * xv[2][0];

      // std::cout << "A = " << A << " , " << "By0 = " << By0 << " , " << " Cz0 = " << Cz0 <<  std::endl;


      if(fabs(A) < epsilon && epsilon < fabs(By0 + Cz0)) { // A = 0 and By0 != -Cz0
        std::cout << "The plane of face " << iface << "and the x-axis don't intersect " << std::endl;

      }
      else {

        std::vector < double > xTilde(dim, 0);
        double r = 0;

        if(fabs(A) < epsilon && fabs(By0 + Cz0) < epsilon) { // A = 0 and By0 = -Cz0
          std::cout << "The plane of face " << iface << "and the x-axis intersect on a line" << std::endl;
          // the marker and the face are already on the same plane so there is no need for further shifting
          lineIntersection = true ;
        }

        else if(epsilon < fabs(A)) {  // A != 0
          std::cout << "The plane of face " << iface << "and the x-axis intersect at a point" << std::endl;

          r = (A * xv[0][0] + B * xv[1][0] + C * xv[2][0]) / A;
          xTilde[0] = r;
          if(fabs(r) < epsilon) { // r = 0
            lineIntersection = true; // the intersection point is the actual marker
          }
        }
        if(r > 0 || fabs(r) < epsilon) {

          for(unsigned i = 0; i < facePointNumber[ifaceType]; i++) {
            for(unsigned k = 0; k < dim; k++) {
              xv[k][i] = xv[k][i] - xTilde[k];     //transate again the reference frame so that the origin is xTilde
            }
          }

          unsigned scalarCount = 0;

          for(unsigned i = 0; i < xv[0].size() - 1; i++) {
            //entries of the vector (xTilde - xi) X ( xi+1 -xi)
            double q0 = xv[1][i] * (xv[2][i] - xv[2][i + 1]) + xv[2][i] * (xv[1][i + 1] - xv[1][i]);
            double q1 = xv[2][i] * (xv[0][i] - xv[0][i + 1]) + xv[0][i] * (xv[2][i + 1] - xv[2][i]);
            double q2 = xv[0][i] * (xv[1][i] - xv[1][i + 1]) + xv[1][i] * (xv[0][i + 1] - xv[0][i]);

            // std::cout << "q0 = " << q0 << " , " << "q1 = " << q1 << " , " << " q2 = " << q2 <<  std::endl;

            double  scalarProduct = q0 * A + q1 * B + q2 * C;

            // std::cout << "scalarProduct = " << scalarProduct << std::endl;

            if(scalarProduct > 0) {  // this can be cancelled once the code is working ok
              std::cout << "the marker is outside face " << iface <<  std::endl;
              break;

            }
            else if(fabs(scalarProduct) < epsilon) {  //scalarProduct == 0
              std::cout << " the marker and the edge are aligned " << std::endl; //check if xTilde is actually on the edge.
              if(xv[0][i]*xv[0][i + 1] < 0 || xv[1][i]*xv[1][i + 1] < 0 || xv[2][i]*xv[2][i + 1] < 0) {
                if(lineIntersection == true) {
                  std::cout << " the marker belongs to an edge of face " << iface << std::endl;
                  markerIsInElement = true;
                  break;
                }
                else {
                  faceIntersectionCounter++;
                  std::cout << " faceIntersectionCounter = " << faceIntersectionCounter << std::endl;
                  break;
                }

              }
              else if((xv[0][i] * xv[0][i]  + xv[1][i] * xv[1][i] + xv[2][i] * xv[2][i]) < epsilon2 ||
                      (xv[0][i + 1]*xv[0][i + 1] + xv[1][i + 1]*xv[1][i + 1] + xv[2][i + 1]*xv[2][i + 1]) < epsilon2) {
                if(lineIntersection == true) {
                  std::cout << " one of the vertexes is the marker" << std::endl;
                  markerIsInElement = true;
                  break;
                }
                else {
                  faceIntersectionCounter++;
                  std::cout << " faceIntersectionCounter " << faceIntersectionCounter << std::endl;
                  break;
                }
              }
            }
            else if(scalarProduct < 0) {
              scalarCount++;
            }
//                     std::cout<<"BBBBBBBBBBBBBBBB "<<scalarCount<<std::endl;
          }

//                 std::cout<<"BBBBBBBBBBBBBBBB "<<scalarCount<<std::endl;

          if(scalarCount == facePointNumber[ifaceType] - 1 && lineIntersection == true) {
            markerIsInElement = true ;
            break;
          }
          else if(scalarCount == facePointNumber[ifaceType] - 1 && lineIntersection == false) {
            faceIntersectionCounter++;
          }
        }

      } //end of else
      if(markerIsInElement == true) {
        break;
      }
    }// end of the loop on iface

    std::cout << "markerIsInElement = " << markerIsInElement << " and faceIntersectionCounter  = " << faceIntersectionCounter  <<  std::endl;

    if(markerIsInElement == true || faceIntersectionCounter % 2 != 0) {
      nextElem = currentElem;
    }
    else if(markerIsInElement == false && faceIntersectionCounter % 2 == 0) {
      std::cout << " The marker doesn't belong to element " << currentElem << std::endl;
      double modulus = 1.e10;

      for(unsigned iface = 0; iface < _mesh->GetElementFaceNumber(currentElem); iface++) {

        unsigned faceNodeNumber = _mesh->GetElementFaceDofNumber(currentElem, iface, 2);
        unsigned i = _mesh->GetLocalFaceVertexIndex(currentElem, iface, faceNodeNumber - 1);
        std::cout << i << " ";

        unsigned faceCentralDof = _mesh->GetSolutionDof(i, currentElem, 2);

        //std::cout << "faceCentralDof =" << faceCentralDof << std::endl;
        double distance2 = 0;
        for(unsigned k = 0; k < dim; k++) {
          double dk = (*_mesh->_topology->_Sol[k])(faceCentralDof) - _x[k];     // global extraction and local storage for the element coordinates
          distance2 += dk * dk;
        }
        double ifaceModulus = sqrt(distance2);

        if(ifaceModulus < modulus) {
          int jel = (_mesh->el->GetFaceElementIndex(currentElem, iface) - 1);
          std::cout << "jel = " << jel << "iface = " << iface <<  std::endl;
          if(jel != previousElem) {
            nextElem = jel;
            modulus = ifaceModulus;
          }
        }
      }
    }
    std::cout << "nextElem = " << nextElem << std::endl;

    return nextElem;
  }
}




