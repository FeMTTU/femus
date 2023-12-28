/*=========================================================================

 Program: FEMUS
 Module: Marker
 Authors: Eugenio Aulisa, Giacomo Capodaglio

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Marker.hpp"
#include "Line.hpp"
#include "NumericVector.hpp"
#include "stationary_or_time_dep.hpp"

#include <cmath>

#include "PolynomialBases.hpp"




const double faceNormal[6][6][3] = {
  {{0., -1., 0.} , {1., 0., 0.}, {0., 1., 0.} , { -1., 0., 0.}, {0., 0., -1.} , {0., 0., 1.}},
  {{0., 0., -1}, {0., -1., 0.}, {1. / sqrt(3.), 1. / sqrt(3.), 1. / sqrt(3.)}, { -1., 0., 0.}},
  {{0., -1., 0.}, {1. / sqrt(2.), 1. / sqrt(2.), 0.}, { -1., 0., 0.}, {0., 0., -1}, {0., 0., 1}},
  {{0., -1.} , {1., 0.}, {0., 1.} , { -1., 0.}},
  {{0., -1.}, {1. / sqrt(2.), 1. / sqrt(2.)}, { -1., 0.}},
  {{ -1}, {1}}
};

const unsigned faceNumber[6] = {6, 4, 5, 4, 3, 2};

unsigned counter = 0;

const unsigned facePointNumber[6][3] = {{}, {}, {}, {5, 9, 9}, {4, 7, 7}, {}}; //facePointNumber[elemType][solType]

const unsigned facePoints[6][3][9] = { //facePointNumber[elemType][solType][localFaceIndex]
  {{ }},
  {{ }},
  {{ }},
  { { 0, 1, 2, 3, 0} , { 0, 4, 1, 5, 2, 6, 3, 7, 0 } , { 0, 4, 1, 5, 2, 6, 3, 7, 0 } },
  { {0, 1, 2, 0} , { 0, 3, 1, 4, 2, 5, 0 } , { 0, 3, 1, 4, 2, 5, 0 } },
  { }
};

const unsigned trianglesPerFace[6][3][6] = { //trainglesPerFace[elementType][solType][numberOfTrianglesPerFace]
  {{4, 4, 4, 4, 4, 4}, {8, 8, 8, 8, 8, 8}, {8, 8, 8, 8, 8, 8}},
  {{1, 1, 1, 1}, {6, 6, 6, 6}, {6, 6, 6, 6}},
  {{4, 4, 4, 1, 1}, {8, 8, 8, 6, 6}, {8, 8, 8, 6, 6}},
  {{4}, {8}, {8}},
  {{1}, {6}, {6}},
  {}
};

const unsigned faceTriangleNodes[6][3][6][8][4] = { // elementType - solType - number of faces - number of triangles per faces - 3 vertices of the triangle (the first is counted twice)

  { { {{0, 1, 20, 0}, {1, 5, 20, 1}, {5, 4, 20, 5}, {4, 0, 20, 4}},  //hex
      {{1, 2, 21, 1}, {2, 6, 21, 2}, {6, 5, 21, 6}, {5, 1, 21, 5}},
      {{2, 3, 22, 2}, {3, 7, 22, 3}, {7, 6, 22, 7}, {6, 2, 22, 6}},
      {{3, 0, 23, 3}, {0, 4, 23, 0}, {4, 7, 23, 4}, {7, 3, 23, 7}},
      {{0, 3, 24, 0}, {3, 2, 24, 3}, {2, 1, 24, 2}, {1, 0, 24, 1}},
      {{4, 5, 25, 4}, {5, 6, 25, 5}, {6, 7, 25, 6}, {7, 4, 25, 7}}
    },
    { {{0, 8, 20, 0}, {8, 1, 20, 8}, {1, 17, 20, 1}, {17, 5, 20, 17}, {5, 12, 20, 5}, {12, 4, 20, 12}, {4, 16, 20, 4 }, {16, 0, 20, 16}},
      {{1, 9, 21, 1}, {9, 2, 21, 9}, {2, 18, 21, 2}, {18, 6, 21, 18}, {6, 13, 21, 6}, {13, 5, 21, 13}, {5, 17, 21, 5}, {17, 1, 21, 17}},
      {{2, 10, 22, 2}, {10, 3, 22, 10}, {3, 19, 22, 3}, {19, 7, 22, 19}, {7, 14, 22, 7}, {14, 6, 22, 14}, {6, 18, 22, 6}, {18, 2, 22, 18}},
      {{3, 11, 23, 3}, {11, 0, 23, 11}, {0, 16, 23, 0}, {16, 4, 23, 16}, {4, 15, 23, 4}, {15, 7, 23, 15}, {7, 19, 23, 7}, {19, 3, 23, 19}},
      {{0, 11, 24, 0}, {11, 3, 24, 11}, {3, 10, 24, 3}, {10, 2, 24, 10}, {2, 9, 24, 2}, {9, 1, 24, 9}, {1, 8, 24, 1}, {8, 0, 24, 8}},
      {{4, 12, 25, 4}, {12, 5, 25, 12}, {5, 13, 25, 5}, {13, 6, 25, 13}, {6, 14, 25, 6}, {14, 7, 25, 14}, {7, 15, 25, 7}, {15, 4, 25, 15}}
    },
    { {{0, 8, 20, 0}, {8, 1, 20, 8}, {1, 17, 20, 1}, {17, 5, 20, 17}, {5, 12, 20, 5}, {12, 4, 20, 12}, {4, 16, 20, 4 }, {16, 0, 20, 16}},
      {{1, 9, 21, 1}, {9, 2, 21, 9}, {2, 18, 21, 2}, {18, 6, 21, 18}, {6, 13, 21, 6}, {13, 5, 21, 13}, {5, 17, 21, 5}, {17, 1, 21, 17}},
      {{2, 10, 22, 2}, {10, 3, 22, 10}, {3, 19, 22, 3}, {19, 7, 22, 19}, {7, 14, 22, 7}, {14, 6, 22, 14}, {6, 18, 22, 6}, {18, 2, 22, 18}},
      {{3, 11, 23, 3}, {11, 0, 23, 11}, {0, 16, 23, 0}, {16, 4, 23, 16}, {4, 15, 23, 4}, {15, 7, 23, 15}, {7, 19, 23, 7}, {19, 3, 23, 19}},
      {{0, 11, 24, 0}, {11, 3, 24, 11}, {3, 10, 24, 3}, {10, 2, 24, 10}, {2, 9, 24, 2}, {9, 1, 24, 9}, {1, 8, 24, 1}, {8, 0, 24, 8}},
      {{4, 12, 25, 4}, {12, 5, 25, 12}, {5, 13, 25, 5}, {13, 6, 25, 13}, {6, 14, 25, 6}, {14, 7, 25, 14}, {7, 15, 25, 7}, {15, 4, 25, 15}}
    }
  },
  { { {{0, 2, 1, 0}},   //tet
      {{0, 1, 3, 0}},
      {{1, 2, 3, 1}},
      {{2, 0, 3, 2}}
    },
    { {{0, 6, 10, 0}, {6, 2, 10, 6}, {2, 5, 10, 2}, {5, 1, 10, 5}, {1, 4, 10, 1}, {4, 0, 10, 4}},
      {{0, 4, 11, 0}, {4, 1, 11, 4}, {1, 8, 11, 1}, {8, 3, 11, 8}, {3, 7, 11, 3}, {7, 0, 11, 7}},
      {{1, 5, 12, 1}, {5, 2, 12, 5}, {2, 9, 12, 2}, {9, 3, 12, 9}, {3, 8, 12, 3}, {8, 1, 12, 8}},
      {{2, 6, 13, 2}, {6, 0, 13, 6}, {0, 7, 13, 0}, {7, 3, 13, 7}, {3, 9, 13, 3}, {9, 2, 13, 9}}
    },
    { {{0, 6, 10, 0}, {6, 2, 10, 6}, {2, 5, 10, 2}, {5, 1, 10, 5}, {1, 4, 10, 1}, {4, 0, 10, 4}},
      {{0, 4, 11, 0}, {4, 1, 11, 4}, {1, 8, 11, 1}, {8, 3, 11, 8}, {3, 7, 11, 3}, {7, 0, 11, 7}},
      {{1, 5, 12, 1}, {5, 2, 12, 5}, {2, 9, 12, 2}, {9, 3, 12, 9}, {3, 8, 12, 3}, {8, 1, 12, 8}},
      {{2, 6, 13, 2}, {6, 0, 13, 6}, {0, 7, 13, 0}, {7, 3, 13, 7}, {3, 9, 13, 3}, {9, 2, 13, 9}}
    }
  },
  { { {{0, 1, 15, 0}, {1, 4, 15, 1}, {4, 3, 15, 4}, {3, 0, 15, 3}},  //wedge
      {{1, 2, 16, 1}, {2, 5, 16, 2}, {5, 4, 16, 5}, {4, 1, 16, 4}},
      {{2, 0, 17, 2}, {0, 3, 17, 0}, {3, 5, 17, 3}, {5, 2, 17, 5}},
      {{0, 2, 1, 0}},
      {{3, 4, 5, 3}}
    },
    { {{0, 6, 15, 0 }, {6, 1, 15, 6}, {1, 13, 15, 1}, {13, 4, 15, 13}, {4, 9, 15, 4}, {9, 3, 15, 9}, {3, 12, 15, 3}, {12, 0, 15, 12}},
      {{1, 7, 16, 1}, {7, 2, 16, 7}, {2, 14, 16, 2}, {14, 5, 16, 14}, {5, 10, 16, 5}, {10, 4, 16, 10}, {4, 13, 16, 4}, {13, 1, 16, 13}},
      {{2, 8, 17, 2}, {8, 0, 17, 8}, {0, 12, 17, 0}, {12, 3, 17, 12}, {3, 11, 17, 3}, {11, 5, 17, 11}, {5, 14, 17, 5}, {14, 2, 17, 14}},
      {{0, 8, 18, 0}, {8, 2, 18, 8}, {2, 7 , 18, 2}, {7, 1, 18, 7} , {1, 6 , 18, 1} , {6, 0, 18, 6}},
      {{3, 9, 19, 3}, {9, 4, 19, 9}, {4, 10, 19, 4}, {10, 5, 19, 10}, {5, 11, 19, 5}, {11, 3, 19, 11}}
    },
    { {{0, 6, 15, 0 }, {6, 1, 15, 6}, {1, 13, 15, 1}, {13, 4, 15, 13}, {4, 9, 15, 4}, {9, 3, 15, 9}, {3, 12, 15, 3}, {12, 0, 15, 12}},
      {{1, 7, 16, 1}, {7, 2, 16, 7}, {2, 14, 16, 2}, {14, 5, 16, 14}, {5, 10, 16, 5}, {10, 4, 16, 10}, {4, 13, 16, 4}, {13, 1, 16, 13}},
      {{2, 8, 17, 2}, {8, 0, 17, 8}, {0, 12, 17, 0}, {12, 3, 17, 12}, {3, 11, 17, 3}, {11, 5, 17, 11}, {5, 14, 17, 5}, {14, 2, 17, 14}},
      {{0, 8, 18, 0}, {8, 2, 18, 8}, {2, 7 , 18, 2}, {7, 1, 18, 7} , {1, 6 , 18, 1} , {6, 0, 18, 6}},
      {{3, 9, 19, 3}, {9, 4, 19, 9}, {4, 10, 19, 4}, {10, 5, 19, 10}, {5, 11, 19, 5}, {11, 3, 19, 11}}
    }
  },
  { { {{0, 1, 8, 0}, {1, 2, 8, 1}, {2, 3, 8, 1}, {3, 0, 8, 3}}, //quad
      {{0, 4, 8, 0}, {4, 1, 8, 4}, {1, 5, 8, 1}, {5, 2, 8, 5}, {2, 6, 8, 2}, {6, 3, 8, 6}, {3, 7, 8, 3}, {7, 0, 8, 7}},
      {{0, 4, 8, 0}, {4, 1, 8, 4}, {1, 5, 8, 1}, {5, 2, 8, 5}, {2, 6, 8, 2}, {6, 3, 8, 6}, {3, 7, 8, 3}, {7, 0, 8, 7}}
    }
  },
  { { {{0, 1, 2, 0}},
      {{0, 3, 6, 0}, {3, 1, 6, 3}, {1, 4, 6, 1}, {4, 2, 6, 4}, {2, 5, 6, 2}, {5, 0, 6, 5}},
      {{0, 3, 6, 0}, {3, 1, 6, 3}, {1, 4, 6, 1}, {4, 2, 6, 4}, {2, 5, 6, 2}, {5, 0, 6, 5}}
    }
  },
  {{{{}}}},
};

namespace femus
{

  const double Marker::_a[4][4][4] = {
    { {} // first order
    },
    { {}, // second order (Heun's)
      {1., 0.}
    },
    { {}, // third-order method
      {0.5},
      { -1., 2}
    },
    { {}, // fourth-order method
      {0.5},
      {0., 0.5},
      {0., 0., 1.}
    }
  };

  const double Marker::_b[4][4] = {
    {1.}, // first order
    {0.5, 0.5}, // second order (Heun's)
    {1. / 6., 2. / 3., 1. / 6.}, // third-order method
    {1. / 6., 1. / 3., 1. / 3., 1. / 6.} // fourth-order method
  };

  const double Marker::_c[4][4] = {
    {0.}, // first order
    {0., 1.}, // second order (Heun's)
    {0., 0.5, 1.}, // third-order method
    {0., 0.5, 0.5, 1.} // fourth-order method
  };

  const double Marker::_localCentralNode[6][3] = {
    {0., 0., 0.},
    {0.25, 0.25, 0.25},
    {1. / 3., 1. / 3., 0},
    {0., 0.},
    {1. / 3., 1. / 3.},
    {0.}
  };

  void Marker::GetElement(const bool &useInitialSearch, const unsigned &initialElem, const Solution* sol, const double &s)
  {

    std::vector < unsigned > processorMarkerFlag(_nprocs, 3);

    unsigned iel = initialElem;

    unsigned ielProc = (initialElem < sol->GetMesh()->GetElementOffset(_nprocs) ) ?
                       sol->GetMesh()->BisectionSearch_find_processor_of_dof(iel, 3) : _nprocs;

    if (useInitialSearch || _iproc != ielProc) {

      //BEGIN SMART search
      // look to the closest element among a restricted list

      iel = UINT_MAX;

      double modulus = 1.e10;

      for (int jel = sol->GetMesh()->GetElementOffset(_iproc); jel < sol->GetMesh()->GetElementOffset(_iproc + 1); jel += 25) {

        unsigned interiorNode = sol->GetMesh()->GetElementDofNumber(jel, 2) - 1;
        unsigned jDof  = sol->GetMesh()->GetSolutionDof(interiorNode, jel, 2);    // global to global mapping between coordinates node and coordinate dof

        double distance2 = 0;

        for (unsigned k = 0; k < _dim; k++) {
          double dk = GetCoordinates(sol, k, jDof, s) - _x[k];   // global extraction and local storage for the element coordinates
          distance2 += dk * dk;
        }

        double modulusKel = sqrt(distance2);

        if (modulusKel < modulus) {
          iel = jel;
          modulus = modulusKel;
        }
      }

      if (iel == UINT_MAX) {
        //  std::cout << "Warning the marker is located on unreasonable distance from the mesh >= 1.e10" << std::endl;
      }
      else {
        //  std::cout << "the smart search starts from element " << iel << std::endl;
      }


      //END SMART search:
    }


    bool elementHasBeenFound = false;
    bool pointIsOutsideThisProcess = false;
    bool pointIsOutsideTheDomain = false;

    std::vector< unsigned > nextElem(_nprocs, UINT_MAX);
    std::vector< unsigned > previousElem(_nprocs, UINT_MAX);
    previousElem[_iproc] = iel;
    unsigned nextProc = _iproc;

    while (!elementHasBeenFound) {

      //BEGIN next element search
      while (elementHasBeenFound + pointIsOutsideThisProcess + pointIsOutsideTheDomain == 0) {
        if (_dim == 2) {
          nextElem[_iproc] = GetNextElement2D(iel, previousElem[_iproc], sol, s);
        }
        else if (_dim == 3) {
          nextElem[_iproc] = GetNextElement3D(iel, previousElem[_iproc], sol, s);
        }

        previousElem[_iproc] = iel;

        if (nextElem[_iproc] == iel) {
          _elem = iel;
          elementHasBeenFound = true;
          processorMarkerFlag[_iproc] = 1;
        }
        else if (nextElem[_iproc] == UINT_MAX) {
          pointIsOutsideTheDomain = true;
          processorMarkerFlag[_iproc] = 0;
        }
        else {
          nextProc = sol->GetMesh()->BisectionSearch_find_processor_of_dof(nextElem[_iproc], 3);

          if (nextProc != _iproc) {
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


      if (elementHasBeenFound) {
        //  std::cout << " The marker belongs to element " << _elem << std::endl;
      }
      else if (pointIsOutsideTheDomain) {
        //  std::cout << " The marker does not belong to this domain" << std::endl;
      }
      else if (pointIsOutsideThisProcess) {
        //  std::cout << "proc " << _iproc << " believes the marker is in proc = " << nextProc << std::endl;
      }


      //END next element search

      //BEGIN process exchange
      //send/receive if any process found the element
      for (unsigned jproc = 0; jproc < _nprocs; jproc++) {
        if (jproc != _iproc) {
          if (processorMarkerFlag[_iproc] == 2 && jproc == nextProc) {
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

      for (unsigned i = 0; i < _nprocs; i++) {
        sumFlag += processorMarkerFlag[i];

        if (processorMarkerFlag[i] == 1) {
          _mproc = i;
          elementHasBeenFound = true;
          break;
        }
      }

      if (sumFlag == 0) {  // all the processes beleive that the marker is outside the domain
        // std::cout << "Marker is outside the domain" << std::endl;
        _mproc = _nprocs;
        _elem = UINT_MAX;


//           for(unsigned j = 0 ; j < _nprocs; j++) {
//             std::cout << " processorMarkerFlag[" << j << "] = " << processorMarkerFlag[j] <<  std::endl;
//           }


        break;
      }

      // _iproc sends its nextElem (which is in jproc) to jproc
      if (!elementHasBeenFound) {
        if (processorMarkerFlag[_iproc] == 2) {
          MPI_Send(&nextElem[_iproc], 1, MPI_UNSIGNED, nextProc, 1 , PETSC_COMM_WORLD);
          MPI_Send(&previousElem[_iproc], 1, MPI_UNSIGNED, nextProc, 2 , PETSC_COMM_WORLD);
        }

        for (unsigned jproc = 0; jproc < _nprocs; jproc++) {
          if (processorMarkerFlag[jproc] == 3) {
            MPI_Recv(&nextElem[jproc], 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&previousElem[jproc], 1, MPI_UNSIGNED, jproc, 2 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          }
        }


//           for(unsigned j = 0 ; j < _nprocs; j++) {
//             std::cout << " processorMarkerFlag[" << j << "] = " << processorMarkerFlag[j] << "  " << nextElem[j] << std::endl;
//           }



        //BEGIN SMART search among all the elements received by _iproc and sent from jprocs
        double modulus = 1.e10;
        iel = sol->GetMesh()->GetElementOffset(_iproc + 1) ;

        for (unsigned jproc = 0; jproc < _nprocs; jproc++) {
          if (processorMarkerFlag[jproc] == 3) {

            unsigned jel = nextElem[jproc];
            unsigned interiorNode = sol->GetMesh()->GetElementDofNumber(jel, 2) - 1;

            unsigned jDof  = sol->GetMesh()->GetSolutionDof(interiorNode, jel, 2);    // global to global mapping between coordinates node and coordinate dof
            double distance2 = 0;

            for (unsigned k = 0; k < _dim; k++) {
              double dk = GetCoordinates(sol, k, jDof, s) - _x[k];   // global extraction and local storage for the element coordinates
              distance2 += dk * dk;
            }

            double modulusKel = sqrt(distance2);

            if (modulusKel < modulus) {
              iel = jel;
              previousElem[_iproc] = previousElem[jproc];
              modulus = modulusKel;

            }
          }

          //END SMART search
        }

        if (iel != sol->GetMesh()->GetElementOffset(_iproc + 1)) {
          //   std::cout << "start element= " << iel << std::endl;
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

    MPI_Bcast(& _elem, 1, MPI_UNSIGNED, _mproc, PETSC_COMM_WORLD);
    //  std::cout << "The marker belongs to process " << _mproc << " and is in element " << _elem  << std::endl;

  }


  void Marker::GetElementSerial(unsigned &previousElem, const Solution* sol, const double &s)
  {

    //std::cout << " SERIALE " << std::endl << std::flush;

    unsigned currentElem = _elem;
    bool elementHasBeenFound = false;
    bool pointIsOutsideThisProcess = false;
    bool pointIsOutsideTheDomain = false;
    unsigned counter = 0;
    std::vector <unsigned> searchHistory(counter + 1);
    searchHistory[counter] = previousElem; //previousElem is always searchHistory[searchHistory.size()-1]

    //BEGIN next element search
    while (elementHasBeenFound + pointIsOutsideThisProcess + pointIsOutsideTheDomain == 0) {
      //std::cout << previousElem << " " << currentElem << std::endl << std::flush;
      if (_dim == 2) {
        _elem = GetNextElement2D(currentElem, previousElem, sol, s);
      }
      else if (_dim == 3) {
        _elem = GetNextElement3D(currentElem, searchHistory, sol, s);

      }
      counter++;
      searchHistory.resize(counter + 1);
      searchHistory[counter] = currentElem;

      //previousElem = currentElem; commented after adding searchHistory

      // std::cout << previousElem << " " << currentElem << std::endl << std::flush;

      // std::cout << previousElem << std::cout << std::flush;

      //std::cout << "previousElem = " << searchHistory[counter] << " " << " nextElem = " << _elem <<  std::endl << std::flush;

      if (_elem == currentElem) {
        elementHasBeenFound = true;
      }
      else if (_elem  == UINT_MAX) {
        pointIsOutsideTheDomain = true;
        break;
      }
      else {
        _mproc = sol->GetMesh()->BisectionSearch_find_processor_of_dof(_elem , 3);
        if (_mproc != _iproc) {
          pointIsOutsideThisProcess = true;
          break;
        }
        else {
          currentElem = _elem;
        }
      }
    }
    //previousElem = _elem;

//     if(elementHasBeenFound) {
//       //  std::cout << " The marker belongs to element " << _elem << std::endl;
//     }
//     else if(pointIsOutsideTheDomain) {
//       //  std::cout << " The marker does not belong to this domain" << std::endl;
//     }
//     else if(pointIsOutsideThisProcess) {
//       //  std::cout << "proc " << _iproc << " believes the marker is in proc = " << _mproc << std::endl;
//     }

    // std::cout << "FINE SERIALE " << std::endl;

    //END next element search

    std::vector < unsigned > ().swap(searchHistory);

  }


  int Marker::FastForward(const unsigned &iel, const unsigned &previousElem, const Solution* sol, const double &s)
  {

    short unsigned linear = 0;

    unsigned nDofs = sol->GetMesh()->GetElementDofNumber(iel, linear);
    short unsigned ielType = sol->GetMesh()->GetElementType(iel);

    //BEGIN extraction nodal coordinate values
    std::vector< std::vector < double > > xv(_dim);

    for (unsigned k = 0; k < _dim; k++) {
      xv[k].resize(nDofs);
    }

    for (unsigned i = 0; i < nDofs; i++) {
      unsigned iDof  = sol->GetMesh()->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < _dim; k++) {
        xv[k][i] = GetCoordinates(sol, k, iDof, s) - _x[k]; // global extraction and local storage for the element coordinates
      }
    }
    //END extraction


    //BEGIN projection nodal to polynomial coefficients
    std::vector < std::vector < double > > a;
    ProjectNodalToPolynomialCoefficients(a, xv, ielType, linear);
    //END projection

    //BEGIN inverse mapping search
    std::vector < double > phi;
    std::vector < std::vector < double > > gradPhi;

    std::vector < double > xi(_dim);
    for (int k = 0; k < _dim; k++) {
      xi[k] = _localCentralNode[ielType][k];
    }

    GetPolynomialShapeFunctionGradient(phi, gradPhi, xi, ielType, linear);

    std::vector < double > v(_dim, 0.);
    std::vector < std::vector < double > > J(_dim);

    for (int k = 0; k < _dim; k++) {
      J[k].assign(_dim, 0.);
    }

    for (int k = 0; k < _dim; k++) {
      for (int i = 0; i < nDofs; i++) {
        v[k] -= a[k][i] * phi[i];

        for (int i1 = 0; i1 < _dim; i1++) {
          J[k][i1] += a[k][i] * gradPhi[i][i1];
        }
      }
    }

    std::vector < std::vector < double > >  Jm1;
    InverseMatrix(J, Jm1);

    std::vector < double > vt(_dim, 0.);
    for (unsigned i = 0; i < _dim; i++) {
      for (unsigned j = 0; j < _dim; j++) {
        vt[i] += Jm1[i][j] * v[j];
      }
    }

    double maxProjection = 0.;
    unsigned faceIndex = 0;
    int nextElem;
    for (unsigned jface = 0; jface < faceNumber[ielType]; jface++) {
      double projection = 0.;
      for (unsigned k = 0; k < _dim; k++) {
        projection += vt[k] * faceNormal[ielType][jface][k];
      }
      int jelement = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(iel, jface) - 1);
      if (projection > maxProjection && jelement != previousElem) {
        maxProjection = projection;
        faceIndex = jface;
        nextElem = jelement;
      }
    }
    return nextElem;

  }



  unsigned Marker::GetNextElement2D(const unsigned &currentElem, const unsigned &previousElem, const Solution* sol, const double &s)
  {


    int nextElem ;
    bool markerIsInElement = false;
    bool nextElementFound = false;
    short unsigned currentElementType = sol->GetMesh()->GetElementType(currentElem);
    double epsilon  = 10.e-10;
    double epsilon2  = epsilon * epsilon;
    double t;

    std::vector<double> xc(_dim, 0); //stores the coordinates of the face node of currentElem
    unsigned faceNodeLocalIndex = (currentElementType == 3) ? 8 : 6;
    unsigned faceNodeDof = sol->GetMesh()->GetSolutionDof(faceNodeLocalIndex, currentElem, 2);
    for (unsigned k = 0; k < _dim; k++) {
      xc[k] = GetCoordinates(sol, k, faceNodeDof, s) - _x[k]; // coordinates are translated so that the marker is the new origin
    }

    if (xc[0]*xc[0] < epsilon2 && xc[1]*xc[1] < epsilon2) {
      //   std::cout << "the marker is the central face node" << std::endl;
      markerIsInElement = true; //the marker is xc
    }
    else {
      unsigned faceNodeNumber = facePointNumber[currentElementType][_solType];
      std::vector< std::vector < double > > xv(_dim);   //stores the coordinates of the vertices and midpoints of the element, the first and the last are the same
      for (unsigned k = 0; k < _dim; k++) {
        xv[k].reserve(faceNodeNumber);
        xv[k].resize(faceNodeNumber - 1);
      }
      for (unsigned i = 0; i < faceNodeNumber - 1; i++) {
        unsigned inodeDof  = sol->GetMesh()->GetSolutionDof(facePoints[currentElementType][_solType][i], currentElem, 2);
        for (unsigned k = 0; k < _dim; k++) {
          xv[k][i] = GetCoordinates(sol, k, inodeDof, s) - _x[k];
        }
      }

      std::vector < double >  xcs;
      double radius;
      GetConvexHullSphere(xv, xcs, radius, 0.1);
      std::vector< std::vector < double > > xe;
      GetBoundingBox(xv, xe, 0.1);

      double radius2 = radius * radius;
      double d2 = 0.;
      for (int d = 0; d < _dim; d++) {
        d2 += xcs[d] * xcs[d];
      }
      bool insideHull = true;
      if (d2 > radius2) {
        insideHull = false;
      }
      for (unsigned k = 0; k < _dim; k++) {
        if (xe[k][0] * xe[k][1] > 0.) {
          insideHull = false;
        }
      }
//       if(!insideHull) {
//         nextElem = FastForward(currentElem, previousElem);
//         nextElementFound = true;
//       }
//
      if (!insideHull) {
        nextElem = FastForward(currentElem, previousElem, sol, s);
        if (nextElem >= 0) {
          nextElementFound = true;
        }
        else {
          insideHull = true;
        }
      }
      if (insideHull) {
        std::vector<double> r(_dim, 0);   //coordinates of the intersection point between the line of the edges and the line that connects the marker and the face node
        for (unsigned k = 0; k < _dim; k++) {
          xv[k].resize(faceNodeNumber);
          xv[k][faceNodeNumber - 1] = xv[k][0];
        }

        //if(true) {
        //BEGIN look for face intersection

        // rescaling coordinates to properly handle different scales of meshes

        double length = 0.;
        double sum = 0.;

        for (unsigned i = 0; i < faceNodeNumber - 1; i++) {
          for (unsigned k = 0; k < _dim; k++) {
            sum += (xv[k][i + 1] - xv[k][i]) * (xv[k][i + 1] - xv[k][i]);
          }

          length += sqrt(sum);
        }

        length /= faceNodeNumber;
        // std::cout << "length= " << length << std::endl;

        for (unsigned k = 0; k < _dim; k++) {
          xc[k] /= length;
          for (unsigned i = 0; i < faceNodeNumber; i++) {
            xv[k][i] /= length;
          }
        }

        for (unsigned i = 0 ; i < faceNodeNumber - 1; i++) {

          // let's find the plane passing through the points xv[][i], xv[][i+1] and xv[][2] = xv[][i] but with z = length .
          double A = (xv[1][i + 1] - xv[1][i]);
          double B = -(xv[0][i + 1] - xv[0][i]);

          //std::cout << "A= " << A << " , " <<"B= " << B <<std::endl;

          double tBottom = (A * xc[0] + B * xc[1]) ;
          double tTop = A * xv[0][i] + B * xv[1][i];
          //std::cout << "tBottom = " << tBottom << " , " << "A= " << A << " , " <<  "B= " << B << " , " << "xv[1][" << i << "] =" << xv[1][i] << " , " <<  "tTop = " <<   tTop << std::endl;

          if (fabs(tBottom) >= epsilon || fabs(tTop) < epsilon) {
            //now let's find the coordinates of the intersection point r
            t = tTop / tBottom ;
            //std::cout << "t = " << t << std::endl;

            for (unsigned k = 0; k < _dim; k++) {
              r[k] = t * xc[k];
              //std::cout << "r[" << k << "] = " << r[k] <<std::endl;
            }

            if (t < 1.) {  //if not, it means the point r is far away from the marker, and we don't want to go in that direction

              std::vector< std::vector < double > > xvr(_dim);

              for (unsigned k = 0; k < _dim; k++) {
                xvr[k].reserve(9);
              }

              for (unsigned k = 0; k < _dim; k++) {
                xvr[k].resize(faceNodeNumber);
              }

              //now we have to determine if r is inside edge i
              for (unsigned j = 0; j < faceNodeNumber; j++) {
                for (unsigned k = 0; k < _dim; k++) {
                  xvr[k][j] = xv[k][j] - r[k];     //transate again the reference frame so that the origin is r
                }
              }


              if ((xvr[0][i] * xvr[0][i]  + xvr[1][i] * xvr[1][i]) < epsilon2 ||
                  (xvr[0][i + 1]*xvr[0][i + 1] + xvr[1][i + 1]*xvr[1][i + 1]) < epsilon2) {
                // std::cout << "intersection on a vertex of the edge" << std::endl;

                if (fabs(t) < epsilon || t < 0) {  //this means the marker is on one of the edges

                  // if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is one of the nodes" << std::endl;

                  // if(t < 0) std::cout << "setting markerIsInElement = true because r is one of the nodes" << std::endl;

                  markerIsInElement = true;
                  break;
                }
                else {
                  unsigned nodeIndex = (_solType == 0) ? i : i / 2;

                  nextElem = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(currentElem, nodeIndex) - 1);
                  if (nextElem != previousElem) {
                    nextElementFound = true;
                  }
                  break;
                }

              }


              else if (xvr[0][i]*xvr[0][i + 1] < 0 || xvr[1][i]*xvr[1][i + 1] < 0) {
                // std::cout << "intersection on an edge" << std::endl;

                if (fabs(t) < epsilon || t < 0) {  //this means the marker is on one of the edges

                  //  if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges " << std::endl;

                  // if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges " << std::endl;

                  markerIsInElement = true;
                  break;
                }
                else {

                  unsigned nodeIndex = (_solType == 0) ? i : i / 2;

                  nextElem = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(currentElem, nodeIndex) - 1);
                  if (nextElem != previousElem) {
                    nextElementFound = true;
                  }
                  break;
                }
              }
            } // closes the " if t < 1 "
          } // closes if
//           else { //questa e' una porcata
//             // std::cout << "The plane of edge " << i << "does not intersect the line" <<std::endl;
//           }
        } //closes the for on the nodes
        //END look for face intersection
      } // closes the else after the sphere check
    }// closes the else before the for

    if (markerIsInElement == true) {
      nextElem = currentElem;
      // std::cout << "The marker belongs to element " << currentElem << std::endl;
    }

    if (nextElementFound == true) {
      // std::cout << "The marker does not belong to element " << currentElem << std::endl;
    }


    // std::cout << "markerIsInElement = " << markerIsInElement << " , " << "nextElementFound= " << nextElementFound << ", " << "nextElem = " << nextElem << std::endl;

    return (nextElem >= 0) ? nextElem : UINT_MAX;

  }


  unsigned Marker::GetNextElement3D(const unsigned & currentElem, const unsigned &previousElem, const Solution* sol, const double &s)
  {

    unsigned nDofs = sol->GetMesh()->GetElementDofNumber(currentElem, _solType);
    unsigned nFaceDofs = (_solType == 2) ? nDofs - 1 : nDofs;
    int nextElem ;
    bool markerIsInElement = false;
    bool nextElementFound = false;
    short unsigned currentElementType = sol->GetMesh()->GetElementType(currentElem);
    double epsilon  = 10.e-10;
    double epsilon2  = epsilon * epsilon;
    double t;

    std::vector<double> xc(_dim, 0); //stores the coordinates of the central node of currentElem
    unsigned centralNodeLocalIndex;

    if (currentElementType == 0) centralNodeLocalIndex = 26;
    else if (currentElementType == 1) centralNodeLocalIndex = 14;
    else if (currentElementType == 2) centralNodeLocalIndex = 20;

    unsigned centralNodeDof = sol->GetMesh()->GetSolutionDof(centralNodeLocalIndex, currentElem, 2);

    for (unsigned k = 0; k < _dim; k++) {
      xc[k] = GetCoordinates(sol, k, centralNodeDof, s)  - _x[k]; // coordinates are translated so that the marker is the new origin
      //std::cout << " xc[" << k << "]= " <<xc[k] <<std::endl;
    }


    //BEGIN ONLY FOR TESTING (TO BE REMOVED ONCE FIXED THE PROBLEM)

//     if(currentElem == 2057){
//       std::cout<<std::endl<<std::flush;
//       for(unsigned iface = 0; iface < sol->GetMesh()->GetElementFaceNumber(currentElem); iface++) {
//         for(unsigned itri = 0; itri < trianglesPerFace[currentElementType][_solType][iface]; itri ++) {
// 	  for(unsigned i = 0; i < 4; i++) {
// 	    unsigned itriDof  = sol->GetMesh()->GetSolutionDof(faceTriangleNodes[currentElementType][_solType][iface][itri][i], currentElem, 2);
// 	    for(unsigned k = 0; k < _dim; k++) {
// 	      std::cout << GetCoordinates(mesh,k,itriDof) <<" ";
// 	    }
// 	    std::cout<<std::endl<<std::flush;
// 	  }
// 	  std::cout<<std::endl<<std::flush;
// 	}
// 	std::cout<<std::endl<<std::flush;
//       }
//
//       unsigned centralNodeDof = sol->GetMesh()->GetSolutionDof(centralNodeLocalIndex, currentElem, 2);
//
//       for(unsigned k = 0; k < _dim; k++) {
// 	std::cout << GetCoordinates(mesh,k,centralNodeDof) <<" ";     // coordinates are translated so that the marker is the new origin
//       }
//       std::cout<<std::endl<<std::flush;
//       for(unsigned k = 0; k < _dim; k++) {
// 	std::cout << _x[k] << " ";    // coordinates are translated so that the marker is the new origin
//       }
//       std::cout<<std::endl<<std::flush;
//       exit(1);
//     }

    //END ONLY FOR TESTING (TO BE REMOVED ONCE FIXED THE PROBLEM)



    if (xc[0]*xc[0] < epsilon2 && xc[1]*xc[1] < epsilon2 && xc[2]*xc[2] < epsilon2) {
      //   std::cout << "the marker is the central element node" << std::endl;
      markerIsInElement = true; //the marker is xc
    }

    else {

      //BEGIN Setting up the fast FastForward

      std::vector< std::vector < double > > xvv(_dim);

      for (unsigned k = 0; k < _dim; k++) {
        xvv[k].reserve(nFaceDofs);
      }

      for (unsigned k = 0; k < _dim; k++) {
        xvv[k].resize(nFaceDofs);
      }

      for (unsigned i = 0; i < nFaceDofs; i++) {
        unsigned inodeDof  = sol->GetMesh()->GetSolutionDof(i, currentElem, _solType);

        for (unsigned k = 0; k < _dim; k++) {
          xvv[k][i] = GetCoordinates(sol, k, inodeDof, s) - _x[k];
        }
      }

      std::vector < double >  xcs;
      double radius;
      GetConvexHullSphere(xvv, xcs, radius, 0.5);
      std::vector< std::vector < double > > xe;
      GetBoundingBox(xvv, xe, 0.5);

      double radius2 = radius * radius;
      double d2 = 0.;
      for (int d = 0; d < _dim; d++) {
        d2 += xcs[d] * xcs[d];
      }
      bool insideHull = true;
      if (d2 > radius2) {
        insideHull = false;
      }
      for (unsigned k = 0; k < _dim; k++) {
        if (xe[k][0] * xe[k][1] > 0.) {
          insideHull = false;
        }
      }

      if (!insideHull) {
        nextElem = FastForward(currentElem, previousElem, sol, s);
        if (nextElem >= 0) {
          nextElementFound = true;
        }
        else {
          insideHull = true;
        }
      }

      if (insideHull) {

        // 	std::cout<<"I am in...."<<std::flush;

        for (unsigned iface = 0; iface < sol->GetMesh()->GetElementFaceNumber(currentElem); iface++) {

          // std::cout << "iface = " << iface << std::endl;

          for (unsigned itri = 0; itri < trianglesPerFace[currentElementType][_solType][iface]; itri ++) {

            //  std::cout << "itri = " << itri << std::endl;

            std::vector<double> xcc(_dim, 0); // will store the coordinates of the center scaled

            unsigned scalarCount = 0;
            std::vector<double> r(_dim, 0);   //coordinates of the intersection point between the plane of itri and the line through the element center point and the marker
            std::vector< std::vector < double > > xv(_dim);   //stores the coordinates of the nodes of the triangle itri

            // fill in the coordinates of the vertices of itri
            for (unsigned k = 0; k < _dim; k++) {
              xv[k].reserve(4);
            }
            for (unsigned k = 0; k < _dim; k++) {
              xv[k].resize(4);
            }

            for (unsigned i = 0; i < 4; i++) {
              unsigned itriDof  = sol->GetMesh()->GetSolutionDof(faceTriangleNodes[currentElementType][_solType][iface][itri][i], currentElem, 2);
              // std::cout << "itriDof = " << itriDof << std::endl;
              for (unsigned k = 0; k < _dim; k++) {
                xv[k][i] = GetCoordinates(sol, k, itriDof, s)  - _x[k];  // coordinates are translated so that the marker is the new origin
              }
            }

            // rescaling coordinates to properly handle different scales of meshes
            double length = 0.;
            double sum = 0.;
            for (unsigned i = 0; i < 3; i++) {
              for (unsigned k = 0; k < _dim; k++) {
                sum += (xv[k][i + 1] - xv[k][i]) * (xv[k][i + 1] - xv[k][i]);
              }
              length += sqrt(sum);
            }

            length /= 4;

            for (unsigned k = 0; k < _dim; k++) {
              xcc[k] = xc[k] / length;
              //std::cout << " xcc[" << k << "]= " <<xcc[k] <<std::endl;
              for (unsigned i = 0; i < 4; i++) {
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

            if (fabs(tBottom) < epsilon && fabs(tTop) >= epsilon) {
              // std::cout << "The plane of face" << itri << "does not intersect the line" <<std::endl;
              break; // must exit the loop on itri
            }

            else { //now let's find the coordinates of the intersection point r
              t = tTop / tBottom ;
              //std::cout << "t = " << t << std::endl;

              for (unsigned k = 0; k < _dim; k++) {
                r[k] = t * xcc[k];
                // std::cout << "r[" << k << "] = " << r[k] <<std::endl;
              }

              if (t < 1) {  //if not, it means the point r is far away from the marker, and we don't want to go in that direction

                for (unsigned i = 0; i < 4; i++) { //now we have to determine if r is inside itri
                  for (unsigned k = 0; k < _dim; k++) {
                    xv[k][i] = xv[k][i] - r[k];     //translate again the reference frame so that the origin is r
                  }
                }

                for (unsigned i = 0; i < 3; i++) {
                  double q0 = xv[1][i] * (xv[2][i] - xv[2][i + 1]) + xv[2][i] * (xv[1][i + 1] - xv[1][i]);
                  double q1 = xv[2][i] * (xv[0][i] - xv[0][i + 1]) + xv[0][i] * (xv[2][i + 1] - xv[2][i]);
                  double q2 = xv[0][i] * (xv[1][i] - xv[1][i + 1]) + xv[1][i] * (xv[0][i + 1] - xv[0][i]);

                  // std::cout << "q0 = " << q0 << " , " << "q1 = " << q1 << " , " << " q2 = " << q2 <<  std::endl;

                  double  scalarProduct = q0 * A + q1 * B + q2 * C;

                  //   std::cout << "fabs(scalarProduct) = " << fabs(scalarProduct) << std::endl;

                  if (scalarProduct > epsilon) {
                    //   std::cout << "r is outside triangle " << itri <<  std::endl;
                    break;

                  }
                  else if (fabs(scalarProduct) < epsilon) {   //scalarProduct == 0

                    if ((xv[0][i] * xv[0][i]  + xv[1][i] * xv[1][i] + xv[2][i] * xv[2][i]) < epsilon2 ||
                        (xv[0][i + 1]*xv[0][i + 1] + xv[1][i + 1]*xv[1][i + 1] + xv[2][i + 1]*xv[2][i + 1]) < epsilon2) {
                      //    std::cout << "intersection on a vertex of itri" << std::endl;
                      if (fabs(t) < epsilon || t < 0) {  //this means the marker is on one of the faces

                        //     if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is one vertex of triangle " << itri << std::endl;
                        //     if(t < 0) std::cout << "setting markerIsInElement = true because r is one vertex of triangle " << itri << std::endl;

                        markerIsInElement = true;
                        break;
                      }
                      else {
                        //     std::cout << "r is in triangle " << itri << std::endl;
                        nextElem = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(currentElem, iface) - 1);
                        if (nextElem != previousElem) {
                          nextElementFound = true;
                        }
                        break;
                      }

                    }


                    else if (xv[0][i]*xv[0][i + 1] < 0 || xv[1][i]*xv[1][i + 1] < 0 || xv[2][i]*xv[2][i + 1] < 0) {
                      //   std::cout << "intersection on an edge of itri" << std::endl;
                      if (fabs(t) < epsilon || t < 0) {  //this means the marker is on one of the faces

                        //    if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges of triangle " << itri << std::endl;
                        //    if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges of triangle " << itri << std::endl;

                        markerIsInElement = true;
                        break;
                      }
                      else {
                        //     std::cout << "r is in triangle " << itri << std::endl;
                        nextElem = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(currentElem, iface) - 1);
                        if (nextElem != previousElem) {
                          nextElementFound = true;
                        }
                        break;
                      }
                    }
                  }
                  else if (scalarProduct < 0) {
                    //    std::cout << " scalarProduct = " << scalarProduct << std::endl;
                    scalarCount++;
                  }
                } // closes the for loop
              } // closes " if t < 1 "
            } // closes the "else" on tBottom = 0


            if (scalarCount == 3) {
              if (fabs(t) < epsilon || t < 0) {  //this means the marker is on one of the faces

                //   if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges of triangle " << itri << std::endl;
                //  if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges of triangle " << itri << std::endl;

                markerIsInElement = true;
                break;
              }
              else {
                //    std::cout << "r is in triangle " << itri << std::endl;
                nextElem = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(currentElem, iface) - 1);
                if (nextElem != previousElem) {
                  nextElementFound = true;
                }
                break;
              }
            }

            if (nextElementFound == true) {
              break;
            }

            if (markerIsInElement == true) {
              break;
            }
          } //end for on itri

          if (nextElementFound == true) {
            break;
          }

          if (markerIsInElement == true) {
            break;
          }
        } //end for on iface
      } //end of else before for on iface
    } // end of the first else

    if (markerIsInElement == true) {
      nextElem = currentElem;
      //  std::cout << "The marker belongs to element " << currentElem << std::endl;
    }
    if (nextElementFound == true) {
      //   std::cout << "The marker does not belong to element " << currentElem << std::endl;
    }
    //  std::cout << "markerIsInElement = " << markerIsInElement << " , " << "nextElementFound= " << nextElementFound << ", " << "nextElem = " << nextElem << std::endl;


    return (nextElem >= 0) ? nextElem : UINT_MAX;

  }



  unsigned Marker::GetNextElement3D(const unsigned & currentElem, const std::vector <unsigned> &searchHistory, const Solution* sol, const double &s)
  {

    unsigned nDofs = sol->GetMesh()->GetElementDofNumber(currentElem, _solType);
    unsigned nFaceDofs = (_solType == 2) ? nDofs - 1 : nDofs;
    int nextElem ;
    bool markerIsInElement = false;
    bool nextElementFound = false;
    short unsigned currentElementType = sol->GetMesh()->GetElementType(currentElem);
    double epsilon  = 10.e-10;
    double epsilon2  = epsilon * epsilon;
    double t;

    std::vector<double> xc(_dim, 0); //stores the coordinates of the central node of currentElem
    unsigned centralNodeLocalIndex;

    if (currentElementType == 0) centralNodeLocalIndex = 26;
    else if (currentElementType == 1) centralNodeLocalIndex = 14;
    else if (currentElementType == 2) centralNodeLocalIndex = 20;

    unsigned centralNodeDof = sol->GetMesh()->GetSolutionDof(centralNodeLocalIndex, currentElem, 2);

    for (unsigned k = 0; k < _dim; k++) {
      xc[k] = GetCoordinates(sol, k, centralNodeDof, s)  - _x[k]; // coordinates are translated so that the marker is the new origin
      //std::cout << " xc[" << k << "]= " <<xc[k] <<std::endl;
    }


    //BEGIN ONLY FOR TESTING (TO BE REMOVED ONCE FIXED THE PROBLEM)

//     if(currentElem == 2057){
//       std::cout<<std::endl<<std::flush;
//       for(unsigned iface = 0; iface < sol->GetMesh()->GetElementFaceNumber(currentElem); iface++) {
//         for(unsigned itri = 0; itri < trianglesPerFace[currentElementType][_solType][iface]; itri ++) {
// 	  for(unsigned i = 0; i < 4; i++) {
// 	    unsigned itriDof  = sol->GetMesh()->GetSolutionDof(faceTriangleNodes[currentElementType][_solType][iface][itri][i], currentElem, 2);
// 	    for(unsigned k = 0; k < _dim; k++) {
// 	      std::cout << GetCoordinates(mesh,k,itriDof) <<" ";
// 	    }
// 	    std::cout<<std::endl<<std::flush;
// 	  }
// 	  std::cout<<std::endl<<std::flush;
// 	}
// 	std::cout<<std::endl<<std::flush;
//       }
//
//       unsigned centralNodeDof = sol->GetMesh()->GetSolutionDof(centralNodeLocalIndex, currentElem, 2);
//
//       for(unsigned k = 0; k < _dim; k++) {
// 	std::cout << GetCoordinates(mesh,k,centralNodeDof) <<" ";     // coordinates are translated so that the marker is the new origin
//       }
//       std::cout<<std::endl<<std::flush;
//       for(unsigned k = 0; k < _dim; k++) {
// 	std::cout << _x[k] << " ";    // coordinates are translated so that the marker is the new origin
//       }
//       std::cout<<std::endl<<std::flush;
//       exit(1);
//     }

    //END ONLY FOR TESTING (TO BE REMOVED ONCE FIXED THE PROBLEM)



    if (xc[0]*xc[0] < epsilon2 && xc[1]*xc[1] < epsilon2 && xc[2]*xc[2] < epsilon2) {
      //   std::cout << "the marker is the central element node" << std::endl;
      markerIsInElement = true; //the marker is xc
    }

    else {

      //BEGIN Setting up the fast FastForward

      std::vector< std::vector < double > > xvv(_dim);

      for (unsigned k = 0; k < _dim; k++) {
        xvv[k].reserve(nFaceDofs);
      }

      for (unsigned k = 0; k < _dim; k++) {
        xvv[k].resize(nFaceDofs);
      }

      for (unsigned i = 0; i < nFaceDofs; i++) {
        unsigned inodeDof  = sol->GetMesh()->GetSolutionDof(i, currentElem, _solType);

        for (unsigned k = 0; k < _dim; k++) {
          xvv[k][i] = GetCoordinates(sol, k, inodeDof, s) - _x[k];
        }
      }

      std::vector < double >  xcs;
      double radius;
      GetConvexHullSphere(xvv, xcs, radius, 0.5);
      std::vector< std::vector < double > > xe;
      GetBoundingBox(xvv, xe, 0.5);

      double radius2 = radius * radius;
      double d2 = 0.;
      for (int d = 0; d < _dim; d++) {
        d2 += xcs[d] * xcs[d];
      }
      bool insideHull = true;
      if (d2 > radius2) {
        insideHull = false;
      }
      for (unsigned k = 0; k < _dim; k++) {
        if (xe[k][0] * xe[k][1] > 0.) {
          insideHull = false;
        }
      }
      if (!insideHull) {
        nextElem = FastForward(currentElem, searchHistory[searchHistory.size() - 1], sol, s);
        if (nextElem >= 0) {
          nextElementFound = true;
        }
        else {
          insideHull = true;
        }
      }

      if (insideHull) {

// 	std::cout<<"I am in...."<<std::flush;

        for (unsigned iface = 0; iface < sol->GetMesh()->GetElementFaceNumber(currentElem); iface++) {

          // std::cout << "iface = " << iface << std::endl;

          for (unsigned itri = 0; itri < trianglesPerFace[currentElementType][_solType][iface]; itri ++) {

            //  std::cout << "itri = " << itri << std::endl;

            std::vector<double> xcc(_dim, 0); // will store the coordinates of the center scaled

            unsigned scalarCount = 0;
            std::vector<double> r(_dim, 0);   //coordinates of the intersection point between the plane of itri and the line through the element center point and the marker
            std::vector< std::vector < double > > xv(_dim);   //stores the coordinates of the nodes of the triangle itri

            // fill in the coordinates of the vertices of itri
            for (unsigned k = 0; k < _dim; k++) {
              xv[k].reserve(4);
            }
            for (unsigned k = 0; k < _dim; k++) {
              xv[k].resize(4);
            }

            for (unsigned i = 0; i < 4; i++) {
              unsigned itriDof  = sol->GetMesh()->GetSolutionDof(faceTriangleNodes[currentElementType][_solType][iface][itri][i], currentElem, 2);
              // std::cout << "itriDof = " << itriDof << std::endl;
              for (unsigned k = 0; k < _dim; k++) {
                xv[k][i] = GetCoordinates(sol, k, itriDof, s)  - _x[k];  // coordinates are translated so that the marker is the new origin
              }
            }

            // rescaling coordinates to properly handle different scales of meshes
            double length = 0.;
            double sum = 0.;
            for (unsigned i = 0; i < 3; i++) {
              for (unsigned k = 0; k < _dim; k++) {
                sum += (xv[k][i + 1] - xv[k][i]) * (xv[k][i + 1] - xv[k][i]);
              }
              length += sqrt(sum);
            }

            length /= 4;

            for (unsigned k = 0; k < _dim; k++) {
              xcc[k] = xc[k] / length;
              //std::cout << " xcc[" << k << "]= " <<xcc[k] <<std::endl;
              for (unsigned i = 0; i < 4; i++) {
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

            if (fabs(tBottom) < epsilon && fabs(tTop) >= epsilon) {
              // std::cout << "The plane of face" << itri << "does not intersect the line" <<std::endl;
              break; // must exit the loop on itri
            }

            else { //now let's find the coordinates of the intersection point r
              t = tTop / tBottom ;
              //std::cout << "t = " << t << std::endl;

              for (unsigned k = 0; k < _dim; k++) {
                r[k] = t * xcc[k];
                // std::cout << "r[" << k << "] = " << r[k] <<std::endl;
              }

              if (t < 1) {  //if not, it means the point r is far away from the marker, and we don't want to go in that direction

                for (unsigned i = 0; i < 4; i++) { //now we have to determine if r is inside itri
                  for (unsigned k = 0; k < _dim; k++) {
                    xv[k][i] = xv[k][i] - r[k];     //translate again the reference frame so that the origin is r
                  }
                }

                for (unsigned i = 0; i < 3; i++) {
                  double q0 = xv[1][i] * (xv[2][i] - xv[2][i + 1]) + xv[2][i] * (xv[1][i + 1] - xv[1][i]);
                  double q1 = xv[2][i] * (xv[0][i] - xv[0][i + 1]) + xv[0][i] * (xv[2][i + 1] - xv[2][i]);
                  double q2 = xv[0][i] * (xv[1][i] - xv[1][i + 1]) + xv[1][i] * (xv[0][i + 1] - xv[0][i]);

                  // std::cout << "q0 = " << q0 << " , " << "q1 = " << q1 << " , " << " q2 = " << q2 <<  std::endl;

                  double  scalarProduct = q0 * A + q1 * B + q2 * C;

                  //   std::cout << "fabs(scalarProduct) = " << fabs(scalarProduct) << std::endl;

                  if (scalarProduct > epsilon) {
                    //   std::cout << "r is outside triangle " << itri <<  std::endl;
                    break;

                  }
                  else if (fabs(scalarProduct) < epsilon) {   //scalarProduct == 0

                    if ((xv[0][i] * xv[0][i]  + xv[1][i] * xv[1][i] + xv[2][i] * xv[2][i]) < epsilon2 ||
                        (xv[0][i + 1]*xv[0][i + 1] + xv[1][i + 1]*xv[1][i + 1] + xv[2][i + 1]*xv[2][i + 1]) < epsilon2) {
                      //    std::cout << "intersection on a vertex of itri" << std::endl;
                      if (fabs(t) < epsilon || t < 0) {  //this means the marker is on one of the faces

                        //     if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is one vertex of triangle " << itri << std::endl;
                        //     if(t < 0) std::cout << "setting markerIsInElement = true because r is one vertex of triangle " << itri << std::endl;

                        markerIsInElement = true;
                        break;
                      }
                      else {
                        //     std::cout << "r is in triangle " << itri << std::endl;
                        nextElem = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(currentElem, iface) - 1);
                        unsigned historyCounter = 0;
                        for(unsigned i = 0; i < searchHistory.size(); i++) {
                          // std::cout << "searchHistory[" << i<< "]"<< searchHistory[i] <<std::endl;
                          if(nextElem != searchHistory[i]) {
                            historyCounter++;
                          }
                        }
                        if(historyCounter == searchHistory.size()) {
                          nextElementFound = true;
                        }
                        break;
                      }

                    }


                    else if (xv[0][i]*xv[0][i + 1] < 0 || xv[1][i]*xv[1][i + 1] < 0 || xv[2][i]*xv[2][i + 1] < 0) {
                      //   std::cout << "intersection on an edge of itri" << std::endl;
                      if (fabs(t) < epsilon || t < 0) {  //this means the marker is on one of the faces

                        //    if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges of triangle " << itri << std::endl;
                        //    if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges of triangle " << itri << std::endl;

                        markerIsInElement = true;
                        break;
                      }
                      else {
                        //     std::cout << "r is in triangle " << itri << std::endl;
                        nextElem = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(currentElem, iface) - 1);
                        unsigned historyCounter = 0;
                        for(unsigned i = 0; i < searchHistory.size(); i++) {
                          //  std::cout << "searchHistory[" << i<< "]"<< searchHistory[i] <<std::endl;
                          if(nextElem != searchHistory[i]) {
                            historyCounter++;
                          }
                        }
                        if(historyCounter == searchHistory.size()) {
                          nextElementFound = true;
                        }
                        break;
                      }
                    }
                  }
                  else if (scalarProduct < 0) {
                    //    std::cout << " scalarProduct = " << scalarProduct << std::endl;
                    scalarCount++;
                  }
                } // closes the for loop
              } // closes " if t < 1 "
            } // closes the "else" on tBottom = 0


            if (scalarCount == 3) {
              if (fabs(t) < epsilon || t < 0) {  //this means the marker is on one of the faces

                //   if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges of triangle " << itri << std::endl;
                //  if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges of triangle " << itri << std::endl;

                markerIsInElement = true;
                break;
              }
              else {
                //    std::cout << "r is in triangle " << itri << std::endl;
                nextElem = (sol->GetMesh()->GetMeshElements()->GetFaceElementIndex(currentElem, iface) - 1);
                unsigned historyCounter = 0;
                for(unsigned i = 0; i < searchHistory.size(); i++) {
                  //  std::cout << "searchHistory[" << i<< "]"<< searchHistory[i] <<std::endl;
                  if(nextElem != searchHistory[i]) {
                    historyCounter++;
                  }
                }
                if(historyCounter == searchHistory.size()) {
                  nextElementFound = true;
                }
                break;
              }
            }

            if (nextElementFound == true) {
              break;
            }

            if (markerIsInElement == true) {
              break;
            }
          } //end for on itri

          if (nextElementFound == true) {
            break;
          }

          if (markerIsInElement == true) {
            break;
          }
        } //end for on iface
      } //end of else before for on iface
    } // end of the first else

    if (markerIsInElement == true) {
      nextElem = currentElem;
      //  std::cout << "The marker belongs to element " << currentElem << std::endl;
    }
    if (nextElementFound == true) {
      //   std::cout << "The marker does not belong to element " << currentElem << std::endl;
    }
    //  std::cout << "markerIsInElement = " << markerIsInElement << " , " << "nextElementFound= " << nextElementFound << ", " << "nextElem = " << nextElem << std::endl;


    return (nextElem >= 0) ? nextElem : UINT_MAX;

  }



  void Marker::InverseMapping(const unsigned & iel, const unsigned & solType,
                              const std::vector< double > &x, std::vector< double > &xi, const Solution* sol, const double &s)
  {

    //std::cout << " ----------------------  Outputs of the inverse mapping ---------------------- " << std::endl;

    unsigned nDofs = sol->GetMesh()->GetElementDofNumber(iel, solType);
    short unsigned ielType = sol->GetMesh()->GetElementType(iel);

    //std::cout << "solType = " << solType << " , " << "nDofs =" <<  nDofs << std::endl;

    //BEGIN extraction nodal coordinate values
    std::vector< std::vector < double > > xv(_dim);

    for (unsigned k = 0; k < _dim; k++) {
      xv[k].resize(nDofs);
    }

    for (unsigned i = 0; i < nDofs; i++) {
      unsigned iDof  = sol->GetMesh()->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < _dim; k++) {
        xv[k][i] = GetCoordinates(sol, k, iDof, s) ;  // global extraction and local storage for the element coordinates
      }
    }
    //END extraction


    //BEGIN projection nodal to polynomial coefficients
    std::vector < std::vector < double > > a;
    ProjectNodalToPolynomialCoefficients(a, xv, ielType, solType);
    //END projection


    //BEGIN inverse mapping search
    bool convergence = false;
    for (unsigned k = 0; k < _dim; k++) {
      //   std::cout << "xi[" << k << "] = " << xi[k] << " ";
    }
    // std::cout << std::endl;

    while (!convergence) {

      std::vector < double > phi;
      std::vector < std::vector < double > > gradPhi;
      std::vector < std::vector < std::vector < double > > > hessPhi;

      GetPolynomialShapeFunctionGradientHessian(phi, gradPhi, hessPhi, xi, ielType, solType);

      if (solType == 0) {
        counter++;
        convergence = GetNewLocalCoordinates(xi, x, phi, gradPhi, a);
      }
      else {
        counter++;
        //convergence = GetNewLocalCoordinates(xi, x, phi, gradPhi, a);
        convergence = GetNewLocalCoordinatesHess(xi, x, phi, gradPhi, hessPhi, a);
      }
      for (unsigned k = 0; k < _dim; k++) {
        //   std::cout << "xi[" << k << "] = " << xi[k] << " ";
      }
      // std::cout << std::endl;
    }
    //END inverse mapping search

    return;


  }

  void Marker::InverseMappingTEST(std::vector < double > &x, const Solution* sol, const double &s)
  {

    for (int solType = 0; solType < NFE_FAMS_C_ZERO_LAGRANGE; solType++) {
      // std::cout << "\n\n--------------------------------------------------" << std::endl;
      // std::cout << "solType = " << solType << std::endl;

      for (int iel = sol->GetMesh()->GetElementOffset(_iproc); iel < sol->GetMesh()->GetElementOffset(_iproc + 1); iel++) {
        //for(int iel = 494; iel < 495; iel++) {
        //  std::cout << "iel = " << iel << std::endl;
        //  std::cout << "--------------------------------------------------\n" << std::endl;

        unsigned nDofs = sol->GetMesh()->GetElementDofNumber(iel, solType);
        short unsigned ielType = sol->GetMesh()->GetElementType(iel);

        std::vector < std::vector < double > > xv(nDofs);

        for (unsigned i = 0; i < nDofs; i++) {
          xv[i].resize(_dim);
          unsigned iDof  = sol->GetMesh()->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof

          for (unsigned k = 0; k < _dim; k++) {
            xv[i][k] = GetCoordinates(sol, k, iDof, s);  // global extraction and local storage for the element coordinates
          }
          //std::cout <<"x"<< i+1 <<"="<< xv[i][2]<<"; " << std::endl;

        }

        //exit(0);
        for (int j = 0; j < nDofs; j++) {
          //   std::cout << j << std::endl;
          std::vector < double > xiT(_dim);

          for (unsigned k = 0; k < _dim; k++) {
            xiT[k] = *(sol->GetMesh()->_finiteElement[ielType][0]->GetBasis()->GetXcoarse(j) + k);
          }

          // This is the test


          for (int k = 0; k < _dim; k++) {
            //   std::cout << "xv[" << k << "]= " << xv[j][k] <<  " ";
          }

          //  std::cout << std::endl;


          std::vector < double > xi;
          GetClosestPointInReferenceElement(xv, xv[j], ielType, xi);


// 	  std::vector < double > xi(_dim);
//           for(int k = 0; k < _dim; k++) {
//              xi[k] = _localCentralNode[ielType][k];
//           }

          for (int itype = 0; itype <= solType; itype++) {
            InverseMapping(iel, itype, xv[j], xi, sol, s);
            //  std::cout << std::endl;
          }

// 	  InverseMapping(iel, 0, xv[j], xi);
// 	  if(solType > 0){
// 	    InverseMapping(iel, solType, xv[j], xi);
// 	    std::cout << std::endl;
// 	  }


          bool test = true;
          for (int k = 0; k < _dim; k++) {
            //    std::cout << "xiT[" << k << "]= " << xiT[k] <<  " xi[" << k << "]= " << xi[k];
            //    std::cout << " error: " << xiT[k] - xi[k] << std::endl;
            if (fabs(xiT[k] - xi[k]) > 1.0e-3) {
              test = false;
            }
          }
          if (test == false) {
            //   std::cout << "Inverse map test failed " << std::endl;
            abort();
          }
          //   std::cout << "--------------------------------------------------\n" << std::endl;
        }
      }
    }

    // std::cout << "total number of calls= " << counter << std::endl;
  }


  //this function returns the position of the marker at time T given the position at time T0 = 0, given the function f and the stepsize h
  void Marker::Advection(const unsigned &n, const double& T, const Solution* sol)
  {

    double s1 = 0.;

    //BEGIN  Initialize the parameters for all processors, to be used when awake

    std::vector < unsigned > solVIndex(_dim);
    solVIndex[0] = sol->GetIndex("U");    // get the position of "U" in the ml_sol object
    solVIndex[1] = sol->GetIndex("V");    // get the position of "V" in the ml_sol object
    if (_dim == 3) solVIndex[2] = sol->GetIndex("W");      // get the position of "V" in the ml_sol object
    unsigned solVType = sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

    std::vector < double > phi;
    std::vector < std::vector<double > > V(2);
    std::vector < std::vector < std::vector < double > > > aV;
    std::vector < std::vector < std::vector < std::vector < double > > > > aX;
    //END

    //BEGIN Numerical integration scheme
    // determine the step size
    double h = T / n;
    bool integrationIsOver = (_elem != UINT_MAX) ? false : true;

    unsigned order = 1;
    unsigned step = 0.;

    if (_iproc == _mproc) {
      if (_elem != UINT_MAX) {
        FindLocalCoordinates(solVType, aX, true, sol, s1);
      }
      _K.resize(order);
      for (unsigned k = 0; k < order; k++) {
        _K[k].resize(_dim);
      }
    }

    while (integrationIsOver == false) {
      unsigned mprocOld = _mproc;


      // std::cout << "INIZIA UN CICLO DI INTEGRAZIONE " << " " << " _elem da cui si parte = " << _elem << " " << " mprocOld = " << mprocOld << " " << "_mproc = " << _mproc <<  std::endl;


      unsigned previousElem;
      if (_iproc == _mproc) { //single process
        bool pcElemUpdate = true ;
        while (step < n * order) {

          unsigned tstep = step / order;
          unsigned istep = step % order;

          if (istep == 0) {
            _x0 = _x;
            for (unsigned k = 0; k < order; k++) {
              _K[k].assign(_dim, 0.);
            }
          }

          //  std::cout << " -----------------------------" << "step = " <<  step << " tstep = " << tstep << " istep = " << istep << " -----------------------------" << std::endl;
          //  std::cout << " _iproc = " << _iproc << std::endl;

          updateVelocity(V, solVIndex, solVType, aV, phi, pcElemUpdate, sol); //send xk

          double s = (tstep + _c[order - 1][istep]) / n;

          for (unsigned i = 0; i < _dim; i++) {
            _K[istep][i] = (s * V[0][i] + (1. - s) * V[1][i]) * h;
          }

          step++;
          istep++;

          if (istep < order) {

            for (unsigned i = 0; i < _dim; i++) {
              _x[i] = _x0[i];
              for (unsigned j = 0; j < order; j++) {
                _x[i] +=  _a[order - 1][istep][j] * _K[j][i];
              }
            }
          }
          else if (istep == order) {

            for (unsigned i = 0; i < _dim; i++) {
              _x[i] = _x0[i];
              for (unsigned j = 0; j < order; j++) {
                _x[i] += _b[order - 1][j] * _K[j][i];
              }
            }
          }


          //BEGIN TO BE REMOVED
//           for(unsigned i = 0; i < _dim; i++) {
//             //    std::cout << "_x[" << i << "] = " << _x[i] ;
//             //    std::cout << " " ;
//           }
          //  std::cout << std::endl;
          //END TO BE REMOVED

          pcElemUpdate = false;
          unsigned iel = _elem;
          previousElem = _elem;
          GetElementSerial(previousElem, sol, s1);

          if (_elem == UINT_MAX) { //out of the domain
            //    std::cout << " the marker has been advected outside the domain " << std::endl;
            break;
          }
          else if (iel != _elem && _iproc != _mproc) {   //different element different process
            break;
          }
          else if (iel != _elem) {   //different element same process
            pcElemUpdate = true;
            FindLocalCoordinates(solVType, aX, pcElemUpdate, sol, s1);
          }
          else {   //same element same process
            FindLocalCoordinates(solVType, aX, pcElemUpdate, sol, s1);
          }
        }
      }
      //   std::cout << step;

      //  std::cout << "prima del bcast _elem = " << _elem << std::endl;

      // all processes
      MPI_Bcast(& _elem, 1, MPI_UNSIGNED, mprocOld, PETSC_COMM_WORLD);
      MPI_Bcast(& step, 1, MPI_UNSIGNED, mprocOld, PETSC_COMM_WORLD);


      // std::cout << "dopo il bcast _elem = " << _elem << std::endl;

      if (_elem == UINT_MAX) {
        //   std::cout << " the marker has been advected outside the domain " << std::endl;
        break;
      }
      else {
        _mproc = sol->GetMesh()->BisectionSearch_find_processor_of_dof(_elem, 3);

        //  std::cout << "_mproc = " << _mproc << " " << " mprocOld = " << mprocOld << std::endl;

        if (_mproc != mprocOld) {

          //   std::cout << "SON QUA previousElem = " << previousElem << std::endl;

          GetElement(previousElem, mprocOld, sol, s1);  //fissato, WARNING QUESTO e' il problema, dovrebbe dire 724 invece dice 716

          //  std::cout << "_elem = " << _elem << std::endl;
          //  std::cout << "_mproc = " << _mproc << " " << " mprocOld = " << mprocOld << std::endl;

          if (_elem == UINT_MAX) break;

          if (_mproc != mprocOld) {
            if (mprocOld == _iproc) {
              unsigned istep = step % order;
              if (istep != 0) {
                for (int i = 0; i < order; i++) {
                  MPI_Send(&_K[i][0], _dim, MPI_DOUBLE, _mproc, i , PETSC_COMM_WORLD);
                }
                MPI_Send(&_x0[0], _dim, MPI_DOUBLE, _mproc, order , PETSC_COMM_WORLD);
              }
              std::vector < double > ().swap(_xi);
              std::vector < double > ().swap(_x0);
              std::vector < std::vector < double > > ().swap(_K);
              std::vector < std::vector < std::vector < std::vector < double > > > >().swap(aX);
            }
            else if (_mproc == _iproc) {
              _x0.resize(_dim);
              _K.resize(order);
              for (unsigned i = 0; i < order; i++) {
                _K[i].resize(_dim);
              }
              unsigned istep = step % order;
              if (istep != 0) {
                for (int i = 0; i < order; i++) {
                  MPI_Recv(&_K[i][0], _dim, MPI_DOUBLE, mprocOld, i , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                MPI_Recv(&_x0[0], _dim, MPI_DOUBLE, mprocOld, order , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
              }
            }
          }
          if (_mproc == _iproc) { //WARNING this now should be outside and was causing the problem with different processes
            FindLocalCoordinates(solVType, aX, true, sol, s1);
          }
        }
      }
//     std::cout << "step = " << step << std::endl;
      if (step == n * order) {
        integrationIsOver = true;
        //     std::cout << "Integration is over, point in proc " << _mproc << std::endl;
      }
    }
  }


  void Marker::GetElement(unsigned &previousElem, const unsigned &previousMproc, const Solution* sol, const double &s)
  {

    unsigned mprocOld = previousMproc;

    if (mprocOld == _iproc) {
      MPI_Send(&previousElem, 1, MPI_UNSIGNED, _mproc, 1 , PETSC_COMM_WORLD);
      MPI_Send(&_x[0], _dim, MPI_DOUBLE, _mproc, 2 , PETSC_COMM_WORLD);
      std::vector < double > ().swap(_x);
    }
    else if (_mproc == _iproc) {
      MPI_Recv(&previousElem, 1, MPI_UNSIGNED, mprocOld, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      _x.resize(_dim);
      MPI_Recv(&_x[0], _dim, MPI_DOUBLE, mprocOld, 2 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    mprocOld = _mproc;
    while (true) {    //Corretto ora dovrebbe funzionare. WARNING credo debba essere while(!found) perche' senno' non fa mai GetElementSerial e non cambia mai 716 in 724, pero' se si mette !found non ranna piu niente
      if (_mproc == _iproc) {
        GetElementSerial(previousElem, sol, s);
      }
      MPI_Bcast(& _elem, 1, MPI_UNSIGNED, mprocOld, PETSC_COMM_WORLD);
      if (_elem == UINT_MAX) {
        //   std::cout << " the marker has been advected outside the domain " << std::endl;
        break;
      }
      else {
        _mproc = sol->GetMesh()->BisectionSearch_find_processor_of_dof(_elem, 3);
        if (_mproc == mprocOld) {
          break;
        }
        else {
          if (mprocOld == _iproc) {
            MPI_Send(&previousElem, 1, MPI_UNSIGNED, _mproc, 1 , PETSC_COMM_WORLD);
            MPI_Send(&_x[0], _dim, MPI_DOUBLE, _mproc, 2 , PETSC_COMM_WORLD);
            std::vector < double > ().swap(_x);
          }
          else if (_mproc == _iproc) {
            MPI_Recv(&previousElem, 1, MPI_UNSIGNED, mprocOld, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
            _x.resize(_dim);
            MPI_Recv(&_x[0], _dim, MPI_DOUBLE, mprocOld, 2 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          }
          mprocOld = _mproc;
        }
      }
    }
  }



  void Marker::ProjectVelocityCoefficients(const std::vector<unsigned> &solVIndex,
      const unsigned & solVType,  const unsigned & nDofsV,
      const unsigned & ielType, std::vector < std::vector < std::vector < double > > > &a, const Solution* sol)
  {

    bool timeDependent = true;
    for (unsigned  k = 0; k < _dim; k++) {
      if (sol->GetSolutionTimeOrder(solVIndex[k]) != TIME_DEPENDENT) {
        timeDependent = false;
      }
    }

    std::vector < std::vector < double > >  solV(_dim);    // local solution
    std::vector < std::vector < double > >  solVold(_dim);    // local solution

    for (unsigned  k = 0; k < _dim; k++) {
      solV[k].resize(nDofsV);
      if (timeDependent) {
        solVold[k].resize(nDofsV);
      }
    }
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = sol->GetMesh()->GetSolutionDof(i, _elem, solVType);    // global to global mapping between solution node and solution dof
      for (unsigned  k = 0; k < _dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        if (timeDependent) {
          solVold[k][i] = (*sol->_SolOld[solVIndex[k]])(solVDof);
        }
      }
    }

    a.resize(2);

    ProjectNodalToPolynomialCoefficients(a[0], solV, ielType, solVType);
    if (timeDependent) {
      ProjectNodalToPolynomialCoefficients(a[1], solVold, ielType, solVType);
    }
    else {
      a[1] = a[0];
    }
  }


  void Marker::updateVelocity(std::vector< std::vector <double> > & V,
                              const std::vector < unsigned > &solVIndex, const unsigned & solVType,
                              std::vector < std::vector < std::vector < double > > > &a,  std::vector < double > &phi,
                              const bool & pcElemUpdate, const Solution* sol)
  {


    unsigned nDofsV = sol->GetMesh()->GetElementDofNumber(_elem, solVType);
    short unsigned ielType = sol->GetMesh()->GetElementType(_elem);

    if (pcElemUpdate) {
      ProjectVelocityCoefficients(solVIndex, solVType, nDofsV, ielType, a, sol);
    }

    GetPolynomialShapeFunction(phi, _xi, ielType, solVType);

    V.resize(2);
    V[0].assign(_dim, 0.);
    V[1].assign(_dim, 0.);
    for (unsigned i = 0; i < _dim; i++) {
      for (unsigned j = 0; j < nDofsV; j++) {
        V[0][i] += a[0][i][j] * phi[j];
        V[1][i] += a[1][i][j] * phi[j];
      }
    }


  }

  void Marker::FindLocalCoordinates(const unsigned & solType, std::vector < std::vector < std::vector < std::vector < double > > > > &aX, const bool & pcElemUpdate, const Solution* sol, const double &s)
  {

    //BEGIN TO BE REMOVED
//     std::cout << "ENTRIAMO NELL' INVERSE MAPPING, _elem =" << _elem << std::endl;
//     for(unsigned i = 0; i < 3; i++) {
//      std::cout << "_x[" << i << "]= " << _x[i] << std::endl;
//     }
    //END TO BE REMOVED

    short unsigned elemType = sol->GetMesh()->GetElementType(_elem);
    unsigned nDofs = sol->GetMesh()->GetElementDofNumber(_elem, solType);

    if (pcElemUpdate) {

      //BEGIN extraction nodal coordinate values
      std::vector <std::vector< std::vector < double > > >xv;



      if (!sol->GetIfFSI()) {
        aX.resize(1);
        xv.resize(1);
      }
      else {
        aX.resize(2);
        xv.resize(2);
      }

      for (unsigned i1 = 0; i1 < aX.size(); i1++) {
        xv[i1].resize(_dim);
        for (unsigned k = 0; k < _dim; k++) {
          xv[i1][k].resize(nDofs);
        }
        double s1 = static_cast<double>(i1);
        for (unsigned i = 0; i < nDofs; i++) {
          unsigned iDof  = sol->GetMesh()->GetSolutionDof(i, _elem, 2);    // global to global mapping between coordinates node and coordinate dof
          for (unsigned k = 0; k < _dim; k++) {
            xv[i1][k][i] = GetCoordinates(sol, k, iDof, s1);  // global extraction and local storage for the element coordinates
          }
        }
        //END extraction


        //BEGIN projection nodal to polynomial coefficients
        aX[i1].resize(solType + 1);
        for (unsigned j = 0; j < solType + 1; j++) {
          ProjectNodalToPolynomialCoefficients(aX[i1][j], xv[i1], elemType, j);
        }
        //END projection nodal to polynomial coefficients
      }


      //BEGIN find initial guess
      if (sol->GetIfFSI()) {
        for (unsigned i = 0; i < nDofs; i++) {
          for (unsigned k = 0; k < _dim; k++) {
            xv[0][k][i] = (1. - s) * xv[0][k][i] + s * xv[1][k][i];
          }
        }
      }
      GetClosestPointInReferenceElement (xv[0], _x, elemType, _xi);
    }


    std::vector < std::vector < std::vector < double > > > aXs;
    if (sol->GetIfFSI()) {
      InterpolatePolynomialCoefficients(aXs, aX[0], aX[1], s);
    }


    //BEGIN Inverse mapping loop
    for (unsigned j = 0; j < solType; j++) {

      std::vector < double > phi;
      std::vector < std::vector < double > > gradPhi;
      bool convergence = false;
      while (!convergence) {
        GetPolynomialShapeFunctionGradient(phi, gradPhi, _xi, elemType, solType);
        if (!sol->GetIfFSI()) {
          convergence = GetNewLocalCoordinates(_xi, _x, phi, gradPhi, aX[0][solType]);
        }
        else {
          convergence = GetNewLocalCoordinates(_xi, _x, phi, gradPhi, aXs[solType]);
        }
      }
    }
    //END Inverse mapping loop


//    std::cout << "ED USCIAMO" << std::endl;

  }



}











