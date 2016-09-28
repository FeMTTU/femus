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

const double initialGuess[6][3] = {
    {0., 0., 0.},
    {0.25, 0.25, 0.25},
    {1. / 3., 1. / 3., 0},
    {0., 0.},
    {1. / 3., 1. / 3.},
    {0.}
};

const double faceNormal[6][6][3] = {
    {{0., -1., 0.} , {1., 0., 0.}, {0., 1., 0.} , { -1., 0., 0.}, {0., 0., -1.} , {0., 0., 1.}},
    {{0., 0., -1}, {0., -1., 0.}, {1. / sqrt(3.), 1. / sqrt(3.), 1. / sqrt(3.)}, {-1., 0., 0.}},
    {{0., -1., 0.}, {1. / sqrt(2.), 1. / sqrt(2.), 0.}, {-1., 0., 0.}, {0., 0., -1}, {0., 0., 1}},
    {{0., -1.} , {1., 0.}, {0., 1.} , { -1., 0.}},
    {{0., -1.}, {1. / sqrt(2.), 1. / sqrt(2.)}, {-1., 0.}},
    {{ -1}, {1}}
};

const unsigned faceNumber[6] = {6,4,5,4,3,2};

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

    {   {   {{0, 1, 20, 0}, {1, 5, 20, 1}, {5, 4, 20, 5}, {4, 0, 20, 4}},  //hex
            {{1, 2, 21, 1}, {2, 6, 21, 2}, {6, 5, 21, 6}, {5, 1, 21, 5}},
            {{2, 3, 22, 2}, {3, 7, 22, 3}, {7, 6, 22, 7}, {6, 2, 22, 6}},
            {{3, 0, 23, 3}, {0, 4, 23, 0}, {4, 7, 23, 4}, {7, 3, 23, 7}},
            {{0, 3, 24, 0}, {3, 2, 24, 3}, {2, 1, 24, 2}, {1, 0, 24, 1}},
            {{4, 5, 25, 4}, {5, 6, 25, 5}, {6, 7, 25, 6}, {7, 4, 25, 7}}
        },
        {   {{0, 8, 20, 0}, {8, 1, 20, 8}, {1, 17, 20, 1}, {17, 5, 20, 17}, {5, 12, 20, 5}, {12, 4, 20, 12}, {4, 16, 20, 4 }, {16, 0, 20, 16}},
            {{1, 9, 21, 1}, {9, 2, 21, 9}, {2, 18, 21, 2}, {18, 6, 21, 18}, {6, 13, 21, 6}, {13, 5, 21, 13}, {5, 17, 21, 5}, {17, 1, 21, 17}},
            {{2, 10, 22, 2}, {10, 3, 22, 10}, {3, 19, 22, 3}, {19, 7, 22, 19}, {7, 14, 22, 7}, {14, 6, 22, 14}, {6, 18, 22, 6}, {18, 2, 22, 18}},
            {{3, 11, 23, 3}, {11, 0, 23, 11}, {0, 16, 23, 0}, {16, 4, 23, 16}, {4, 15, 23, 4}, {15, 7, 23, 15}, {7, 19, 23, 7}, {19, 3, 23, 19}},
            {{0, 11, 24, 0}, {11, 3, 24, 11}, {3, 10, 24, 3}, {10, 2, 24, 10}, {2, 9, 24, 2}, {9, 1, 24, 9}, {1, 8, 24, 1}, {8, 0, 24, 8}},
            {{4, 12, 25, 4}, {12, 5, 25, 12}, {5, 13, 25, 5}, {13, 6, 25, 13}, {6, 14, 25, 6}, {14, 7, 25, 14}, {7, 15, 25, 7}, {15, 4, 25, 15}}
        },
        {   {{0, 8, 20, 0}, {8, 1, 20, 8}, {1, 17, 20, 1}, {17, 5, 20, 17}, {5, 12, 20, 5}, {12, 4, 20, 12}, {4, 16, 20, 4 }, {16, 0, 20, 16}},
            {{1, 9, 21, 1}, {9, 2, 21, 9}, {2, 18, 21, 2}, {18, 6, 21, 18}, {6, 13, 21, 6}, {13, 5, 21, 13}, {5, 17, 21, 5}, {17, 1, 21, 17}},
            {{2, 10, 22, 2}, {10, 3, 22, 10}, {3, 19, 22, 3}, {19, 7, 22, 19}, {7, 14, 22, 7}, {14, 6, 22, 14}, {6, 18, 22, 6}, {18, 2, 22, 18}},
            {{3, 11, 23, 3}, {11, 0, 23, 11}, {0, 16, 23, 0}, {16, 4, 23, 16}, {4, 15, 23, 4}, {15, 7, 23, 15}, {7, 19, 23, 7}, {19, 3, 23, 19}},
            {{0, 11, 24, 0}, {11, 3, 24, 11}, {3, 10, 24, 3}, {10, 2, 24, 10}, {2, 9, 24, 2}, {9, 1, 24, 9}, {1, 8, 24, 1}, {8, 0, 24, 8}},
            {{4, 12, 25, 4}, {12, 5, 25, 12}, {5, 13, 25, 5}, {13, 6, 25, 13}, {6, 14, 25, 6}, {14, 7, 25, 14}, {7, 15, 25, 7}, {15, 4, 25, 15}}
        }
    },
    {   {   {{0, 2, 1, 0}},   //tet
            {{0, 1, 3, 0}},
            {{1, 2, 3, 1}},
            {{2, 0, 3, 2}}
        },
        {   {{0, 6, 10, 0}, {6, 2, 10, 6}, {2, 5, 10, 2}, {5, 1, 10, 5}, {1, 4, 10, 1}, {4, 0, 10, 4}},
            {{0, 4, 11, 0}, {4, 1, 11, 4}, {1, 8, 11, 1}, {8, 3, 11, 8}, {3, 7, 11, 3}, {7, 0, 11, 7}},
            {{1, 5, 12, 1}, {5, 2, 12, 5}, {2, 9, 12, 2}, {9, 3, 12, 9}, {3, 8, 12, 3}, {8, 1, 12, 8}},
            {{2, 6, 13, 2}, {6, 0, 13, 6}, {0, 7, 13, 0}, {7, 3, 13, 7}, {3, 9, 13, 3}, {9, 2, 13, 9}}
        },
        {   {{0, 6, 10, 0}, {6, 2, 10, 6}, {2, 5, 10, 2}, {5, 1, 10, 5}, {1, 4, 10, 1}, {4, 0, 10, 4}},
            {{0, 4, 11, 0}, {4, 1, 11, 4}, {1, 8, 11, 1}, {8, 3, 11, 8}, {3, 7, 11, 3}, {7, 0, 11, 7}},
            {{1, 5, 12, 1}, {5, 2, 12, 5}, {2, 9, 12, 2}, {9, 3, 12, 9}, {3, 8, 12, 3}, {8, 1, 12, 8}},
            {{2, 6, 13, 2}, {6, 0, 13, 6}, {0, 7, 13, 0}, {7, 3, 13, 7}, {3, 9, 13, 3}, {9, 2, 13, 9}}
        }
    },
    {   {   {{0, 1, 15, 0}, {1, 4, 15, 1}, {4, 3, 15, 4}, {3, 0, 15, 3}},  //wedge
            {{1, 2, 16, 1}, {2, 5, 16, 2}, {5, 4, 16, 5}, {4, 1, 16, 4}},
            {{2, 0, 17, 2}, {0, 3, 17, 0}, {3, 5, 17, 3}, {5, 2, 17, 5}},
            {{0, 2, 1, 0}},
            {{3, 4, 5, 3}}
        },
        {   {{0, 6, 15, 0 }, {6, 1, 15, 6}, {1, 13, 15, 1}, {13, 4, 15, 13}, {4, 9, 15, 4}, {9, 3, 15, 9}, {3, 12, 15, 3}, {12, 0, 15, 12}},
            {{1, 7, 16, 1}, {7, 2, 16, 7}, {2, 14, 16, 2}, {14, 5, 16, 14}, {5, 10, 16, 5}, {10, 4, 16, 10}, {4, 13, 16, 4}, {13, 1, 16, 13}},
            {{2, 8, 17, 2}, {8, 0, 17, 8}, {0, 12, 17, 0}, {12, 3, 17, 12}, {3, 11, 17, 3}, {11, 5, 17, 11}, {5, 14, 17, 5}, {14, 2, 17, 14}},
            {{0, 8, 18, 0}, {8, 2, 18, 8}, {2, 7 , 18, 2}, {7, 1, 18, 7} , {1, 6 , 18, 1} , {6, 0, 18, 6}},
            {{3, 9, 19, 3}, {9, 4, 19, 9}, {4, 10, 19, 4}, {10, 5, 19, 10}, {5, 11, 19, 5}, {11, 3, 19, 11}}
        },
        {   {{0, 6, 15, 0 }, {6, 1, 15, 6}, {1, 13, 15, 1}, {13, 4, 15, 13}, {4, 9, 15, 4}, {9, 3, 15, 9}, {3, 12, 15, 3}, {12, 0, 15, 12}},
            {{1, 7, 16, 1}, {7, 2, 16, 7}, {2, 14, 16, 2}, {14, 5, 16, 14}, {5, 10, 16, 5}, {10, 4, 16, 10}, {4, 13, 16, 4}, {13, 1, 16, 13}},
            {{2, 8, 17, 2}, {8, 0, 17, 8}, {0, 12, 17, 0}, {12, 3, 17, 12}, {3, 11, 17, 3}, {11, 5, 17, 11}, {5, 14, 17, 5}, {14, 2, 17, 14}},
            {{0, 8, 18, 0}, {8, 2, 18, 8}, {2, 7 , 18, 2}, {7, 1, 18, 7} , {1, 6 , 18, 1} , {6, 0, 18, 6}},
            {{3, 9, 19, 3}, {9, 4, 19, 9}, {4, 10, 19, 4}, {10, 5, 19, 10}, {5, 11, 19, 5}, {11, 3, 19, 11}}
        }
    },
    {   {   {{0, 1, 8, 0}, {1, 2, 8, 1}, {2, 3, 8, 1}, {3, 0, 8, 3}}, //quad
            {{0, 4, 8, 0}, {4, 1, 8, 4}, {1, 5, 8, 1}, {5, 2, 8, 5}, {2, 6, 8, 2}, {6, 3, 8, 6}, {3, 7, 8, 3}, {7, 0, 8, 7}},
            {{0, 4, 8, 0}, {4, 1, 8, 4}, {1, 5, 8, 1}, {5, 2, 8, 5}, {2, 6, 8, 2}, {6, 3, 8, 6}, {3, 7, 8, 3}, {7, 0, 8, 7}}
        }
    },
    {   {   {{0, 1, 2, 0}},
            {{0, 3, 6, 0}, {3, 1, 6, 3}, {1, 4, 6, 1}, {4, 2, 6, 4}, {2, 5, 6, 2}, {5, 0, 6, 5}},
            {{0, 3, 6, 0}, {3, 1, 6, 3}, {1, 4, 6, 1}, {4, 2, 6, 4}, {2, 5, 6, 2}, {5, 0, 6, 5}}
        }
    },
    {{{{}}}},
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
                nextElem[_iproc] = GetNextElement2D(iel);
            }
            else if(dim == 3) {
                nextElem[_iproc] = GetNextElement3D(iel);
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



const std::vector< std::vector <double> > inverseMatrix(const std::vector< std::vector <double> > &A, const unsigned dim) {

    std::vector < std::vector <double> > invA(dim);
    
    for(int i = 0; i < dim; i++) {
        invA[i].resize(dim);
    }
    
    double detA;

    if(dim == 2) {

        detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];

        if(detA == 0) {
            std::cout<< " ERROR: the matrix is singular " << std::endl;
            abort();
        }
        else {
            invA[0][0] = A[1][1] / detA;
            invA[0][1] = -A[0][1] / detA;
            invA[1][0] = -A[1][0] / detA;
            invA[1][1] = A[0][0] / detA;
        }
    }
    else if(dim == 3) {

        detA = (A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1])
               - (A[2][0] * A[1][1] * A[0][2] + A[2][1] * A[1][2] * A[0][0] + A[2][2] * A[1][0] * A[0][1]) ;

        if(detA == 0) {
            std::cout<< " ERROR: the matrix is singular " << std::endl;
            abort();
        }
        else {

            invA[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) / detA ;
            invA[0][1] = (A[0][2] * A[2][1] - A[2][2] * A[0][1]) / detA ;
            invA[0][2] = (A[0][1] * A[1][2] - A[1][1] * A[0][2]) / detA ;
            invA[1][0] = (A[1][2] * A[2][0] - A[2][2] * A[1][0]) / detA ;
            invA[1][1] = (A[0][0] * A[2][2] - A[2][0] * A[0][2]) / detA ;
            invA[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / detA ;
            invA[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) / detA ;
            invA[2][1] = (A[0][1] * A[2][0] - A[2][1] * A[0][0]) / detA ;
            invA[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) / detA ;

        }
    }
    else {
        std::cout << " ERROR: the matrix is neither 2x2 nor 3x3 so we cannot use this function " <<std::endl;
        abort();
    }

    return invA;
}


unsigned Marker::FastForward(const unsigned &iel){
              //std::cout<<"aaaaaaaaaaaaaaaaaaa\n";

            unsigned dim =  _mesh->GetDimension();
            unsigned nDofs = _mesh->GetElementDofNumber(iel, _solType);
            short unsigned ielType = _mesh->GetElementType(iel);

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
            ProjectNodalToPolynomialCoefficients(a, xv, ielType, _solType);

            //END projection


            //BEGIN inverse mapping search
            //std::vector <double> xi(dim, 0.);

            std::vector < double > phi;
            std::vector < std::vector < double > > gradPhi;

            std::vector < double > xi(dim);
            for(int k = 0; k < dim; k++) {
                xi[k] = initialGuess[ielType][k];
            }

            GetPolynomialShapeFunctionGradient(phi, gradPhi, xi, ielType, _solType);


            std::vector < double > v(dim, 0.);
            std::vector < std::vector < double > > J(dim);

            for(int k = 0; k < dim; k++) {
                J[k].assign(dim, 0.);
            }

            for(int k = 0; k < dim; k++) {
                for(int i = 0; i < nDofs; i++) {
                    v[k] -= a[k][i] * phi[i];   // why a - sign ?????

                    for(int i1 = 0; i1 < dim; i1++) {
                        J[k][i1] += a[k][i] * gradPhi[i][i1];
                    }
                }
                v[k] += _x[k];  // why - _x ??????
            }

           // std::cout << " v[0] = " << v[0] << " "<< " v[1] = " << v[1] <<std::endl;

            std::vector < std::vector < double > >  Jm1(dim);

            for(int i1 = 0; i1 < dim; i1++) {
                Jm1[i1].resize(dim);
            }

            Jm1 = inverseMatrix(J, dim);

            std::vector < double > vt(dim, 0.);
            for(unsigned i = 0; i < dim; i++) {
                for(unsigned j = 0; j < dim; j++) {
                    vt[i] += Jm1[i][j] * v[j];
                }
            }


           // std::cout << " vt[0] " <<  vt[0] << " " << " vt[1] " << vt[1] <<std::endl;

            double maxProjection = 0.;
            unsigned faceIndex = 0;

            for(unsigned jface=0; jface < faceNumber[ielType]; jface++) {
                double projection = 0.;
                for( unsigned k = 0; k < dim; k++) {
                    projection += vt[k] * faceNormal[ielType][jface][k];
                }
                if(projection > maxProjection) {
		  maxProjection = projection;
                  //std::cout<< " jface = " << jface << " projection= " << projection <<std::endl;  
		  faceIndex = jface;
                }
            }

            
            unsigned nextElem = (_mesh->el->GetFaceElementIndex(iel, faceIndex) -1);
            //nextElementFound = true;

           // std::cout<<"I want to go to "<< nextElem <<std::endl;
	    
	    return nextElem;
  
}



unsigned Marker::GetNextElement2D(const unsigned &currentElem) {
  


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
            xv[k].resize(facePointNumber[currentElementType][_solType]);
        }

        for(unsigned i = 0; i < facePointNumber[currentElementType][_solType]; i++) {
            unsigned inodeDof  = _mesh->GetSolutionDof(facePoints[currentElementType][_solType][i], currentElem, 2);
           // std::cout << "inodeDof = " << inodeDof << std::endl;

            for(unsigned k = 0; k < dim; k++) {
                xv[k][i] = (*_mesh->_topology->_Sol[k])(inodeDof) - _x[k];
            }
        }

        // rescaling coordinates to properly handle different scales of meshes
        double length = 0.;
        double sum = 0.;

        for(unsigned i = 0; i < facePointNumber[currentElementType][_solType] - 1; i++) {
            for(unsigned k = 0; k < dim; k++) {
                sum += (xv[k][i + 1] - xv[k][i]) * (xv[k][i + 1] - xv[k][i]);
            }

            length += sqrt(sum);
        }

        length /= facePointNumber[currentElementType][_solType];
       // std::cout << "length= " << length << std::endl;

        for(unsigned k = 0; k < dim; k++) {
            xcc[k] = xc[k] / length;

            for(unsigned i = 0; i < facePointNumber[currentElementType][_solType]; i++) {
                xv[k][i] /= length;
            }
        }


        //Find a ball centered in (xcc[0] , xcc[1]) that inscribes the element
        double radius2 = 0.;
        for(unsigned i = 0; i < facePointNumber[currentElementType][_solType] - 1; i++) {
            double iradius2 = 0.;
            for(unsigned k = 0; k < dim; k++) {
                iradius2 += (xcc[k] - xv[k][i]) * (xcc[k] - xv[k][i]);
            }
            if(radius2 < iradius2) radius2 = iradius2;
        }
        radius2 *= 1.2;

      //  std::cout << "radius2 = " <<radius2 << " " << " xcc[0]*xcc[0] + xcc[1]*xcc[1] = " << xcc[0]*xcc[0] + xcc[1]*xcc[1] << " " ;

//       //TEST if the marker is inside a ball centered in (xcc[0] , xcc[1]) and given radius
        //std::cout<<"aaaaaaaaaaaaaaaaaaa\n";
        if(radius2 < xcc[0]*xcc[0] + xcc[1]*xcc[1]) {  // project direction

//             //std::cout<<"aaaaaaaaaaaaaaaaaaa\n";
// 
//             unsigned iel = currentElem;
//             unsigned solType = _solType;
// 
//             unsigned dim =  _mesh->GetDimension();
//             unsigned nDofs = _mesh->GetElementDofNumber(iel, solType);
//             short unsigned ielType = _mesh->GetElementType(iel);
// 
//             //BEGIN extraction nodal coordinate values
//             std::vector< std::vector < double > > xv(dim);
// 
//             for(unsigned k = 0; k < dim; k++) {
//                 xv[k].resize(nDofs);
//             }
// 
//             for(unsigned i = 0; i < nDofs; i++) {
//                 unsigned iDof  = _mesh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
// 
//                 for(unsigned k = 0; k < dim; k++) {
//                     xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
//                 }
//             }
//             //END extraction
// 
// 
//             //BEGIN projection nodal to polynomial coefficients
//             std::vector < std::vector < double > > a;
//             ProjectNodalToPolynomialCoefficients(a, xv, ielType, solType);
// 
//             //END projection
// 
// 
//             //BEGIN inverse mapping search
//             //std::vector <double> xi(dim, 0.);
// 
//             std::vector < double > phi;
//             std::vector < std::vector < double > > gradPhi;
// 
//             std::vector < double > xi(dim);
//             for(int k = 0; k < dim; k++) {
//                 xi[k] = initialGuess[ielType][k];
//             }
// 
//             GetPolynomialShapeFunctionGradient(phi, gradPhi, xi, ielType, solType);
// 
// 
//             std::vector < double > v(dim, 0.);
//             std::vector < std::vector < double > > J(dim);
// 
//             for(int k = 0; k < dim; k++) {
//                 J[k].assign(dim, 0.);
//             }
// 
//             for(int k = 0; k < dim; k++) {
//                 for(int i = 0; i < nDofs; i++) {
//                     v[k] -= a[k][i] * phi[i];   // why a - sign ?????
// 
//                     for(int i1 = 0; i1 < dim; i1++) {
//                         J[k][i1] += a[k][i] * gradPhi[i][i1];
//                     }
//                 }
//                 v[k] += _x[k];  // why - _x ??????
//             }
// 
//            // std::cout << " v[0] = " << v[0] << " "<< " v[1] = " << v[1] <<std::endl;
// 
//             std::vector < std::vector < double > >  Jm1(dim);
// 
//             for(int i1 = 0; i1 < dim; i1++) {
//                 Jm1[i1].resize(dim);
//             }
// 
//             Jm1 = inverseMatrix(J, dim);
// 
//             std::vector < double > vt(dim, 0.);
//             for(unsigned i = 0; i < dim; i++) {
//                 for(unsigned j = 0; j < dim; j++) {
//                     vt[i] += Jm1[i][j] * v[j];
//                 }
//             }
// 
// 
//            // std::cout << " vt[0] " <<  vt[0] << " " << " vt[1] " << vt[1] <<std::endl;
// 
//             double maxProjection = 0.;
//             unsigned faceIndex = 0;
// 
//             for(unsigned jface=0; jface < faceNumber[ielType]; jface++) {
//                 double projection = 0.;
//                 for( unsigned k = 0; k < dim; k++) {
//                     projection += vt[k] * faceNormal[ielType][jface][k];
//                 }
//                 if(projection > maxProjection) {
// 		  maxProjection = projection;
//                   //std::cout<< " jface = " << jface << " projection= " << projection <<std::endl;  
// 		  faceIndex = jface;
//                 }
//             }
// 
//             
//             nextElem = (_mesh->el->GetFaceElementIndex(currentElem, faceIndex) -1);
// 
//            // std::cout<<"I want to go to "<< nextElem <<std::endl;
	     
	
	  nextElem = FastForward(currentElem);
	  nextElementFound = true;
	  
        }

        else {
        //if(true) {
            //BEGIN look for face intersection

            for(unsigned i = 0 ; i < facePointNumber[currentElementType][_solType] - 1; i++) {

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
                            xvr[k].resize(facePointNumber[currentElementType][_solType]);
                        }

                        //now we have to determine if r is inside edge i
                        for(unsigned j = 0; j < facePointNumber[currentElementType][_solType]; j++) {
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
                                unsigned nodeIndex = (_solType == 0) ? i : i / 2;


                                /*
                                             if(i % 2 == 0 && i != facePointNumber[currentElementType][_solType]) nodeIndex = i / 2 ;
                                             else if(i == facePointNumber[currentElementType][_solType]) nodeIndex = (i - 2) / 2 ;
                                             else if(i % 2 != 0) nodeIndex = (i - 1) / 2 ;*/

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

                                unsigned nodeIndex = (_solType == 0) ? i : i / 2;
//                   unsigned nodeIndex;
//
// //                   if(i % 2 == 0 && i != facePointNumber[currentElementType][solType]) nodeIndex = i / 2 ;
// //                   else if(i == facePointNumber[currentElementType][solType]) nodeIndex = (i - 2) / 2 ;
// //                   else if(i % 2 != 0) nodeIndex = (i - 1) / 2 ;

                                nextElem = (_mesh->el->GetFaceElementIndex(currentElem, nodeIndex) - 1);
                                nextElementFound = true;
                                break;
                            }
                        }
                    } // closes the " if t < 1 "
                } // closes the else
            } //closes the for on the nodes
            //END look for face intersection
        }
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





unsigned Marker::GetNextElement3D(const unsigned & currentElem) {


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

            for(unsigned itri = 0; itri < trianglesPerFace[currentElementType][_solType][iface]; itri ++) {

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
                    unsigned itriDof  = _mesh->GetSolutionDof(faceTriangleNodes[currentElementType][_solType][iface][itri][i], currentElem, 2);
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


bool SPDCheck2D(const std::vector< std::vector <double> > &A) {
    bool SPD = true;

    if(A[0][1] != A[1][0]) {
        SPD = false;
        std::cout << "The 2D matrix is not symmetric" << std::endl;
    }

    else if(A[0][0] < 1.0e-8) {
        SPD = false;
    }
    else {
        double lambda = A[1][1] - ((A[0][1] * A[0][1]) / A[0][0]);
        if(lambda < 0 || fabs(lambda) < 1.0e-8) {
            SPD = false;
        }
    }

    if(SPD == true) {
        std::cout << "The 2D matrix is SPD" << std::endl;
    }

    else {
        std::cout << "The 2D matrix is not SPD" << std::endl;
    }

    return SPD;

}




bool SPDCheck3D(const std::vector< std::vector <double> > &A) {

    bool SPD = true;
    bool notSymm = false;

    if(A[0][2] != A[2][0] || A[1][2] != A[2][1] || A[0][1] != A[1][0]) {
        notSymm = true;
        std::cout << "The 3D matrix is not symmetric" << std::endl;
    }

    std::vector < std::vector <double> > B(2);
    for(int i = 0; i < 2; i++) {
        B[i].resize(2);
    }

    B[0][0] = A[0][0];
    B[0][1] = A[0][1];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1];



    if(SPDCheck2D(B) == false || notSymm == true) {
        SPD = false;
    }

    else {
        double l00 = sqrt(A[0][0]);
        double l11 = sqrt(A[1][1] - ((A[0][1] * A[0][1]) / A[0][0]));
        double l01 = A[0][1] / l00;

        double detL =  l00 * l11;
        std::vector < std::vector < double > > Lm1(2);
        for(int i = 0; i < 2; i++) {
            Lm1[i].resize(2);
        }

        Lm1[0][0] = (1 / detL) * l11;
        Lm1[1][0] = -(1 / detL) * l01;
        Lm1[0][1] = 0. ;
        Lm1[1][1] = (1 / detL) * l00 ;

        std::vector < double > K(2);

        K[0] = Lm1[0][0] * A[0][2] + Lm1[0][1] * A[1][2];
        K[1] = Lm1[1][0] * A[0][2] + Lm1[1][1] * A[1][2];

        double KK = sqrt(K[0] * K[0] + K[1] * K[1]) ;

        if(A[2][2] - KK < 0 || fabs(A[2][2] - KK) < 1.0e-8) {
            SPD = false;
        }
    }

    if(SPD == true) {
        std::cout << "The 3D matrix is SPD" << std::endl;
    }

    else {
        std::cout << "The 3D matrix is not SPD" << std::endl;
    }

    return SPD;

}



bool GetNewLocalCoordinatesHess(std::vector <double> &xi, const std::vector< double > &x, const std::vector <double> &phi,
                                const std::vector < std::vector <double > > &gradPhi, const std::vector < std::vector < std::vector <double> > > hessPhi,
                                const std::vector < std::vector <double > > &a, const unsigned & dim, const unsigned & nDofs) {

    bool convergence = false;
    std::vector < double > xp(dim, 0.);
    std::vector < std::vector < double > > gradXp(dim);
    std::vector < std::vector < std::vector < double > > > hessXp(dim);

//     for(int k = 0; k < nDofs; k++) {
//       bool isSPD;
//       if(dim == 2) {
// 	isSPD = SPDCheck2D(hessPhi[k]);
//       }
//       else if(dim == 3) {
// 	isSPD = SPDCheck3D(hessPhi[k]);
//       }
//     }



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

//     for(int k = 0; k < dim; k++) {
//       bool isSPD;
//       if(dim == 2) {
// 	isSPD = SPDCheck2D(hessXp[k]);
//       }
//       else if(dim == 3) {
// 	isSPD = SPDCheck3D(hessXp[k]);
//       }
//     }



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



//     bool isSPD;
//     if(dim == 2) {
//       isSPD = SPDCheck2D(hessF);
//     }
//     else if(dim == 3) {
//       isSPD = SPDCheck3D(hessF);
//     }


//     if(isSPD == false) {
//       double norm = 0.;
//       for(int i1 = 0; i1 < dim; i1++) {
//         double rowSum = 0.;
//         for(int i2 = 0; i2 < dim; i2++) {
//           rowSum += fabs(hessF[i1][i2]);
//         }
//         if(rowSum > norm) norm  = rowSum;
//       }
//
//       for(int i1 = 0; i1 < dim; i1++) {
//         hessF[i1][i1] += norm;
//       }
//     }
//
//     if(dim == 2) {
//       isSPD = SPDCheck2D(hessF);
//     }
//     else if(dim == 3) {
//       isSPD = SPDCheck3D(hessF);
//     }


    std::vector < std::vector < double > >  hessFm1(dim);

    for(int i1 = 0; i1 < dim; i1++) {
        hessFm1[i1].resize(dim);
    }

    
    hessFm1 = inverseMatrix(hessF, dim);
    
    
//     if(dim == 2) {
//         double det = hessF[0][0] * hessF[1][1] - hessF[0][1] * hessF[1][0];
//         hessFm1[0][0] = hessF[1][1] / det;
//         hessFm1[0][1] = -hessF[0][1] / det;
//         hessFm1[1][0] = -hessF[1][0] / det;
//         hessFm1[1][1] = hessF[0][0] / det;
//     }
//     else if(dim == 3) {
//         double det = (hessF[0][0] * hessF[1][1] * hessF[2][2] + hessF[0][1] * hessF[1][2] * hessF[2][0] + hessF[0][2] * hessF[1][0] * hessF[2][1])
//                      - (hessF[2][0] * hessF[1][1] * hessF[0][2] + hessF[2][1] * hessF[1][2] * hessF[0][0] + hessF[2][2] * hessF[1][0] * hessF[0][1]) ;
// 
//         hessFm1[0][0] = (hessF[1][1] * hessF[2][2] - hessF[2][1] * hessF[1][2]) / det ;
//         hessFm1[0][1] = (hessF[0][2] * hessF[2][1] - hessF[2][2] * hessF[0][1]) / det ;
//         hessFm1[0][2] = (hessF[0][1] * hessF[1][2] - hessF[1][1] * hessF[0][2]) / det ;
//         hessFm1[1][0] = (hessF[1][2] * hessF[2][0] - hessF[2][2] * hessF[1][0]) / det ;
//         hessFm1[1][1] = (hessF[0][0] * hessF[2][2] - hessF[2][0] * hessF[0][2]) / det ;
//         hessFm1[1][2] = (hessF[0][2] * hessF[1][0] - hessF[0][0] * hessF[1][2]) / det ;
//         hessFm1[2][0] = (hessF[1][0] * hessF[2][1] - hessF[2][0] * hessF[1][1]) / det ;
//         hessFm1[2][1] = (hessF[0][1] * hessF[2][0] - hessF[2][1] * hessF[0][0]) / det ;
//         hessFm1[2][2] = (hessF[0][0] * hessF[1][1] - hessF[1][0] * hessF[0][1]) / det ;
//     }

    double delta2 = 0.;

    for(int i1 = 0; i1 < dim; i1++) {
        double deltak = 0.;

        for(int i2 = 0; i2 < dim; i2++) {
            deltak -= hessFm1[i1][i2] * gradF[i2];
        }

        xi[i1] += deltak;
        delta2 += deltak * deltak;
    }

//     for(int k = 0; k < dim; k++) {
//       std::cout << "xT[" << k << "]= " << xp[k] <<  " ";
//     }
//     std::cout << std::endl;

    if(delta2 < 1.0e-9) {
        convergence = true;
    }

    return convergence;
}



bool GetNewLocalCoordinates(std::vector <double> &xi, const std::vector< double > &x, const std::vector <double> &phi,
                            const std::vector < std::vector <double > > &gradPhi, const std::vector < std::vector < std::vector <double> > > hessPhi,
                            const std::vector < std::vector <double > > &a, const unsigned & dim, const unsigned & nDofs) {

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
            }
        }
        F[k] -= x[k];
    }


    std::vector < std::vector < double > >  Jm1(dim);

    for(int i1 = 0; i1 < dim; i1++) {
        Jm1[i1].resize(dim);
    }

    Jm1 = inverseMatrix(J,dim);

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

    if(delta2 < 1.0e-9) {
        convergence = true;
    }

    return convergence;
}



void Marker::InverseMapping(const unsigned & iel, const unsigned & solType,
                            const std::vector< double > &x, std::vector< double > &xi) {



    unsigned dim =  _mesh->GetDimension();
    unsigned nDofs = _mesh->GetElementDofNumber(iel, solType);
    short unsigned ielType = _mesh->GetElementType(iel);

    std::cout << solType << " " << nDofs << std::endl;

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

//     for(int k=0; k<dim; k++) {
//         for(int i=0; i<nDofs; i++) {
//             std::cout << a[k][i] <<std::endl;
//         }
//         std::cout << std::endl;
//     }

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

        if(solType == 0) {
            convergence = GetNewLocalCoordinates(xi, x, phi, gradPhi, hessPhi, a, dim, nDofs);
        }
        else {
            counter++;
            //convergence = GetNewLocalCoordinates(xi, x, phi, gradPhi, hessPhi, a, dim, nDofs);
            convergence = GetNewLocalCoordinatesHess(xi, x, phi, gradPhi, hessPhi, a, dim, nDofs);
        }
        for(unsigned k = 0; k < dim; k++) {
            std::cout << "xi[" << k << "] = " << xi[k] << " ";
        }
        std::cout << std::endl;
    }
    //END inverse mapping search

    return;


}











void Marker::InverseMappingTEST(std::vector < double > &x) {
    unsigned dim = _mesh->GetDimension();

    for(int solType = 0; solType < 3; solType++) {
        std::cout << "\n\n--------------------------------------------------" << std::endl;
        std::cout << "solType = " << solType << std::endl;

        for(int iel = _mesh->_elementOffset[_iproc]; iel < _mesh->_elementOffset[_iproc + 1]; iel++) {
            //for(int iel = 494; iel < 495; iel++) {
            std::cout << "iel = " << iel << std::endl;
            std::cout << "--------------------------------------------------\n" << std::endl;

            unsigned nDofs = _mesh->GetElementDofNumber(iel, solType);
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
                std::cout << j << std::endl;
                std::vector < double > xiT(dim);

                for(unsigned k = 0; k < dim; k++) {
                    xiT[k] = *(_mesh->_finiteElement[ielType][0]->GetBasis()->GetXcoarse(j) + k);
                }

                // This is the test
                std::vector < double > xi(dim);
                for(int k = 0; k < dim; k++) {
                    xi[k] = initialGuess[ielType][k];
                }

                for(int k = 0; k < dim; k++) {
                    std::cout << "xv[" << k << "]= " << xv[j][k] <<  " ";
                }

                std::cout << std::endl;


                for(int itype = 0; itype <= solType; itype++) {
                    InverseMapping(iel, itype, xv[j], xi);
                    std::cout << std::endl;
                }

// 	  InverseMapping(iel, 0, xv[j], xi);
// 	  if(solType > 0){
// 	    InverseMapping(iel, solType, xv[j], xi);
// 	    std::cout << std::endl;
// 	  }


                bool test = true;
                for(int k = 0; k < dim; k++) {
                    std::cout << "xiT[" << k << "]= " << xiT[k] <<  " xi[" << k << "]= " << xi[k];
                    std::cout << " error: " << xiT[k] - xi[k] << std::endl;
                    if(fabs(xiT[k] - xi[k]) > 1.0e-3) {
                        test = false;
                    }
                }
                if(test == false) {
                    std::cout << "Inverse map test failed " << std::endl;
                    abort();
                }
                std::cout << "--------------------------------------------------\n" << std::endl;
            }
        }
    }

    std::cout << "total number of calls= " << counter << std::endl;
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

}







