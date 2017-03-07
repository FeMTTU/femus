/*=========================================================================

 Program: FEMUS
 Module: Line
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
#include <cmath>
#include "PolynomialBases.hpp"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using boost::math::ellint_1;
using boost::math::ellint_2;

namespace femus
{

  const double Line::_a[4][4][4] = {
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

  const double Line::_b[4][4] = {
    {1.}, // first order
    {0.5, 0.5}, // second order (Heun's)
    {1. / 6., 2. / 3., 1. / 6.}, // third-order method
    {1. / 6., 1. / 3., 1. / 3., 1. / 6.} // fourth-order method
  };

  const double Line::_c[4][4] = {
    {0.}, // first order
    {0., 1.}, // second order (Heun's)
    {0., 0.5, 1.}, // third-order method
    {0., 0.5, 0.5, 1.} // fourth-order method
  };


  Line::Line(const std::vector < std::vector < double > > x,
             const std::vector <MarkerType> &markerType,
             Solution *sol, const unsigned & solType)
  {

    _sol = sol;
    _mesh = _sol->GetMesh();

    _time.assign(10, 0);

    _size = x.size();

    std::vector < Marker*> particles(_size);

    _dim = _mesh->GetDimension();

    _markerOffset.resize(_nprocs + 1);
    _markerOffset[_nprocs] = _size;

    _particles.resize(_size);
    _printList.resize(_size);

    for (unsigned j = 0; j < _size; j++) {
      particles[j] = new Marker(x[j], markerType[j], _sol, solType, true);
    }



    //BEGIN TEST ASSIGNATION

//     std::vector <double> xtry(_dim);
//     std::vector <double> xitry(_dim);
//
//
//     for(unsigned j = 0; j < _size; j++) {
//
//       unsigned proc;
//       unsigned elem;
//       proc = particles[j]->GetMarkerProc();
//       elem = particles[j]->GetMarkerElement();
//
//       std::cout << "proc = " << proc << " , " << "elem = " << elem <<
//                 " , " << std::endl;
//       particles[j]->GetMarkerCoordinates(xtry);
//       for(unsigned i = 0 ; i < _dim; i++) {
//         std::cout << "x[" << i << "]=" << xtry[i] << std::endl;
//       }
//       particles[j]->GetMarkerLocalCoordinatesLine(xitry);
//       for(unsigned i = 0 ; i < _dim; i++) {
//         std::cout << "xi[" << i << "]=" << xitry[i] << std::endl;
//       }
//     }
//
//
//     for(unsigned j = 0; j < _size; j++) {
//       _particles[j] = particles[j];
//     }
//
//     std::cout << " ------------------------------------------------------------------------------------------------ " << std::endl;
//
//     for(unsigned j = 0; j < _size; j++) {
//
//       unsigned proc;
//       unsigned elem;
//       proc = _particles[j]->GetMarkerProc();
//       elem = _particles[j]->GetMarkerElement();
//       std::cout << "proc = " << proc << " , " << "elem = " << elem <<
//                 " , " << std::endl;
//       _particles[j]->GetMarkerCoordinates(xtry);
//       for(unsigned i = 0 ; i < _dim; i++) {
//         std::cout << "x[" << i << "]=" << xtry[i] << std::endl;
//       }
//       _particles[j]->GetMarkerLocalCoordinatesLine(xitry);
//       for(unsigned i = 0 ; i < _dim; i++) {
//         std::cout << "xi[" << i << "]=" << xitry[i] << std::endl;
//       }
//     }


    //END TEST ASSIGNATION


//     std::cout << "FIRST OF ALL" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc;
//       particles[i]->GetMarkerProcLine(proc);
//       unsigned elem;
//       particles[i]->GetMarkerElementLine(elem);
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " << std::endl;
//     }

    //BEGIN reorder the markers by proc
    unsigned counter = 0;
    for (unsigned iproc = 0; iproc < _nprocs; iproc++) {
      _markerOffset[iproc] = counter;
      //unsigned  offsetCounter = counter;
      for (unsigned j = 0; j < _size; j++) {
        unsigned markerProc = particles[j]->GetMarkerProc(_sol);
        if (markerProc == iproc) {
          _particles[counter] = particles[j];
          _printList[j] = counter;
          counter++;
        }
      }
    }
    //END reorder the markers by proc




//     std::cout << "AFTER THE REORDERING BY PROC" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc;
//       _particles[i]->GetMarkerProcLine(proc);
//       unsigned elem;
//       _particles[i]->GetMarkerElementLine(elem);
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " << std::endl;
//     }
//
//     for(unsigned iproc = 0; iproc <= _nprocs; iproc++) {
//       std::cout << "_markerOffset[" << iproc << "]= " << _markerOffset[iproc] << std::endl;
//     }

//         for(unsigned i = 0; i < _size; i++) {
//           std::cout << "_printList[ " << i << "] = " << _printList[i] << std::endl;
//         }
//
//
//         std::cout << " ------------------------------------------------------------------------------------------------ " << std::endl;

//     //BEGIN reorder markers also by element
//
//     for(unsigned j = 0; j < _size; j++) {
//       particles[j] = _particles[j];
//     }
//
//     unsigned* elementList = new unsigned [_sol->GetMesh()->GetNumberOfElements()];
//
//     std::cout << "number of elements = " <<  _sol->GetMesh()->GetNumberOfElements() << std::endl;
//
//     std::vector < unsigned> printList(_size);
//
//     //flags to see if we already ordered by that element, 0 = not considered yet
//     for(unsigned iFlag = 0; iFlag < _mesh->GetNumberOfElements(); iFlag++) {
//       elementList[iFlag] = 0;
//     }
//
//     printList = _printList;
//
//     for(unsigned iproc = 0; iproc < _nprocs; iproc++) {
//
//       counter = 0;
//
//       for(unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {
//
//         unsigned jel;
//         jel = particles[jp]->GetMarkerElement();
//
// 	std::cout<< " Marker in element jel = " << jel << std::endl;
//
//         if(elementList[jel] == 0) {
//
//           elementList[jel] = 1;
//
//           _particles[_markerOffset[iproc] + counter] = particles[jp];
//           for(unsigned iList = 0; iList < _size; iList++) {
//             if(printList[iList] == jp) {
//               _printList[iList] = _markerOffset[iproc] + counter;
//               break;
//             }
//           }
//           counter++;
//
//
//           for(unsigned ip = _markerOffset[iproc]; ip < _markerOffset[iproc + 1]; ip++) {
//             unsigned iel;
//             iel = particles[ip]->GetMarkerElement();
//             if(ip != jp && iel == jel) {
//               //std::cout << " jel =" << jel << " , " << "ip = " << ip << " , " << "jp = " << jp <<  std::endl;
//               elementList[iel] = 1;
//               _particles[_markerOffset[iproc] + counter] = particles[ip];
//               for(unsigned iList = 0; iList < _size; iList++) {
//                 if(printList[iList] == ip) {
//                   _printList[iList] = _markerOffset[iproc] + counter;
//                   break;
//                 }
//               }
//               counter++;
//             }
//           }
//         }
//       }
//     }
//
//     //END reorder markers also by element

//         for(unsigned i = 0; i < _size; i++) {
//           std::cout << "_printList[ " << i << "] = " << _printList[i] << std::endl;
//         }


//     std::cout << "AFTER THE REORDERING BY ELEM" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc;
//       _particles[i]->GetMarkerProcLine(proc);
//       unsigned elem;
//       _particles[i]->GetMarkerElementLine(elem);
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " << "Pointer: " << _particles[i] << std::endl;
//     }


    _line.resize(_size + 1);

    for (unsigned j = 0; j < _size; j++) {
      _particles[_printList[j]]->GetMarkerCoordinates(_line[j]);
    }
    _particles[_printList[0]]->GetMarkerCoordinates(_line[_size]);

    // delete[] elementList;
    // std::vector < unsigned > ().swap(printList);


//     for(unsigned i=0;i<_size;i++){
//       if(_printList[i] == 3395){
// 	std::cout << i<< std::endl;
//       }
//     }
//
//     exit(0);
  };

  Line::~Line()
  {
    for (unsigned j = 0; j < _size; j++) {
      //std::cout << j << " " << _particles[j] << std::endl<< std::flush;
      delete _particles[j];
    }
  }

  void Line::UpdateLine()
  {

    std::vector < Marker*> particles(_size);
    std::vector < unsigned> printList(_size);

    for (unsigned j = 0; j < _size; j++) {
      particles[j] = _particles[j];
//       printList[j] = _printList[j];
    }

    printList = _printList;


    //BEGIN TEST ASSIGNATION

//     std::vector <double> xtry(_dim);
//     std::vector <double> xitry(_dim);
//
//
//     for(unsigned j = 0; j < _size; j++) {
//
//       unsigned proc;
//       unsigned elem;
//       proc = particles[j]->GetMarkerProc();
//       elem = particles[j]->GetMarkerElement();
//
//       std::cout << "proc = " << proc << " , " << "elem = " << elem <<
//                 " , " << std::endl;
//       particles[j]->GetMarkerCoordinates(xtry);
//       for(unsigned i = 0 ; i < _dim; i++) {
//         std::cout << "x[" << i << "]=" << xtry[i] << std::endl;
//       }
//       particles[j]->GetMarkerLocalCoordinatesLine(xitry);
//       for(unsigned i = 0 ; i < _dim; i++) {
//         std::cout << "xi[" << i << "]=" << xitry[i] << std::endl;
//       }
//     }
//
//
//     for(unsigned j = 0; j < _size; j++) {
//       _particles[j] = particles[j];
//     }
//
//     std::cout << " ------------------------------------------------------------------------------------------------ " << std::endl;
//
//     for(unsigned j = 0; j < _size; j++) {
//
//       unsigned proc;
//       unsigned elem;
//       proc = _particles[j]->GetMarkerProc();
//       elem = _particles[j]->GetMarkerElement();
//       std::cout << "proc = " << proc << " , " << "elem = " << elem <<
//                 " , " << std::endl;
//       _particles[j]->GetMarkerCoordinates(xtry);
//       for(unsigned i = 0 ; i < _dim; i++) {
//         std::cout << "x[" << i << "]=" << xtry[i] << std::endl;
//       }
//       _particles[j]->GetMarkerLocalCoordinatesLine(xitry);
//       for(unsigned i = 0 ; i < _dim; i++) {
//         std::cout << "xi[" << i << "]=" << xitry[i] << std::endl;
//       }
//     }


    //END TEST ASSIGNATION







//     std::cout << "FIRST OF ALL ----------------------- UPDATE LINE" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc;
//       unsigned elem;
//       std::vector <double> chiappe;
//       chiappe = _particles[_printList[i]]->GetIprocMarkerCoordinates();
//       proc = _particles[_printList[i]]->GetMarkerProc();
//       elem = _particles[_printList[i]]->GetMarkerElement();
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " ";
//       for(unsigned i = 0; i < _dim; i++) {
//         std::cout << "x[" << i << "]=" << chiappe[i] << " ";
//       }
//       std::cout << std::endl;
//     }

    //BEGIN reorder the markers by proc
    unsigned counter = 0;
    for (unsigned iproc = 0; iproc < _nprocs; iproc++) {
      _markerOffset[iproc] = counter;
      for (unsigned j = 0; j < _size; j++) {
        unsigned markerProc = particles[j]->GetMarkerProc(_sol);
        if (markerProc == iproc) {
          _particles[counter] = particles[j];
          for (unsigned iList = 0; iList < _size; iList++) {
            if (printList[iList] == j) {
              _printList[iList] = counter;
              break;
            }
          }
          counter++;
        }
      }
    }
    //END reorder the markers by proc

//     for(unsigned iproc = 0; iproc <= _nprocs; iproc++) {
//       std::cout << "_markerOffset[" << iproc << "]= " << _markerOffset[iproc] << std::endl;
//     }

//     std::cout << "-----------------------------------------  UPDATE LINE  BY PROCS  -----------------------------------------" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc;
//       unsigned elem;
//       proc = _particles[i]->GetMarkerProc();
//       elem = _particles[i]->GetMarkerElement();
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " << std::endl;
//     }



//     //BEGIN reorder markers also by element
//
//     for(unsigned j = 0; j < _size; j++) {
//       particles[j] = _particles[j];
//       //printList[j] = _printList[j];
//     }
//
//     printList = _printList;
//
//     unsigned meshElements;
//     particles[0]->GetNumberOfMeshElements(meshElements);
//     unsigned* elementList = new unsigned [meshElements];
//     bool someMarkersOutsideDomain = false;
//
//     //flags to see if we already ordered by that element, 0 = not considered yet
//     for(unsigned iFlag = 0; iFlag < meshElements; iFlag++) {
//       elementList[iFlag] = 0;
//     }
//
//     for(unsigned iproc = 0; iproc < _nprocs; iproc++) {
//
//       counter = 0;
//
//       for(unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {
//
//         unsigned jel;
//         jel = particles[jp]->GetMarkerElement();
//
//         if(jel != UINT_MAX) {
//
//           if(elementList[jel] == 0) {
//
//             elementList[jel] = 1;
//
//             _particles[_markerOffset[iproc] + counter] = particles[jp];
//             for(unsigned iList = 0; iList < _size; iList++) {
//               if(printList[iList] == jp) {
//                 _printList[iList] = _markerOffset[iproc] + counter;
//                 break;
//               }
//             }
//             counter++;
//
//
//             for(unsigned ip = _markerOffset[iproc]; ip < _markerOffset[iproc + 1]; ip++) {
//               unsigned iel;
//               iel = particles[ip]->GetMarkerElement();
//               if(ip != jp && iel == jel) {
//                 elementList[iel] = 1;
//
//                 _particles[_markerOffset[iproc] + counter] = particles[ip];
//                 for(unsigned iList = 0; iList < _size; iList++) {
//                   if(printList[iList] == ip) {
//                     _printList[iList] = _markerOffset[iproc] + counter;
//                     break;
//                   }
//                 }
//                 counter++;
//               }
//             }
//           }
//         }
//         else {
//           someMarkersOutsideDomain = true;
//         }
//       }
//
//       if(someMarkersOutsideDomain == true) {
//         for(unsigned i = 0; i < _size; i++) {
//           unsigned jel;
//           jel = particles[i]->GetMarkerElement();
//           if(jel == UINT_MAX) {
//
//             _particles[_markerOffset[0] + counter] = particles[i];
//             for(unsigned iList = 0; iList < _size; iList++) {
//               if(printList[iList] == i) {
//                 _printList[iList] = _markerOffset[0] + counter;
//                 break;
//               }
//             }
//             counter++;
//           }
//         }
//         someMarkersOutsideDomain = false;
//       }
//     }
//
//
//
//
//     //END reorder markers also by element

    for (unsigned j = 0; j < _size; j++) {

      _particles[_printList[j]]->GetMarkerCoordinates(_line[j]);
    }
    _particles[_printList[0]]->GetMarkerCoordinates(_line[_size]);

    //delete[] elementList;
    std::vector < unsigned > ().swap(printList);


//     std::cout << "-----------------------------------------  UPDATE LINE  BY ELEMENTS-----------------------------------------" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc;
//       unsigned elem;
//       //std::vector <double> chiappe;
//       //chiappe = _particles[_printList[i]]->GetIprocMarkerCoordinates();
//       proc = _particles[i]->GetMarkerProc();
//       elem = _particles[i]->GetMarkerElement();
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " <<  " "  ;
// //       for(unsigned j = 0; j < _dim; j++) {
// //         std::cout << "x[" << j << "]=" << chiappe[j] << " " ;
// //       }
//       std::cout << std::endl;
//     }


  }

  void Line::AdvectionParallel(const unsigned &n, const double& T, const unsigned &order)
  {

    //BEGIN  Initialize the parameters for all processors

    double s;
    double PI = acos(-1.);
    std::vector <double> Fm(3, 0); // magnetic force initialization
    vector < unsigned > solVIndex(_dim);
    solVIndex[0] = _sol->GetIndex("U");    // get the position of "U" in the ml_sol object
    solVIndex[1] = _sol->GetIndex("V");    // get the position of "V" in the ml_sol object
    if (_dim == 3) solVIndex[2] = _sol->GetIndex("W");     // get the position of "V" in the ml_sol object
    unsigned solVType = _sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

    std::vector < double > phi;
    std::vector < std::vector<double > > V(2);
    std::map<unsigned, std::vector < std::vector < std::vector < double > > > > aV;
    std::map<unsigned, std::vector < std::vector < std::vector < std::vector < double > > > > > aX;
    double h = T / n;

    //END


    //BEGIN declare marker instances
    unsigned step;
    std::vector < double > x(_dim);
    std::vector < double > x0(_dim);
    std::vector < std::vector < double > > K(order);
    for (unsigned j = 0; j < order; j++) {
      K[j].resize(_dim);
    }
    //END

    unsigned integrationIsOverCounter = 0; // when integrationIsOverCounter = _size (which is the number of particles) it means all particles have been advected;
    // when the step of each marker is = n * order and/or marker element is UINT_MAX we increase by one such counter

    //BEGIN Numerical integration scheme

    for (unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {
      _particles[iMarker]->InitializeMarkerForAdvection(order);
    }
//     unsigned maxload = 0;
//     for(unsigned jproc=0; jproc<_nprocs;jproc++){
//       unsigned temp = (_markerOffset[jproc + 1] - _markerOffset[jproc]);
//       maxload = ( temp > maxload) ? temp : maxload;
//     }
//     maxload *= (n * order)/(_nprocs * _nprocs);

    while (integrationIsOverCounter != _size) {

      MyVector <unsigned> integrationIsOverCounterProc(1, 0);
      integrationIsOverCounterProc.stack();

      //BEGIN LOCAL ADVECTION INSIDE IPROC
      clock_t startTime = clock();
      unsigned counter = 0;
      for (unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {

        //std::cout << _printList[iMarker] <<" "<<std::flush;

        unsigned currentElem = _particles[iMarker]->GetMarkerElement();
        bool markerOutsideDomain = (currentElem != UINT_MAX) ? false : true;

        step = _particles[iMarker]->GetIprocMarkerStep();

        if (!markerOutsideDomain) {

          while (step < n * order) {

            x = _particles[iMarker]->GetIprocMarkerCoordinates();
            x0 = _particles[iMarker]->GetIprocMarkerOldCoordinates();
            K = _particles[iMarker]->GetIprocMarkerK();

            bool elementUpdate = (aX.find(currentElem) != aX.end()) ? false : true; //update if currentElem was never updated

            clock_t localTime = clock();
            _particles[iMarker]->GetMarkerS(n, order, s);
            _particles[iMarker]->FindLocalCoordinates(solVType, aX[currentElem], elementUpdate, _sol, s);
            _particles[iMarker]->updateVelocity(V, solVIndex, solVType, aV[currentElem], phi, elementUpdate, _sol); // we put pcElemUpdate instead of true but it wasn't running
            _time[3] += static_cast<double>((clock() - localTime)) / CLOCKS_PER_SEC;

            unsigned tstep = step / order;
            unsigned istep = step % order;

            if (istep == 0) {
              x0 = x;
              for (unsigned j = 0; j < order; j++) {
                K[j].assign(_dim, 0.);
              }
            }

            // double s = (tstep + _c[order - 1][istep]) / n;


            if (_sol->GetIfFSI()) {
              unsigned material = _sol->GetMesh()->GetElementMaterial(currentElem);
              MagneticForce(x, Fm, material, 1);
            }

//             for(unsigned l = 0; l < Fm.size(); l++) {
//               std::cout << "Fm[" << l << "]=" << Fm[l] << std::endl;
//             }

            for (unsigned k = 0; k < _dim; k++) {
              K[istep][k] = (s * V[0][k] + (1. - s) * V[1][k] + Fm[k]) * h;
            }

            counter++;

            step++;
            istep++;


            if (istep < order) {
              for (unsigned k = 0; k < _dim; k++) {
                x[k] = x0[k];
                for (unsigned j = 0; j < order; j++) {
                  x[k] +=  _a[order - 1][istep][j] * K[j][k];
                }
              }
            } else {
              for (unsigned i = 0; i < _dim; i++) {
                x[i] = x0[i];
                for (unsigned j = 0; j < order; j++) {
                  x[i] += _b[order - 1][j] * K[j][i];
                }
              }
            }

            _particles[iMarker]->SetIprocMarkerOldCoordinates(x0);
            _particles[iMarker]->SetIprocMarkerCoordinates(x);
            _particles[iMarker]->SetIprocMarkerK(K);

            _particles[iMarker]->SetIprocMarkerStep(step);
            _particles[iMarker]->GetMarkerS(n, order, s);

            //std::cout << step <<" "<< istep<< " AAAAAAAAAAAAA"<<std::endl<<std::flush;
            unsigned previousElem = currentElem;
            localTime = clock();
            _particles[iMarker]->GetElementSerial(previousElem, _sol, s);
            _time[4] += static_cast<double>((clock() - localTime)) / CLOCKS_PER_SEC;
            localTime = clock();
            // std::cout << step <<" "<< istep<< " BBBBBBBBBBBBB"<<std::endl<<std::flush;

            _particles[iMarker]->SetIprocMarkerPreviousElement(previousElem);

            currentElem = _particles[iMarker]->GetMarkerElement();
            unsigned mproc = _particles[iMarker]->GetMarkerProc(_sol);

            if (currentElem == UINT_MAX) { // the marker has been advected outise the domain
              markerOutsideDomain = true;
              step = UINT_MAX;
              _particles[iMarker]->SetIprocMarkerStep(step);
              break;
            } else if (_iproc != mproc) { // the marker has been advected outise the process
              break;
            }
          }
          if (step == n * order) {
            step = UINT_MAX;
            _particles[iMarker]->SetIprocMarkerStep(step);
          }
        } else { // the marker started outise the domain
          step = UINT_MAX;
          _particles[iMarker]->SetIprocMarkerStep(step);
        }

        if (step == UINT_MAX || markerOutsideDomain) {
          integrationIsOverCounterProc[_iproc] += 1;
        }
        //if(counter > maxload) break;
      }
      _time[5] += static_cast<double>((clock() - startTime)) / CLOCKS_PER_SEC;
      MPI_Barrier(PETSC_COMM_WORLD);
      _time[0] += static_cast<double>((clock() - startTime)) / CLOCKS_PER_SEC;
      startTime = clock();
      //END LOCAL ADVECTION INSIDE IPROC

      integrationIsOverCounter = 0;
      for (unsigned jproc = 0; jproc < _nprocs; jproc++) {
        integrationIsOverCounterProc.broadcast(jproc);
        integrationIsOverCounter += integrationIsOverCounterProc[jproc];
        integrationIsOverCounterProc.clearBroadcast();
      }

//       MPI_Barrier( PETSC_COMM_WORLD );
//       _time[1] += static_cast<double>((clock() - startTime)) / CLOCKS_PER_SEC;
//       startTime = clock();

      // std::cout << "DOPO integrationIsOverCounter = " << integrationIsOverCounter << std::endl;

      //BEGIN exchange on information

      //std::cout << " ----------------------------------PROCESSES EXCHANGE INFO ---------------------------------- " << std::endl;
      for (unsigned jproc = 0; jproc < _nprocs; jproc++) {
        for (unsigned iMarker = _markerOffset[jproc]; iMarker < _markerOffset[jproc + 1]; iMarker++) {
          unsigned elem =  _particles[iMarker]->GetMarkerElement();
          MPI_Bcast(& elem, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
          _particles[iMarker]->SetMarkerElement(elem);

          unsigned step =  _particles[iMarker]->GetIprocMarkerStep();
          MPI_Bcast(& step, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
          _particles[iMarker]->SetIprocMarkerStep(step);

          if (elem != UINT_MAX) { // if it is outside jproc, ACTUALLY IF WE ARE HERE IT COULD STILL BE IN JPROC but not outside the domain
            unsigned mproc = _particles[iMarker]->GetMarkerProc(_sol);
            _particles[iMarker]->SetMarkerProc(mproc);
            if (mproc != jproc) {
              unsigned prevElem = _particles[iMarker]->GetIprocMarkerPreviousElement();
              _particles[iMarker]->GetMarkerS(n, order, s);
              _particles[iMarker]->GetElement(prevElem, jproc, _sol, s);
              _particles[iMarker]->SetIprocMarkerPreviousElement(prevElem);
            }
            elem = _particles[iMarker]->GetMarkerElement();
            if (elem != UINT_MAX) { // if it is not outside the domain
              unsigned mproc = _particles[iMarker]->GetMarkerProc(_sol);
              if (mproc != jproc) {
                if (jproc == _iproc) {

                  unsigned step =  _particles[iMarker]->GetIprocMarkerStep();
                  MPI_Send(& step, 1, MPI_UNSIGNED, mproc, order + 1, PETSC_COMM_WORLD);

                  unsigned istep = step % order;
                  if (istep != 0) {
                    K = _particles[iMarker]->GetIprocMarkerK();
                    for (int i = 0; i < order; i++) {
                      MPI_Send(&K[i][0], _dim, MPI_DOUBLE, mproc, i , PETSC_COMM_WORLD);
                    }
                    x0 = _particles[iMarker]->GetIprocMarkerOldCoordinates();
                    MPI_Send(&x0[0], _dim, MPI_DOUBLE, mproc, order , PETSC_COMM_WORLD);
                  }
                  _particles[iMarker]->FreeXiX0andK();

                } else if (mproc == _iproc) {

                  _particles[iMarker]->InitializeX0andK(order);

                  unsigned step;
                  MPI_Recv(& step, 1, MPI_UNSIGNED, jproc, order + 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  _particles[iMarker]->SetIprocMarkerStep(step);

                  unsigned istep = step % order;
                  if (istep != 0) {
                    for (int i = 0; i < order; i++) {
                      MPI_Recv(&K[i][0], _dim, MPI_DOUBLE, jproc, i , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                      _particles[iMarker]->SetIprocMarkerK(K);
                    }
                    MPI_Recv(&x0[0], _dim, MPI_DOUBLE, jproc, order , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                    _particles[iMarker]->SetIprocMarkerOldCoordinates(x0);

                  }
                }
              }
            }
          }
          if (elem == UINT_MAX && jproc != 0) { // elem = UINT_MAX, but not yet in jproc = 0
            if (jproc == _iproc) {
              x = _particles[iMarker]->GetIprocMarkerCoordinates();
              MPI_Send(&x[0], _dim, MPI_DOUBLE, 0, 1 , PETSC_COMM_WORLD);
              _particles[iMarker]->FreeXiX0andK();
            } else if (_iproc == 0) {
              _particles[iMarker]->InitializeX();
              MPI_Recv(&x[0], _dim, MPI_DOUBLE, jproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
              _particles[iMarker]->SetIprocMarkerCoordinates(x);
            }
          }
        }
      }

      MPI_Barrier(PETSC_COMM_WORLD);
      _time[1] += static_cast<double>((clock() - startTime)) / CLOCKS_PER_SEC;
      startTime = clock();

      //std::cout << " ----------------------------------END PROCESSES EXCHANGE INFO ---------------------------------- " << std::endl;

      //END exchange of information


//       //BEGIN to be removed
//       std::vector<double> tr;
//       unsigned step;
//       for(unsigned j = 0; j < _size; j++) {
//         _particles[j]->GetMarkerCoordinates(tr);
//         step = _particles[j]->GetIprocMarkerStep();
//         for(unsigned i = 0; i < _dim; i++) {
//           std::cout << "x[ " << j << " ][ " << i << " ] = " << tr[i] ;
//         }
//         tr = _particles[j]->GetIprocMarkerOldCoordinates();
//         for(unsigned i = 0; i < _dim; i++) {
//           std::cout << "x0[ " << j << " ][ " << i << " ] = " << tr[i] ;
//         }
//
//
//         std::cout << " step = " << step << "currentElem = " << _particles[j]->GetMarkerElement() << " previousElem = " <<  _particles[j]->GetIprocMarkerPreviousElement() << std::endl;
//       }
//       //END to be removed



      UpdateLine();

      MPI_Barrier(PETSC_COMM_WORLD);
      _time[2] += static_cast<double>((clock() - startTime)) / CLOCKS_PER_SEC;
      startTime = clock();

//       //BEGIN to removed
//       for(unsigned j = 0; j < _size; j++) {
//         _particles[j]->GetMarkerCoordinates(tr);
//         step = _particles[j]->GetIprocMarkerStep();
//         for(unsigned i = 0; i < _dim; i++) {
//           std::cout << "x[ " << j << " ][ " << i << " ] = " << tr[i] ;
//         }
//         tr = _particles[j]->GetIprocMarkerOldCoordinates();
//         for(unsigned i = 0; i < _dim; i++) {
//           std::cout << "x0[ " << j << " ][ " << i << " ] = " << tr[i] ;
//         }
//         std::cout << " step = " << step << "currentElem = " << _particles[j]->GetMarkerElement() << " previousElem = " <<  _particles[j]->GetIprocMarkerPreviousElement() << std::endl;
//       }
//
//       //END to be remove

      // std::cout << " END integrationIsOverCounter = " << integrationIsOverCounter <<  std::endl;

    }

    //std::cout<<_time[0]<<" "<<_time[1]<<" "<<_time[2]<<" "<<_time[3]<<" "<<_time[4]<<" "<<_time[5]<<std::endl<<std::flush;
  }


  void Line::MagneticForce(const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned &material, unsigned forceType)
  {

    //case 0: infinitely long wire with a current I flowing modelled by the line identified by x and v
    //case 1: current loop of radius a, center x and for z-axis the line identified by x and v

    //BEGIN magnetic and electric parameters

    double PI = acos(-1.);
    double I = 1.e5; // electric current intensity
    double Msat = 1.e6;  //  magnetic saturation
    double  chi = 3.; //magnetic susceptibility
    double mu0 = 4 * PI * 1.e-7;  //magnetic permeability of the vacuum
    double H;
    std::vector <double> gradH(3);
    std::vector <double> gradHSquared(3);
    std::vector<double> vectorH(3);
    double H0 = Msat / chi;

    //END


    //BEGIN fluid viscosity

    double muf = (material == 2) ? 3.5 * 1.0e-3 : 1.0e100; // fluid viscosity

    //END


    //BEGIN geometric parameters

    double D = 500 * 1.e-9;       //diameter of the particle

    double a = 0.005; //radius of the circular current loop in cm

    std::vector <double> v(3);   //case 0: direction vector of the line that identifies the infinite wire
    //case 1: z axis of the current loop (symmetry axis)


//     aortic bifurcation wire
//     v[0] = 0.;
//     v[1] = 0.;
//     v[2] = 1.;


    
//     aortic bifurcation current loop
    v[0] = -1.;
    v[1] = 0.;
    v[2] = 0.;
    
    //TODO put the def of u here otherwise we will forget



    //bent tube no FSI wire
//     v[0] = 0.;
//     v[1] = 1.;
//     v[2] = 0.;

    std::vector <double> x(3);   //case 0: point that with v identifies the line of the wire
    //case 1: center of the current loop
//aortic bifurcation wire
//     x[0] = 0.015;
//     x[1] = 0.;
//     x[2] = 0.;


    //aortic bifurcation current loop
     x[0] = 0.02;
     x[1] = 0.;
     x[2] = 0.;

    
    //bent tube no FSI wire
//     x[0] = 9.;
//     x[1] = 0.;
//     x[2] = 3.;


    //END


    //BEGIN extraction of the coordinates of the particle

    std::vector <double> xM(3);
    xM[0] = xMarker[0];
    xM[1] = xMarker[1];
    xM[2] = (xMarker.size() == 3) ? xMarker[2] : 0. ;

    //END


    //BEGIN evaluate H

    switch (forceType) { 

      case 0: {  // infinite wire

          double Gamma;
          double Omega;
          std::vector<double> gradOmega(3);

          Gamma = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
          Omega = (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) +
                  (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) +
                  (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1])) * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1])) ;

          gradOmega[0] = 2 * v[2] * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) + 2 * v[1] * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]));

          gradOmega[1] = 2 * v[2] * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) - 2 * v[0] * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]));

          gradOmega[2] = - 2 * v[1] * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) - 2 * v[0] * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2]));

          H = (I / (2 * PI)) * Gamma * (1 / sqrt(Omega));

          gradHSquared[0] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[0];
          gradHSquared[1] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[1];
          gradHSquared[2] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[2];

          gradH[0] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[0]) / sqrt(Omega * Omega * Omega);
          gradH[1] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[1]) / sqrt(Omega * Omega * Omega);
          gradH[2] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[2]) / sqrt(Omega * Omega * Omega);

        }
        break;

      case 1: {  //current loop

          for (unsigned i = 0; i < 3 ; i++) {
            xM[i] -= x[i];
          }

          std::vector < double > vectorHp(3);
          std::vector < std::vector < double > > jacobianVectorHp(3);
          for (unsigned i = 0; i < 3; i++) {
            jacobianVectorHp[i].resize(3);
          }
          std::vector < std::vector <double> > jacobianVectorH(3);
          for (unsigned i = 0; i < 3; i++) {
            jacobianVectorH[i].resize(3);
          }

          double rhoSquared = xM[0] * xM[0] + xM[1] * xM[1];
          double rSquared = rhoSquared + xM[2] * xM[2];
          double alphaSquared = a * a + rSquared - 2 * a * sqrt(rhoSquared);
          double betaSquared = a * a + rSquared + 2 * a * sqrt(rhoSquared);
          double kSquared = 1 - (alphaSquared / betaSquared);
          double gamma = xM[0] * xM[0] - xM[1] * xM[1];
          double C = (I * mu0) / PI;


          double v2 = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] );
          v[0] /= v2;
          v[1] /= v2;
          v[2] /= v2;

          std::vector <double> u(3); // will be the x'-axis
          std::vector <double> w(3); // will be the y'-axis

	  u[0] = 1.;
	  u[1] = 1.;
	  u[2] = 0.;
	  	  
	  
//           unsigned imax = 0;
//           double vmax = fabs(v[0]);
//           for (unsigned i = 1; i < 3; i++) {
//             if ( fabs(v[i]) > vmax ) {
//               imax = i;
//               vmax = fabs(v[i]);
//             }
//           }

//           u[0] = 1. - (imax == 0) * (1 + (v[1] + v[2]) / v[0] );
//           u[1] = 1. - (imax == 1) * (1 + (v[2] + v[0]) / v[1] );
//           u[2] = 1. - (imax == 2) * (1 + (v[0] + v[1]) / v[2] );

          double u2 = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );
          u[0] /= u2;
          u[1] /= u2;
          u[2] /= u2;

          w[0] = v[1] * u[2] - v[2] * u[1];
          w[1] = v[2] * u[0] - v[0] * u[2];
          w[2] = v[0] * u[1] - v[1] * u[0];

          double w2 = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2] );
          if (fabs(w2 - 1.) > 1.0e-12) {
            std::cout << "ERRORE AAAAAAAAAAAAAAAAAAA" << std::endl;
          }

          std::vector< std::vector <double> > R(3);
          for (unsigned i = 0; i < 3; i++) {
            R[i].resize(3);
            R[i][0] = u[i];
            R[i][1] = w[i];
            R[i][2] = v[i];
          }

          vectorHp[0] = (1 / mu0) * ((C * xM[0] * xM[2]) / (2 * alphaSquared * sqrt(betaSquared) * rhoSquared)) * ((a * a + rSquared) * ellint_2(kSquared) - alphaSquared * ellint_1(kSquared));
          vectorHp[1] = (xM[1] / xM[0]) * vectorH[0];
          vectorHp[2] = (1 / mu0) * (C / (2 * alphaSquared * sqrt(betaSquared))) * ((a * a - rSquared) * ellint_2(kSquared) + alphaSquared * ellint_1(kSquared));


          jacobianVectorHp[0][0] = ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared)) * (
                                     (a * a * a * a * (- gamma * (3 * xM[2] * xM[2] + a * a) + rhoSquared * (8 * xM[0] * xM[0] - xM[1] * xM[1])) -
                                      a * a * (rhoSquared * rhoSquared * (5 * xM[0] * xM[0] + xM[1] * xM[1]) -
                                               2 * rhoSquared * xM[2] * xM[2] * (2 * xM[0] * xM[0] + xM[1] * xM[1]) + 3 * xM[2] * xM[2] * xM[2] * xM[2] * xM[1]) -
                                      rSquared * rSquared * (2 * xM[0] * xM[0] * xM[0] * xM[0] + gamma * (xM[1] * xM[1] + xM[2] * xM[2])))  * ellint_2(kSquared) +
                                     (a * a * (gamma * (a * a + 2 * xM[2] * xM[2]) - rhoSquared * (3 * xM[0] * xM[0] - 2 * xM[1] * xM[1])) +
                                      rSquared * (2 * xM[0] * xM[0] * xM[0] * xM[0] + gamma * (xM[1] * xM[1] + xM[2] * xM[2]))) * alphaSquared * ellint_1(kSquared));


          jacobianVectorHp[0][1] = ((C * xM[0] * xM[1] * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared)) * (
                                     (3 * a * a * a * a * (3 * rhoSquared - 2 * xM[2] * xM[2]) - rSquared * rSquared * (2 * rSquared + rhoSquared) - 2 * a * a * a * a * a * a  -
                                      2 * a * a * (2 * rhoSquared * rhoSquared - rhoSquared * xM[2] * xM[2] + 3 * xM[2] * xM[2] * xM[2] * xM[2])) * ellint_2(kSquared) +
                                     (rSquared * (2 * rSquared + rhoSquared) - a * a * (5 * rhoSquared - 4 * xM[2] * xM[2]) + 2 * a * a * a * a) * alphaSquared * ellint_1(kSquared));


          jacobianVectorHp[0][2] = ((C * xM[0]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared)) * (
                                     ((rhoSquared - a * a) * (rhoSquared - a * a) * (rhoSquared + a * a) + 2 * xM[2] * xM[2] * (a * a * a * a - 6 * a * a * rhoSquared + rhoSquared * rhoSquared) +
                                      xM[2] * xM[2] * xM[2] * xM[2] * (a * a + rhoSquared)) *  ellint_2(kSquared) -
                                     ((rhoSquared - a * a) * (rhoSquared - a * a) + xM[2] * xM[2] * (rhoSquared + a * a)) * alphaSquared *  ellint_1(kSquared));


          jacobianVectorHp[1][0] = jacobianVectorHp[0][1];


          jacobianVectorHp[1][1] = ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared)) * (
                                     (a * a * a * a * (gamma * (3 * xM[2] * xM[2] + a * a) + rhoSquared * (8 * xM[1] * xM[1] - xM[0] * xM[0])) -
                                      a * a * (rhoSquared * rhoSquared * (5 * xM[1] * xM[1] + xM[0] * xM[0]) - 2 * rhoSquared * xM[2] * xM[2] * (2 * xM[1] * xM[1] + xM[0] * xM[0]) -
                                               3 * xM[2] * xM[2] * xM[2] * xM[2] * xM[1])  - rSquared * rSquared * (2 * xM[1] * xM[1] * xM[1] * xM[1] - gamma * (xM[0] * xM[0] + xM[2] * xM[2]))) * ellint_2(kSquared) +
                                     (a * a * (- gamma * (a * a + 2 * xM[2] * xM[2])  - rhoSquared * (3 * xM[1] * xM[1] - 2 * xM[0] * xM[0])) + rSquared * (2 * xM[1] * xM[1] * xM[1] * xM[1] -
                                         gamma * (xM[0] * xM[0] + xM[2] * xM[2]))) * alphaSquared * ellint_1(kSquared));


          jacobianVectorHp[1][2] = (xM[1] / xM[0]) * jacobianVectorHp[0][2];


          jacobianVectorHp[2][0] = jacobianVectorHp[0][2];


          jacobianVectorHp[2][1] = jacobianVectorHp[1][2];


          jacobianVectorHp[2][2] = ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared))) * (
                                     (6 * a * a * (rhoSquared - xM[2] * xM[2]) - 7 * a * a * a * a + rSquared * rSquared) * ellint_2(kSquared) +
                                     (a * a - rSquared) * alphaSquared * ellint_1(kSquared));


          for (unsigned i = 0; i < 3; i++) {
            vectorH[i] = 0.;
            for (unsigned j = 0; j < 3; j++) {
              vectorH[i] += R[i][j] * vectorHp[j];
              jacobianVectorH[i][j] = 0.;
	      //std::cout << R[i][j] <<" ";
	      //std::cout << jacobianVectorHp[i][j] <<" ";
              for (unsigned k = 0; k < 3; k++) {
                jacobianVectorH[i][j] += R[i][k] * jacobianVectorHp[k][j];
              }
            }
          }


          for (unsigned i = 0; i < 3; i++) {
            for (unsigned j = 0; j < 3; j++) {
              jacobianVectorHp[i][j] = 0. ;
              for (unsigned k = 0; k < 3; k++) {
                jacobianVectorHp[i][j] += jacobianVectorH[i][k] * R[j][k];
              }
            }
          }

          for (unsigned i = 0; i < 3; i++) {
            for (unsigned j = 0; j < 3; j++) {
	      jacobianVectorH[i][j] = jacobianVectorHp[i][j];
	    }
          }


          H = sqrt(vectorH[0] * vectorH[0] + vectorH[1] * vectorH[1] * vectorH[2] * vectorH[2]);


          gradHSquared[0] = 2 * vectorH[0] * jacobianVectorH[0][0] + 2 * vectorH[1] * jacobianVectorH[1][0] + 2 * vectorH[2] * jacobianVectorH[2][0];
          gradHSquared[1] = 2 * vectorH[0] * jacobianVectorH[0][1] + 2 * vectorH[1] * jacobianVectorH[1][1] + 2 * vectorH[2] * jacobianVectorH[2][1];
          gradHSquared[2] = 2 * vectorH[0] * jacobianVectorH[0][2] + 2 * vectorH[1] * jacobianVectorH[1][2] + 2 * vectorH[2] * jacobianVectorH[2][2];

          gradH[0] = 0.5 * (1 / H) * gradHSquared[0];
          gradH[1] = 0.5 * (1 / H) * gradHSquared[1];
          gradH[2] = 0.5 * (1 / H) * gradHSquared[2];

	  
	  
        }
        
        //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
        
        break;

    }
    //END valuate H


    //BEGIN evaluate Fm

    if (H < H0) {


      double C1 = (PI * D * D * D * mu0 * chi) / 12;

      Fm[0] = C1 * gradHSquared[0];
      Fm[1] = C1 * gradHSquared[1];
      Fm[2] = C1 * gradHSquared[2];
    }

    else {

      double C2 = (PI * D * D * D * mu0 * Msat) / 6;

      Fm[0] = C2 * gradH[0];
      Fm[1] = C2 * gradH[1];
      Fm[2] = C2 * gradH[2];

    }


    for (unsigned i = 0 ; i < Fm.size(); i++) {
      Fm[i] = Fm[i] / (3 * PI * D * muf) ;
    }

    //END


    //BEGIN cheating to have attractive force

    for (unsigned i = 0 ; i < Fm.size(); i++) {
      Fm[i] = - Fm[i] ;
      //std::cout << Fm[i] << " " <<std::flush;
    }

    //END cheating

    
    
    
  }


}





