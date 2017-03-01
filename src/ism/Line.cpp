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
#include <math.h>

#include "PolynomialBases.hpp"

namespace femus {

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
             Solution *sol, const unsigned & solType) {

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

    for(unsigned j = 0; j < _size; j++) {
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
    for(unsigned iproc = 0; iproc < _nprocs; iproc++) {
      _markerOffset[iproc] = counter;
      //unsigned  offsetCounter = counter;
      for(unsigned j = 0; j < _size; j++) {
        unsigned markerProc = particles[j]->GetMarkerProc(_sol);
        if(markerProc == iproc) {
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

    for(unsigned j = 0; j < _size; j++) {
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

  Line::~Line() {
    for(unsigned j = 0; j < _size; j++) {
      //std::cout << j << " " << _particles[j] << std::endl<< std::flush;
      delete _particles[j];
    }
  }

  void Line::UpdateLine() {

    std::vector < Marker*> particles(_size);
    std::vector < unsigned> printList(_size);

    for(unsigned j = 0; j < _size; j++) {
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
    for(unsigned iproc = 0; iproc < _nprocs; iproc++) {
      _markerOffset[iproc] = counter;
      for(unsigned j = 0; j < _size; j++) {
        unsigned markerProc = particles[j]->GetMarkerProc(_sol);
        if(markerProc == iproc) {
          _particles[counter] = particles[j];
          for(unsigned iList = 0; iList < _size; iList++) {
            if(printList[iList] == j) {
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

    for(unsigned j = 0; j < _size; j++) {

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

  void Line::AdvectionParallel(const unsigned &n, const double& T, const unsigned &order) {

    //BEGIN  Initialize the parameters for all processors

    double s;
    double PI = acos(-1.);
    std::vector <double> Fm(3,0); // magnetic force initialization
    vector < unsigned > solVIndex(_dim);
    solVIndex[0] = _sol->GetIndex("U");    // get the position of "U" in the ml_sol object
    solVIndex[1] = _sol->GetIndex("V");    // get the position of "V" in the ml_sol object
    if(_dim == 3) solVIndex[2] = _sol->GetIndex("W");      // get the position of "V" in the ml_sol object
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
    for(unsigned j = 0; j < order; j++) {
      K[j].resize(_dim);
    }
    //END

    unsigned integrationIsOverCounter = 0; // when integrationIsOverCounter = _size (which is the number of particles) it means all particles have been advected;
    // when the step of each marker is = n * order and/or marker element is UINT_MAX we increase by one such counter

    //BEGIN Numerical integration scheme

    for(unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {
      _particles[iMarker]->InitializeMarkerForAdvection(order);
    }
//     unsigned maxload = 0;
//     for(unsigned jproc=0; jproc<_nprocs;jproc++){
//       unsigned temp = (_markerOffset[jproc + 1] - _markerOffset[jproc]);
//       maxload = ( temp > maxload) ? temp : maxload;
//     }
//     maxload *= (n * order)/(_nprocs * _nprocs);

    while(integrationIsOverCounter != _size) {

      MyVector <unsigned> integrationIsOverCounterProc(1, 0);
      integrationIsOverCounterProc.stack();

      //BEGIN LOCAL ADVECTION INSIDE IPROC
      clock_t startTime = clock();
      unsigned counter = 0;
      for(unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {

        //std::cout << _printList[iMarker] <<" "<<std::flush;

        unsigned currentElem = _particles[iMarker]->GetMarkerElement();
        bool markerOutsideDomain = (currentElem != UINT_MAX) ? false : true;

        step = _particles[iMarker]->GetIprocMarkerStep();

        if(!markerOutsideDomain) {

          while(step < n * order) {

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

            if(istep == 0) {
              x0 = x;
              for(unsigned j = 0; j < order; j++) {
                K[j].assign(_dim, 0.);
              }
            }

            // double s = (tstep + _c[order - 1][istep]) / n;
            
            
            MagneticForceWire(x,Fm);
	    

	    for(unsigned l=0; l<Fm.size(); l++){
	      std::cout << "Fm[" << l << "]=" << Fm[l] << std::endl;
	    }
	    
	    for(unsigned k = 0; k < _dim; k++) {
              K[istep][k] = (s * V[0][k] + (1. - s) * V[1][k] + Fm[k]) * h; 
            }

            counter++;

            step++;
            istep++;


            if(istep < order) {
              for(unsigned k = 0; k < _dim; k++) {
                x[k] = x0[k];
                for(unsigned j = 0; j < order; j++) {
                  x[k] +=  _a[order - 1][istep][j] * K[j][k];
                }
              }
            }
            else {
              for(unsigned i = 0; i < _dim; i++) {
                x[i] = x0[i];
                for(unsigned j = 0; j < order; j++) {
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

            if(currentElem == UINT_MAX) { // the marker has been advected outise the domain
              markerOutsideDomain = true;
              step = UINT_MAX;
              _particles[iMarker]->SetIprocMarkerStep(step);
              break;
            }
            else if(_iproc != mproc) { // the marker has been advected outise the process
              break;
            }
          }
          if(step == n * order) {
            step = UINT_MAX;
            _particles[iMarker]->SetIprocMarkerStep(step);
          }
        }
        else { // the marker started outise the domain
          step = UINT_MAX;
          _particles[iMarker]->SetIprocMarkerStep(step);
        }

        if(step == UINT_MAX || markerOutsideDomain) {
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
      for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
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
      for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
        for(unsigned iMarker = _markerOffset[jproc]; iMarker < _markerOffset[jproc + 1]; iMarker++) {
          unsigned elem =  _particles[iMarker]->GetMarkerElement();
          MPI_Bcast(& elem, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
          _particles[iMarker]->SetMarkerElement(elem);

          unsigned step =  _particles[iMarker]->GetIprocMarkerStep();
          MPI_Bcast(& step, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
          _particles[iMarker]->SetIprocMarkerStep(step);

          if(elem != UINT_MAX) {  // if it is outside jproc, ACTUALLY IF WE ARE HERE IT COULD STILL BE IN JPROC but not outside the domain
            unsigned mproc = _particles[iMarker]->GetMarkerProc(_sol);
            _particles[iMarker]->SetMarkerProc(mproc);
            if(mproc != jproc) {
              unsigned prevElem = _particles[iMarker]->GetIprocMarkerPreviousElement();
              _particles[iMarker]->GetMarkerS(n, order, s);
              _particles[iMarker]->GetElement(prevElem, jproc, _sol, s);
              _particles[iMarker]->SetIprocMarkerPreviousElement(prevElem);
            }
            elem = _particles[iMarker]->GetMarkerElement();
            if(elem != UINT_MAX) {  // if it is not outside the domain
              unsigned mproc = _particles[iMarker]->GetMarkerProc(_sol);
              if(mproc != jproc) {
                if(jproc == _iproc) {

                  unsigned step =  _particles[iMarker]->GetIprocMarkerStep();
                  MPI_Send(& step, 1, MPI_UNSIGNED, mproc, order + 1, PETSC_COMM_WORLD);

                  unsigned istep = step % order;
                  if(istep != 0) {
                    K = _particles[iMarker]->GetIprocMarkerK();
                    for(int i = 0; i < order; i++) {
                      MPI_Send(&K[i][0], _dim, MPI_DOUBLE, mproc, i , PETSC_COMM_WORLD);
                    }
                    x0 = _particles[iMarker]->GetIprocMarkerOldCoordinates();
                    MPI_Send(&x0[0], _dim, MPI_DOUBLE, mproc, order , PETSC_COMM_WORLD);
                  }
                  _particles[iMarker]->FreeXiX0andK();

                }
                else if(mproc == _iproc) {

                  _particles[iMarker]->InitializeX0andK(order);

                  unsigned step;
                  MPI_Recv(& step, 1, MPI_UNSIGNED, jproc, order + 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  _particles[iMarker]->SetIprocMarkerStep(step);

                  unsigned istep = step % order;
                  if(istep != 0) {
                    for(int i = 0; i < order; i++) {
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
          if(elem == UINT_MAX && jproc != 0) { // elem = UINT_MAX, but not yet in jproc = 0
            if(jproc == _iproc) {
              x = _particles[iMarker]->GetIprocMarkerCoordinates();
              MPI_Send(&x[0], _dim, MPI_DOUBLE, 0, 1 , PETSC_COMM_WORLD);
              _particles[iMarker]->FreeXiX0andK();
            }
            else if(_iproc == 0) {
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


  void Line::MagneticForceWire(const std::vector <double> & xMarker, std::vector <double> &Fm) {

    double PI = acos(-1.);
    
    double mu0 = 4 * PI * 1.e-7;  //magnetic permeability of the vacuum
    
    double muf = 3.5 * 1.0e-3; // fluid viscosity
    
    double D = 500 * 1.e-9;       //diameter of the particle
    
    std::vector <double> v(3);   //direction vector of the line that identifies the infinite wire
    
//     aortic bifurcation
    v[0] = 0.;
    v[1] = 0.;
    v[2] = 1.;
    
    //bent tube no FSI
//     v[0] = 0.;
//     v[1] = 1.;
//     v[2] = 0.;
    
    std::vector <double> x(3);   //point that with v identifies the line
    
//aortic bifurcation
    x[0] = -0.016;
    x[1] = 0.06;
    x[2] = 0.;

//bent tube no FSI
//     x[0] = 9.;
//     x[1] = 0.;
//     x[2] = 3.;
    
    double I = 2.e5; // electric current intensity
    double Msat = 1.e6;  //  magnetic saturation
    double  chi = 3.; //magnetic susceptibility

    std::vector <double> xM(3);
    xM[0] = xMarker[0];
    xM[1] = xMarker[1];
    xM[2] = (xMarker.size() == 3) ? xMarker[2] : 0. ;

    double Gamma = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    double Omega = v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2]) * v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2]) +
                   v[0] * (x[2] - xM[2]) - v[2] * (x[0] - xM[0]) * v[0] * (x[2] - xM[2]) - v[2] * (x[0] - xM[0]) +
                   v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]) * v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]) ;

    std::vector<double> gradOmega(3);

    gradOmega[0] = 2 * v[2] * (v[0] * (x[2] - xM[2]) - v[2] * (x[0] - xM[0])) - 2 * v[1] * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]));

    gradOmega[1] = -2 * v[2] * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) + 2 * v[0] * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]));
    
    gradOmega[2] = 2 * v[1] * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) - 2 * v[0] * (v[0] * (x[2] - xM[2]) - v[2] * (x[0] - xM[0]));
    
    double H = (I / PI) * Gamma * (1 / sqrt(Omega)); //absolute value of the magnetic field at x = xMarker

    double H0 = Msat / chi;

    if(H < H0) {

      double C1 = (PI * D * D *D * mu0 * chi) / 12;

      Fm[0] = -C1 * I *  I * Gamma * Gamma * gradOmega[0] /(Omega*Omega*4*PI*PI);
      Fm[1] = -C1 * I *  I * Gamma * Gamma * gradOmega[1] /(Omega*Omega*4*PI*PI);
      Fm[2] = -C1 * I *  I * Gamma * Gamma * gradOmega[2] /(Omega*Omega*4*PI*PI);

    }
    
    else{
      
      double C2 = (PI * D * D *D * mu0 * Msat) / 6;
      
      Fm[0] = - C2 * I *  Gamma * gradOmega[0] / sqrt(Omega*Omega*Omega*2*PI);
      Fm[1] = - C2 * I *  Gamma * gradOmega[1] / sqrt(Omega*Omega*Omega*2*PI); 
      Fm[2] = - C2 * I *  Gamma * gradOmega[2] / sqrt(Omega*Omega*Omega*2*PI);
      
    }
    
    for(unsigned i = 0 ; i<Fm.size(); i++){
      Fm[i] = Fm[i] / (3 * PI * D *muf) ;
    }
    

  }


}




