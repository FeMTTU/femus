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
             Mesh *mesh, const unsigned & solType) {

    _size = x.size();

    std::vector < Marker*> particles(_size);

    _dim = mesh->GetDimension();

    _markerOffset.resize(_nprocs + 1);
    _markerOffset[_nprocs] = _size;

    _particles.resize(_size);
    _printList.resize(_size);

    for(unsigned j = 0; j < _size; j++) {
      particles[j] = new Marker(x[j], markerType[j], mesh, solType, true);
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
        unsigned markerProc = particles[j]->GetMarkerProc();
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

    //BEGIN reorder markers also by element

    for(unsigned j = 0; j < _size; j++) {
      particles[j] = _particles[j];
    }

    unsigned* elementList = new unsigned [mesh->GetNumberOfElements()];
    std::vector < unsigned> printList(_size);

    //flags to see if we already ordered by that element, 0 = not considered yet
    for(unsigned iFlag = 0; iFlag < mesh->GetNumberOfElements(); iFlag++) {
      elementList[iFlag] = 0;
    }

    printList = _printList;

    for(unsigned iproc = 0; iproc < _nprocs; iproc++) {

      counter = 0;

      for(unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {

        unsigned jel;
        particles[jp]->GetMarkerElementLine(jel);

        if(elementList[jel] == 0) {

          elementList[jel] = 1;

          _particles[_markerOffset[iproc] + counter] = particles[jp];
          for(unsigned iList = 0; iList < _size; iList++) {
            if(printList[iList] == jp) {
              _printList[iList] = _markerOffset[iproc] + counter;
              break;
            }
          }
          counter++;


          for(unsigned ip = _markerOffset[iproc]; ip < _markerOffset[iproc + 1]; ip++) {
            unsigned iel;
            particles[ip]->GetMarkerElementLine(iel);
            if(ip != jp && iel == jel) {
              //std::cout << " jel =" << jel << " , " << "ip = " << ip << " , " << "jp = " << jp <<  std::endl;
              elementList[iel] = 1;
              _particles[_markerOffset[iproc] + counter] = particles[ip];
              for(unsigned iList = 0; iList < _size; iList++) {
                if(printList[iList] == ip) {
                  _printList[iList] = _markerOffset[iproc] + counter;
                  break;
                }
              }
              counter++;
            }
          }
        }
      }
    }

    //END reorder markers also by element

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


    _line.resize(1);
    _line[0].resize(_size + 1);

    for(unsigned j = 0; j < _size; j++) {
      _particles[_printList[j]]->GetMarkerCoordinates(_line[0][j]);
    }
    _particles[_printList[0]]->GetMarkerCoordinates(_line[0][_size]);

    delete[] elementList;
    std::vector < unsigned > ().swap(printList);

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
        unsigned markerProc = particles[j]->GetMarkerProc();
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
//       _particles[i]->GetMarkerProcLine(proc);
//       _particles[i]->GetMarkerElementLine(elem);
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " << std::endl;
//     }



    //BEGIN reorder markers also by element

    for(unsigned j = 0; j < _size; j++) {
      particles[j] = _particles[j];
      //printList[j] = _printList[j];
    }

    printList = _printList;

    unsigned meshElements;
    particles[0]->GetNumberOfMeshElements(meshElements);
    unsigned* elementList = new unsigned [meshElements];
    bool someMarkersOutsideDomain = false;

    //flags to see if we already ordered by that element, 0 = not considered yet
    for(unsigned iFlag = 0; iFlag < meshElements; iFlag++) {
      elementList[iFlag] = 0;
    }

    for(unsigned iproc = 0; iproc < _nprocs; iproc++) {

      counter = 0;

      for(unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {

        unsigned jel;
        particles[jp]->GetMarkerElementLine(jel);

        if(jel != UINT_MAX) {

          if(elementList[jel] == 0) {

            elementList[jel] = 1;

            _particles[_markerOffset[iproc] + counter] = particles[jp];
            for(unsigned iList = 0; iList < _size; iList++) {
              if(printList[iList] == jp) {
                _printList[iList] = _markerOffset[iproc] + counter;
                break;
              }
            }
            counter++;


            for(unsigned ip = _markerOffset[iproc]; ip < _markerOffset[iproc + 1]; ip++) {
              unsigned iel;
              particles[ip]->GetMarkerElementLine(iel);
              if(ip != jp && iel == jel) {
                elementList[iel] = 1;
                _particles[_markerOffset[iproc] + counter] = particles[ip];
                for(unsigned iList = 0; iList < _size; iList++) {
                  if(printList[iList] == ip) {
                    _printList[iList] = _markerOffset[iproc] + counter;
                    break;
                  }
                }
                counter++;
              }
            }
          }
        }
        else {
          someMarkersOutsideDomain = true;
        }
      }
      if(someMarkersOutsideDomain == true) {
        for(unsigned i = 0; i < _size; i++) {
          unsigned jel;
          particles[i]->GetMarkerElementLine(jel);
          if(jel == UINT_MAX) {
            _particles[_markerOffset[_nprocs - 1] + counter] = particles[i];
            for(unsigned iList = 0; iList < _size; iList++) {
              if(printList[iList] == i) {
                _printList[iList] = _markerOffset[_nprocs - 1] + counter;
                break;
              }
            }
            counter++;
          }
        }
      }
    }




    //END reorder markers also by element

    for(unsigned j = 0; j < _size; j++) {

      _particles[_printList[j]]->GetMarkerCoordinates(_line[0][j]);
    }
    _particles[_printList[0]]->GetMarkerCoordinates(_line[0][_size]);

    delete[] elementList;
    std::vector < unsigned > ().swap(printList);


//     std::cout << "-----------------------------------------  UPDATE LINE  BY ELEMENTS-----------------------------------------" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc;
//       unsigned elem;
//       std::vector <double> chiappe;
//       chiappe = _particles[_printList[i]]->GetIprocMarkerCoordinates();
//       proc = _particles[_printList[i]]->GetMarkerProc();
//       elem = _particles[_printList[i]]->GetMarkerElement();
//       std::cout << "Particle: " << _printList[i] << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " << "_printList[ " << i << "] = " << _printList[i] << " "  ;
//       for(unsigned j = 0; j < _dim; j++) {
//         std::cout << "x[" << j << "]=" << chiappe[j] << " " ;
//       }
//       std::cout << std::endl;
//     }


  }


  void Line::AdvectionParallel(Solution* sol, const unsigned &n, const double& T, const unsigned &order) {

    //BEGIN  Initialize the parameters for all processors

    vector < unsigned > solVIndex(_dim);
    solVIndex[0] = sol->GetIndex("U");    // get the position of "U" in the ml_sol object
    solVIndex[1] = sol->GetIndex("V");    // get the position of "V" in the ml_sol object
    if(_dim == 3) solVIndex[2] = sol->GetIndex("W");       // get the position of "V" in the ml_sol object
    unsigned solVType = sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

    std::vector < double > phi;
    std::vector < std::vector<double > > V(2);
    std::vector < std::vector < std::vector < double > > > aV;
    std::vector < std::vector < std::vector < double > > > aX;
    double h = T / n;

    //END


    //std::vector<unsigned> previousElem(_size, 0); //sbagliato

    //BEGIN declare marker instances
    unsigned step;
    std::vector < double > x;
    std::vector < double > x0;
    std::vector < std::vector < double > > K;
    //END

    unsigned integrationIsOverCounter = 0; // when integrationIsOverCounter = _size (which is the number of particles) it means all particles have been advected;
    // when the step of each marker is = n * order and/or marker element is UINT_MAX we increase by one such counter



    //BEGIN Numerical integration scheme

    for(unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {
      _particles[iMarker]->InitializeMarkerForAdvection(order);
    }

    while(integrationIsOverCounter != _size) {

      std::vector <unsigned> integrationIsOverCounterProc(_nprocs, 0);

      // std::cout << " BEGIN integrationIsOverCounter = " << integrationIsOverCounter << std::endl;

      unsigned initializedElem = UINT_MAX;
      unsigned currentElem;


      //BEGIN LOCAL ADVECTION INSIDE IPROC
      for(unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {



        //BEGIN extraction of the marker instances
        x = _particles[iMarker]->GetIprocMarkerCoordinates();
        x0 = _particles[iMarker]->GetIprocMarkerOldCoordinates();
        K = _particles[iMarker]->GetIprocMarkerK();
        step = _particles[iMarker]->GetIprocMarkerStep();
        currentElem = _particles[iMarker]->GetMarkerElement();
        //END


        for(unsigned i = 0; i < _dim; i++) {
          std::cout << "Coordinates as we begin x[ " << i << " ] = " << x[_printList[i]] << std::endl;
        }


        bool markerOutsideDomain = (currentElem != UINT_MAX) ? false : true;

        if(markerOutsideDomain == false) {

          std::cout << currentElem << " " << initializedElem << std::endl;

          bool pcElemUpdate = (initializedElem == currentElem) ? false : true; //update only if the marker is in a different element
          std::cout << " PRIMA DI FIND LOCAL COORDINATES pcElemUpdate =" << pcElemUpdate << std::endl;
          _particles[iMarker]->FindLocalCoordinates(solVType, aX, pcElemUpdate);
          std::cout << " PRIMA DI UPDATE VELOCITY pcElemUpdate =" << pcElemUpdate << std::endl;
          initializedElem = currentElem;

          while(step < n * order) {

            unsigned tstep = step / order;
            unsigned istep = step % order;

            if(istep == 0) {
              x0 = x;
              K.resize(order);
              for(unsigned j = 0; j < order; j++) {
                K[j].assign(_dim, 0.);
              }
            }

            
            _particles[iMarker]->updateVelocity(V, sol, solVIndex, solVType, aV, phi, true); //switch back to pcElemUpdate al posto di true
	    std::cout << " DOPO UPDATE VELOCITY pcElemUpdate =" << pcElemUpdate << std::endl;
	    pcElemUpdate = false;
            
            double s = (tstep + _c[order - 1][istep]) / n;

            for(unsigned k = 0; k < _dim; k++) {
              K[istep][k] = (s * V[0][k] + (1. - s) * V[1][k]) * h;
            }

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

            else if(istep == order) {

              for(unsigned i = 0; i < _dim; i++) {
                x[i] = x0[i];
                for(unsigned j = 0; j < order; j++) {
                  x[i] += _b[order - 1][j] * K[j][i];
                }
              }
            }

            _particles[iMarker]->SetIprocMarkerOldCoordinates(x0);
            _particles[iMarker]->SetIprocMarkerCoordinates(x);
            _particles[iMarker]->SetIprocMarkerStep(step);
            _particles[iMarker]->SetIprocMarkerK(K);

            _particles[iMarker]->GetElementSerial(currentElem);
            currentElem = _particles[iMarker]->GetMarkerElement(); // probablilmente non c'e' ne bisogno

            unsigned markProc = _particles[iMarker]->GetMarkerProc();

            _particles[iMarker]->SetIprocMarkerPreviousElement(initializedElem);

            if(currentElem == UINT_MAX) {  //out of the domain
              // std::cout << " OUT OF DOMAIN " <<std::endl;
              markerOutsideDomain = true;
              break;
            }
            else if(initializedElem != currentElem && _iproc != markProc) {    //different element different process
              // std::cout << " DIFFERENT ELEMENT DIFFERENT PROCESS " <<std::endl;
              _particles[iMarker]->SetIprocMarkerPreviousElement(initializedElem);
              break;
            }
            else if(initializedElem != currentElem) {    //different element same process
              //std::cout << " DIFFERENT ELEM SAME PROCESS " <<std::endl;
              break;
            }
            else {   //same element same process
              // std::cout << " SAME ELEM SAME PROCESS " << std::endl;
              _particles[iMarker]->FindLocalCoordinates(solVType, aX, false); //inverse mapping to continue
            }
          }

          if(step == n * order) {
            integrationIsOverCounterProc[_iproc] += 1;
            step = UINT_MAX;
            _particles[iMarker]->SetIprocMarkerStep(step);
            //     std::cout << "Integration is over, point in proc " << _mproc << std::endl;
          }
        }

        if(step != UINT_MAX && markerOutsideDomain == true) {  //prima era else if (step != UINT_MAX)
          integrationIsOverCounterProc[_iproc] += 1;
          step = UINT_MAX;
          _particles[iMarker]->SetIprocMarkerStep(step);
        }

      }
      //END LOCAL ADVECTION INSIDE IPROC



      // std::cout << "PRIMA integrationIsOverCounter = " << integrationIsOverCounter << std::endl;


      for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
        MPI_Bcast(& integrationIsOverCounterProc[jproc], 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
        integrationIsOverCounter += integrationIsOverCounterProc[jproc];
      }

      // std::cout << "DOPO integrationIsOverCounter = " << integrationIsOverCounter << std::endl;

      //BEGIN exchange on information

      for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
        for(unsigned iMarker = _markerOffset[jproc]; iMarker < _markerOffset[jproc + 1]; iMarker++) {
          unsigned elem =  _particles[iMarker]->GetMarkerElement();
          MPI_Bcast(& elem, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
          _particles[iMarker]->SetMarkerElement(elem);

          unsigned step =  _particles[iMarker]->GetIprocMarkerStep();
          MPI_Bcast(& step, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
          _particles[iMarker]->SetIprocMarkerStep(step);

          if(elem != UINT_MAX) {  // if it is outside jproc //TODO ACTUALLY IF WE ARE HERE IT COULD STILL BE IN JPROC but no outside the domain
            unsigned mproc = _particles[iMarker]->GetMarkerProc();
            //          _particles[iMarker]->SetMarkerProc(mproc); //TODO perche' facciamo sta cosa?
            if(mproc != jproc) {
              unsigned prevElem = _particles[iMarker]->GetIprocMarkerPreviousElement();
              _particles[iMarker]->GetElement(prevElem, jproc);
              _particles[iMarker]->SetIprocMarkerPreviousElement(prevElem);
            }
            elem = _particles[iMarker]->GetMarkerElement(); //TODO don't we have to resend elem with a broadcast? we need it for update line
            if(elem != UINT_MAX) {  // if it is not outside the domain
              unsigned mproc = _particles[iMarker]->GetMarkerProc();
              if(mproc != jproc) {
                if(jproc == _iproc) {
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
                  std::vector < std::vector < std::vector < double > > >().swap(aX);

                }
                else if(mproc == _iproc) {
                  _particles[iMarker]->InitializeX0andK(order);
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
//               if(mproc == _iproc) { //WARNING this now should be outside and was causing the problem with different processes
//                 _particles[iMarker]-> FindLocalCoordinates(solVType, aX, true);
//               }
            }
          }
        }
      }

      //END exchange of information


      //BEGIN to be remove
      std::vector<double> tr;
      unsigned step;
      for(unsigned j = 0; j < _size; j++) {
        _particles[_printList[j]]->GetMarkerCoordinates(tr);
        step = _particles[_printList[j]]->GetIprocMarkerStep();
        for(unsigned i = 0; i < _dim; i++) {
          std::cout << "x[ " << j << " ][ " << i << " ] = " << tr[i] ;
        }
        tr = _particles[_printList[j]]->GetIprocMarkerOldCoordinates();
        for(unsigned i = 0; i < _dim; i++) {
          std::cout << "x0[ " << j << " ][ " << i << " ] = " << tr[i] ;
        }


        std::cout << " step = " << step << "currentElem = " << _particles[_printList[j]]->GetMarkerElement() << " previousElem = " <<  _particles[_printList[j]]->GetIprocMarkerPreviousElement() << std::endl;
      }
      //END to be removed

      UpdateLine();

      for(unsigned j = 0; j < _size; j++) {
        _particles[_printList[j]]->GetMarkerCoordinates(tr);
        step = _particles[_printList[j]]->GetIprocMarkerStep();
        for(unsigned i = 0; i < _dim; i++) {
          std::cout << "x[ " << j << " ][ " << i << " ] = " << tr[i] ;
        }
        tr = _particles[_printList[j]]->GetIprocMarkerOldCoordinates();
        for(unsigned i = 0; i < _dim; i++) {
          std::cout << "x0[ " << j << " ][ " << i << " ] = " << tr[i] ;
        }
        std::cout << " step = " << step << "currentElem = " << _particles[j]->GetMarkerElement() << " previousElem = " <<  _particles[j]->GetIprocMarkerPreviousElement() << std::endl;
      }


      // std::cout << " END integrationIsOverCounter = " << integrationIsOverCounter <<  std::endl;

    }
  }
}

