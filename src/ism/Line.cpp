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
#include "Line.hpp"
#include "Marker.hpp"
#include "NumericVector.hpp"
#include "PolynomialBases.hpp"

#include <cmath>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using boost::math::ellint_1;
using boost::math::ellint_2;

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


  Line::Line (const std::vector < std::vector < double > > x, const std::vector < double > &mass,
              const std::vector <MarkerType>& markerType,
              const Solution* sol, const unsigned& solType)
     :
    _mesh ( sol->GetMesh() ),
    _sol ( sol)
              {

    _time.assign (10, 0);

    _size = x.size();

    std::vector < Marker*> particles (_size);

    _dim = _mesh->GetDimension();

    _markerOffset.resize (_nprocs + 1);
    _markerOffset[_nprocs] = _size;

    _particles.resize (_size);
    _printList.resize (_size);

    for (unsigned j = 0; j < _size; j++) {
      particles[j] = new Marker (x[j], mass[j], markerType[j], _sol, solType, true);
    }
    Reorder (particles);

  }

  Line::Line (const std::vector < std::vector < double > > x,
              const std::vector <MarkerType>& markerType,
              const Solution* sol, const unsigned& solType)
  :
    _mesh ( sol->GetMesh() ),
    _sol ( sol)
 {

    _time.assign (10, 0);

    _size = x.size();

    std::vector < Marker*> particles (_size);

    _dim = _mesh->GetDimension();

    _markerOffset.resize (_nprocs + 1);
    _markerOffset[_nprocs] = _size;

    _particles.resize (_size);
    _printList.resize (_size);

    for (unsigned j = 0; j < _size; j++) {
      particles[j] = new Marker (x[j], 0., markerType[j], _sol, solType, true);
    }
    Reorder (particles);
  }

  Line::Line (const std::vector < std::vector < double > > x,
              const std::vector < std::vector < std::vector < double > > > &tangent,
              const std::vector <MarkerType>& markerType,
              const Solution* sol, const unsigned& solType)
    :
    _mesh ( sol->GetMesh() ),
    _sol ( sol)
 {


    _time.assign (10, 0);

    _size = x.size();

    std::vector < Marker*> particles (_size);

    _dim = _mesh->GetDimension();

    _markerOffset.resize (_nprocs + 1);
    _markerOffset[_nprocs] = _size;

    _particles.resize (_size);
    _printList.resize (_size);

    for (unsigned j = 0; j < _size; j++) {
      particles[j] = new Marker (x[j], 0., markerType[j], _sol, solType, true);
      particles[j]->SetMarkerTangentGlobal (tangent[j]);
    }
    Reorder (particles);
  }

  void Line::Reorder (std::vector < Marker*> &particles) {

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
        unsigned markerProc = particles[j]->GetMarkerProc (_sol);

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


    //BEGIN (commented because there is another reorder function down below) reorder markers also by element

//     for (unsigned j = 0; j < _size; j++) {
//       particles[j] = _particles[j];
//     }
//
//     std::vector < unsigned> printList (_size);
//     printList = _printList;
//
//     unsigned meshElements = _mesh->GetNumberOfElements();
//     std::vector< bool > elementList (meshElements, false);
//
//
//     //flags to see if we already ordered by that element, 0 = not considered yet
// //     for(unsigned iFlag = 0; iFlag < meshElements; iFlag++) {
// //       elementList[iFlag] = 0;
// //     }
//
//     for (unsigned iproc = 0; iproc < _nprocs; iproc++) {
//
//       bool someMarkersOutsideDomain = false;
//       counter = 0;
//
//       for (unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {
//
//         unsigned jel;
//         jel = particles[jp]->GetMarkerElement();
//
//         if (jel != UINT_MAX) {
//
//           if (elementList[jel] == false) {
//
//             elementList[jel] = true;
//
//             _particles[_markerOffset[iproc] + counter] = particles[jp];
//
//             for (unsigned iList = 0; iList < _size; iList++) {
//               if (printList[iList] == jp) {
//                 _printList[iList] = _markerOffset[iproc] + counter;
//                 break;
//               }
//             }
//
//             counter++;
//
//
//             for (unsigned ip = jp + 1; ip < _markerOffset[iproc + 1]; ip++) {
//               unsigned iel;
//               iel = particles[ip]->GetMarkerElement();
//
//               if (iel == jel) {
//                 _particles[_markerOffset[iproc] + counter] = particles[ip];
//
//                 for (unsigned iList = 0; iList < _size; iList++) {
//                   if (printList[iList] == ip) {
//                     _printList[iList] = _markerOffset[iproc] + counter;
//                     break;
//                   }
//                 }
//
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
//       if (someMarkersOutsideDomain == true) {
//         if (iproc == 0) {
//           for (unsigned i = 0; i < _size; i++) {
//             unsigned iel = particles[i]->GetMarkerElement();
//
//             if (iel == UINT_MAX) {
//
//               _particles[_markerOffset[0] + counter] = particles[i];
//
//               for (unsigned iList = 0; iList < _size; iList++) {
//                 if (printList[iList] == i) {
//                   _printList[iList] = _markerOffset[0] + counter;
//                   break;
//                 }
//               }
//
//               counter++;
//             }
//           }
//         }
//         else {
//           std::cout << "warning" << std::endl;
//           std::cout << "the marker is outside the domain but not in iproc 0" << std::endl;
//           abort();
//         }
//       }
//
//     }

    //END (commented because there is another reorder function down below) reorder markers also by element



//     //BEGIN OLD (don't use) reorder markers also by element
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
//  std::cout<< " Marker in element jel = " << jel << std::endl;
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
//     //END OLD (don't use) reorder markers also by element

//         for(unsigned i = 0; i < _size; i++) {
//           std::cout << "_printList[ " << i << "] = " << _printList[i] << std::endl;
//         }


//     std::cout << "AFTER THE REORDERING BY ELEM" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc =_particles[i]->GetMarkerProc (_sol);
//       unsigned elem;
//       _particles[i]->GetMarkerElementLine(elem);
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " << "Pointer: " << _particles[i] << std::endl;
//     }

    //BEGIN reorder markers by element but in order  (i.e. first the first element of the proc and so on)

    for (unsigned j = 0; j < _size; j++) {
      particles[j] = _particles[j];
    }

    std::vector < unsigned> printList2 (_size);
    printList2 = _printList;

    std::vector < bool > particleHasBeenChecked (_size, false);

    for (unsigned iproc = 0; iproc < _nprocs; iproc++) {

      bool someMarkersOutsideDomain = false;
      counter = 0;
      unsigned iprocParticles = _markerOffset[iproc + 1] - _markerOffset[iproc];

      for (unsigned iel = _mesh->GetElementOffset(iproc); iel < _mesh->GetElementOffset(iproc + 1); iel++) {

        for (unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {

          if (!particleHasBeenChecked[jp]) {

            unsigned jel;
            jel = particles[jp]->GetMarkerElement();

            if (jel == UINT_MAX) {
              someMarkersOutsideDomain = true;
              particleHasBeenChecked[jp] = true;
            }

            else {
              if (iel == jel) {

                particleHasBeenChecked[jp] = true;
                _particles[_markerOffset[iproc] + counter] = particles[jp];

                for (unsigned iList = 0; iList < _size; iList++) {
                  if (printList2[iList] == jp) {
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

      if (someMarkersOutsideDomain == true) {
        if (iproc == 0) {
          for (unsigned i = 0; i < _size; i++) {
            unsigned iel = particles[i]->GetMarkerElement();
            if (iel == UINT_MAX) {

              _particles[_markerOffset[0] + counter] = particles[i];

              for (unsigned iList = 0; iList < _size; iList++) {
                if (printList2[iList] == i) {
                  _printList[iList] = _markerOffset[0] + counter;
                  break;
                }
              }
              counter++;
            }
          }
        }
        else {
          std::cout << "warning" << std::endl;
          std::cout << "the marker is outside the domain but not in iproc 0" << std::endl;
          abort();
        }
      }

      //check that all particles have been considered
      unsigned check2 = 0;
      for (unsigned icheck = _markerOffset[iproc]; icheck < _markerOffset[iproc + 1]; icheck++) {
        if (particleHasBeenChecked[icheck]) check2++;
      }
//        std::cout<<"check2 = " << check2 <<std::endl;
//        std::cout<<"iprocParticles = " << iprocParticles <<std::endl;
      if (check2 != iprocParticles) std::cout << "------------------------- WARNING NOT ALL PARTICLES HAVE BEEN CHECKED IN PROC " << iproc << std::endl;
    }

//     std::cout << "AFTER THE REORDERING BY ELEM IN ORDER (this is line 523 in Line.cpp)" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//         unsigned proc =_particles[i]->GetMarkerProc (_sol);
//         unsigned elem;
//         _particles[i]->GetMarkerElementLine(elem);
//         std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , " << "Element: " << elem << " " << "Pointer: " << _particles[i] << std::endl;
//     }


    //END reorder markers by element but in order  (i.e. first the first element of the proc and so on)

    _line.resize (_size + 1);

    for (unsigned j = 0; j < _size; j++) {
      _particles[_printList[j]]->GetMarkerCoordinates (_line[j]);
    }

    _particles[_printList[0]]->GetMarkerCoordinates (_line[_size]);

    // delete[] elementList;
    // std::vector < unsigned > ().swap(printList);


//     for(unsigned i=0;i<_size;i++){
//       if(_printList[i] == 3395){
//  std::cout << i<< std::endl;
//       }
//     }
//
//     exit(0);
  };

  Line::~Line() {
    for (unsigned j = 0; j < _size; j++) {
      //std::cout << j << " " << _particles[j] << std::endl<< std::flush;
      delete _particles[j];
    }
  }

  void Line::UpdateLine() {

    std::vector < Marker*> particles (_size);
    std::vector < unsigned> printList (_size);

    for (unsigned j = 0; j < _size; j++) {
      particles[j] = _particles[j];
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
        unsigned markerProc = particles[j]->GetMarkerProc (_sol);

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



    //BEGIN (commented because there is another reorder function down below) reorder markers also by element

//     for (unsigned j = 0; j < _size; j++) {
//       particles[j] = _particles[j];
//     }
//
//     printList = _printList;
//
//     unsigned meshElements = _mesh->GetNumberOfElements();
//     std::vector< bool > elementList (meshElements, false);
//
//
//     //flags to see if we already ordered by that element, 0 = not considered yet
// //     for(unsigned iFlag = 0; iFlag < meshElements; iFlag++) {
// //       elementList[iFlag] = 0;
// //     }
//
//     for (unsigned iproc = 0; iproc < _nprocs; iproc++) {
//
//       bool someMarkersOutsideDomain = false;
//       counter = 0;
//
//       for (unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {
//
//         unsigned jel;
//         jel = particles[jp]->GetMarkerElement();
//
//         if (jel != UINT_MAX) {
//
//           if (elementList[jel] == false) {
//
//             elementList[jel] = true;
//
//             _particles[_markerOffset[iproc] + counter] = particles[jp];
//
//             for (unsigned iList = 0; iList < _size; iList++) {
//               if (printList[iList] == jp) {
//                 _printList[iList] = _markerOffset[iproc] + counter;
//                 break;
//               }
//             }
//
//             counter++;
//
//
//             for (unsigned ip = jp + 1; ip < _markerOffset[iproc + 1]; ip++) {
//               unsigned iel;
//               iel = particles[ip]->GetMarkerElement();
//
//               if (iel == jel) {
//                 _particles[_markerOffset[iproc] + counter] = particles[ip];
//
//                 for (unsigned iList = 0; iList < _size; iList++) {
//                   if (printList[iList] == ip) {
//                     _printList[iList] = _markerOffset[iproc] + counter;
//                     break;
//                   }
//                 }
//
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
//       if (someMarkersOutsideDomain == true) {
//         if (iproc == 0) {
//           for (unsigned i = 0; i < _size; i++) {
//             unsigned iel = particles[i]->GetMarkerElement();
//
//             if (iel == UINT_MAX) {
//
//               _particles[_markerOffset[0] + counter] = particles[i];
//
//               for (unsigned iList = 0; iList < _size; iList++) {
//                 if (printList[iList] == i) {
//                   _printList[iList] = _markerOffset[0] + counter;
//                   break;
//                 }
//               }
//
//               counter++;
//             }
//           }
//         }
//         else {
//           std::cout << "warning" << std::endl;
//           std::cout << "the marker is outside the domain but not in iproc 0" << std::endl;
//           abort();
//         }
//       }
//
//     }

    //END (commented because there is another reorder function down below) reorder markers also by element


    //BEGIN reorder markers by element but in order  (i.e. first the first element of the proc and so on)

    for (unsigned j = 0; j < _size; j++) {
      particles[j] = _particles[j];
    }

    std::vector < unsigned> printList2 (_size);
    printList2 = _printList;

    std::vector < bool > particleHasBeenChecked (_size, false);

    for (unsigned iproc = 0; iproc < _nprocs; iproc++) {

      bool someMarkersOutsideDomain = false;
      counter = 0;
      unsigned iprocParticles = _markerOffset[iproc + 1] - _markerOffset[iproc];

      for (unsigned iel = _mesh->GetElementOffset(iproc); iel < _mesh->GetElementOffset(iproc + 1); iel++) {

        for (unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {

          if (!particleHasBeenChecked[jp]) {

            unsigned jel;
            jel = particles[jp]->GetMarkerElement();

            if (jel == UINT_MAX) {
              someMarkersOutsideDomain = true;
              particleHasBeenChecked[jp] = true;
            }

            else {
              if (iel == jel) {

                particleHasBeenChecked[jp] = true;
                _particles[_markerOffset[iproc] + counter] = particles[jp];

                for (unsigned iList = 0; iList < _size; iList++) {
                  if (printList2[iList] == jp) {
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

      if (someMarkersOutsideDomain == true) {
        if (iproc == 0) {
          for (unsigned i = 0; i < _size; i++) {
            unsigned iel = particles[i]->GetMarkerElement();
            if (iel == UINT_MAX) {

              _particles[_markerOffset[0] + counter] = particles[i];

              for (unsigned iList = 0; iList < _size; iList++) {
                if (printList2[iList] == i) {
                  _printList[iList] = _markerOffset[0] + counter;
                  break;
                }
              }
              counter++;
            }
          }
        }
        else {
          std::cout << "warning" << std::endl;
          std::cout << "the marker is outside the domain but not in iproc 0" << std::endl;
          abort();
        }
      }

      //check that all particles have been considered
      unsigned check2 = 0;
      for (unsigned icheck = _markerOffset[iproc]; icheck < _markerOffset[iproc + 1]; icheck++) {
        if (particleHasBeenChecked[icheck]) check2++;
      }
//        std::cout<<"check2 = " << check2 <<std::endl;
//        std::cout<<"iprocParticles = " << iprocParticles <<std::endl;
      if (check2 != iprocParticles) std::cout << "------------------------- WARNING NOT ALL PARTICLES HAVE BEEN CHECKED IN PROC " << iproc << std::endl;
    }

    //END reorder markers by element but in order  (i.e. first the first element of the proc and so on)

    for (unsigned j = 0; j < _size; j++) {

      _particles[_printList[j]]->GetMarkerCoordinates (_line[j]);
    }

    _particles[_printList[0]]->GetMarkerCoordinates (_line[_size]);

    //delete[] elementList;
    //std::vector < unsigned > ().swap(printList);
    //std::vector < unsigned > ().swap(elementList);


//     std::cout << "-----------------------------------------  UPDATE LINE  BY ELEMENTS-----------------------------------------" << std::endl;
//     for(unsigned i = 0; i < _size; i++) {
//       unsigned proc;
//       unsigned elem;
//       //std::vector <double> chiappe;
//       //chiappe = _particles[_printList[i]]->GetIprocMarkerCoordinates();
//       proc = _particles[i]->GetMarkerProc(_sol);
//       elem = _particles[i]->GetMarkerElement();
//       std::cout << "Particle: " << i << " , " << "Processor: " << proc << " , "
//                 << "Element: " << elem << " " <<  " "  ;
// //       for(unsigned j = 0; j < _dim; j++) {
// //         std::cout << "x[" << j << "]=" << chiappe[j] << " " ;
// //       }
//       std::cout << std::endl;
//     }


  }

  void Line::AdvectionParallel (const unsigned& n, const double& T, const unsigned& order, ForceFunction force) {

    //BEGIN  Initialize the parameters for all processors

    double s;
    double PI = acos (-1.);
    std::vector <double> Fm (3, 0);  // magnetic force initialization
    std::vector < unsigned > solVIndex (_dim);
    solVIndex[0] = _sol->GetIndex ("U");   // get the position of "U" in the ml_sol object
    solVIndex[1] = _sol->GetIndex ("V");   // get the position of "V" in the ml_sol object

    if (_dim == 3) solVIndex[2] = _sol->GetIndex ("W");    // get the position of "V" in the ml_sol object

    unsigned solVType = _sol->GetSolutionType (solVIndex[0]);   // get the finite element type for "u"

    std::vector < double > phi;
    std::vector < std::vector<double > > V (2);
    std::map<unsigned, std::vector < std::vector < std::vector < double > > > > aV;
    std::map<unsigned, std::vector < std::vector < std::vector < std::vector < double > > > > > aX;
    double h = T / n;

    //END


    //BEGIN declare marker instances
    unsigned step;
    std::vector < double > x (_dim);
    std::vector < double > x0 (_dim);
    std::vector < std::vector < double > > K (order);

    for (unsigned j = 0; j < order; j++) {
      K[j].resize (_dim);
    }

    //END

    unsigned integrationIsOverCounter = 0; // when integrationIsOverCounter = _size (which is the number of particles) it means all particles have been advected;
    // when the step of each marker is = n * order and/or marker element is UINT_MAX we increase by one such counter

    //BEGIN Numerical integration scheme

    for (unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {
      _particles[iMarker]->InitializeMarkerForAdvection (order);
    }

//     unsigned maxload = 0;
//     for(unsigned jproc=0; jproc<_nprocs;jproc++){
//       unsigned temp = (_markerOffset[jproc + 1] - _markerOffset[jproc]);
//       maxload = ( temp > maxload) ? temp : maxload;
//     }
//     maxload *= (n * order)/(_nprocs * _nprocs);

    while (integrationIsOverCounter != _size) {

      MyVector <unsigned> integrationIsOverCounterProc (1, 0);
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

            bool elementUpdate = (aX.find (currentElem) != aX.end()) ? false : true;    //update if currentElem was never updated

            clock_t localTime = clock();
            _particles[iMarker]->GetMarkerS (n, order, s);
            _particles[iMarker]->FindLocalCoordinates (solVType, aX[currentElem], elementUpdate, _sol, s);
            _particles[iMarker]->updateVelocity (V, solVIndex, solVType, aV[currentElem], phi, elementUpdate, _sol);  // we put pcElemUpdate instead of true but it wasn't running
            _time[3] += static_cast<double> ( (clock() - localTime)) / CLOCKS_PER_SEC;

            unsigned tstep = step / order;
            unsigned istep = step % order;

            if (istep == 0) {
              x0 = x;

              for (unsigned j = 0; j < order; j++) {
                K[j].assign (_dim, 0.);
              }
            }

            // double s = (tstep + _c[order - 1][istep]) / n;


            //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

            if (force != NULL) {
              //else if (_sol->GetIfFSI()) {
              unsigned material = _sol->GetMesh()->GetElementMaterial (currentElem);
//               MagneticForce(x, Fm, material, 0);
              force (x, Fm, material);
            }

            //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::flush;

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
            }
            else {
              for (unsigned i = 0; i < _dim; i++) {
                x[i] = x0[i];

                for (unsigned j = 0; j < order; j++) {
                  x[i] += _b[order - 1][j] * K[j][i];
                }
              }
            }

            _particles[iMarker]->SetIprocMarkerOldCoordinates (x0);
            _particles[iMarker]->SetIprocMarkerCoordinates (x);
            _particles[iMarker]->SetIprocMarkerK (K);

            _particles[iMarker]->SetIprocMarkerStep (step);
            _particles[iMarker]->GetMarkerS (n, order, s);

            //std::cout << step <<" "<< istep<< " AAAAAAAAAAAAA"<<std::endl<<std::flush;
            unsigned previousElem = currentElem;
            localTime = clock();
            _particles[iMarker]->GetElementSerial (previousElem, _sol, s);
            _time[4] += static_cast<double> ( (clock() - localTime)) / CLOCKS_PER_SEC;
            localTime = clock();
            // std::cout << step <<" "<< istep<< " BBBBBBBBBBBBB"<<std::endl<<std::flush;

            _particles[iMarker]->SetIprocMarkerPreviousElement (previousElem);

            currentElem = _particles[iMarker]->GetMarkerElement();
            unsigned mproc = _particles[iMarker]->GetMarkerProc (_sol);

            if (currentElem == UINT_MAX) {   // the marker has been advected outside the domain
              markerOutsideDomain = true;
              step = UINT_MAX;
              _particles[iMarker]->SetIprocMarkerStep (step);
              break;
            }
            else if (_iproc != mproc) {   // the marker has been advected outise the process
              break;
            }
          }

          if (step == n * order) {
            step = UINT_MAX;
            _particles[iMarker]->SetIprocMarkerStep (step);
          }
        }
        else {   // the marker started outise the domain
          step = UINT_MAX;
          _particles[iMarker]->SetIprocMarkerStep (step);
        }

        if (step == UINT_MAX || markerOutsideDomain) {
          integrationIsOverCounterProc[_iproc] += 1;
        }

        //if(counter > maxload) break;
      }

      _time[5] += static_cast<double> ( (clock() - startTime)) / CLOCKS_PER_SEC;
      MPI_Barrier (PETSC_COMM_WORLD);
      _time[0] += static_cast<double> ( (clock() - startTime)) / CLOCKS_PER_SEC;
      startTime = clock();
      //END LOCAL ADVECTION INSIDE IPROC

      integrationIsOverCounter = 0;

      for (unsigned jproc = 0; jproc < _nprocs; jproc++) {
        integrationIsOverCounterProc.broadcast (jproc);
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
          MPI_Bcast (& elem, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
          _particles[iMarker]->SetMarkerElement (elem);

          unsigned step =  _particles[iMarker]->GetIprocMarkerStep();
          MPI_Bcast (& step, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
          _particles[iMarker]->SetIprocMarkerStep (step);

          if (elem != UINT_MAX) {   // if it is outside jproc, ACTUALLY IF WE ARE HERE IT COULD STILL BE IN JPROC but not outside the domain
            unsigned mproc = _particles[iMarker]->GetMarkerProc (_sol);
            _particles[iMarker]->SetMarkerProc (mproc);

            if (mproc != jproc) {
              unsigned prevElem = _particles[iMarker]->GetIprocMarkerPreviousElement();
              _particles[iMarker]->GetMarkerS (n, order, s);
              _particles[iMarker]->GetElement (prevElem, jproc, _sol, s);
              _particles[iMarker]->SetIprocMarkerPreviousElement (prevElem);
            }

            elem = _particles[iMarker]->GetMarkerElement();

            if (elem != UINT_MAX) {   // if it is not outside the domain
              unsigned mproc = _particles[iMarker]->GetMarkerProc (_sol);

              if (mproc != jproc) {
                if (jproc == _iproc) {

                  unsigned step =  _particles[iMarker]->GetIprocMarkerStep();
                  MPI_Send (& step, 1, MPI_UNSIGNED, mproc, order + 1, PETSC_COMM_WORLD);

                  unsigned istep = step % order;

                  if (istep != 0) {
                    K = _particles[iMarker]->GetIprocMarkerK();

                    for (int i = 0; i < order; i++) {
                      MPI_Send (&K[i][0], _dim, MPI_DOUBLE, mproc, i, PETSC_COMM_WORLD);
                    }

                    x0 = _particles[iMarker]->GetIprocMarkerOldCoordinates();
                    MPI_Send (&x0[0], _dim, MPI_DOUBLE, mproc, order, PETSC_COMM_WORLD);
                  }

                  _particles[iMarker]->FreeVariables();

                }
                else if (mproc == _iproc) {

                  _particles[iMarker]->InitializeVariables (order);

                  unsigned step;
                  MPI_Recv (& step, 1, MPI_UNSIGNED, jproc, order + 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  _particles[iMarker]->SetIprocMarkerStep (step);

                  unsigned istep = step % order;

                  if (istep != 0) {
                    for (int i = 0; i < order; i++) {
                      MPI_Recv (&K[i][0], _dim, MPI_DOUBLE, jproc, i, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                      _particles[iMarker]->SetIprocMarkerK (K);
                    }

                    MPI_Recv (&x0[0], _dim, MPI_DOUBLE, jproc, order, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                    _particles[iMarker]->SetIprocMarkerOldCoordinates (x0);

                  }
                }
              }
            }
          }

          if (elem == UINT_MAX && jproc != 0) {   // elem = UINT_MAX, but not yet in jproc = 0
            if (jproc == _iproc) {
              x = _particles[iMarker]->GetIprocMarkerCoordinates();
              MPI_Send (&x[0], _dim, MPI_DOUBLE, 0, 1, PETSC_COMM_WORLD);
              _particles[iMarker]->FreeVariables();
            }
            else if (_iproc == 0) {
              _particles[iMarker]->InitializeX();
              MPI_Recv (&x[0], _dim, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
              _particles[iMarker]->SetIprocMarkerCoordinates (x);
            }
          }
        }
      }

      MPI_Barrier (PETSC_COMM_WORLD);
      _time[1] += static_cast<double> ( (clock() - startTime)) / CLOCKS_PER_SEC;
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

      MPI_Barrier (PETSC_COMM_WORLD);
      _time[2] += static_cast<double> ( (clock() - startTime)) / CLOCKS_PER_SEC;
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


  unsigned Line::NumberOfParticlesOutsideTheDomain() {

    unsigned counter = 0;

    for (unsigned iMarker = _markerOffset[0]; iMarker < _markerOffset[1]; iMarker++) {
      unsigned elem =  _particles[iMarker]->GetMarkerElement();

      if (elem == UINT_MAX) {
        counter++;
      }
    }

    return counter;
  }



  void Line::GetParticlesToGridMaterial (const bool updateMat) {

    unsigned solIndexM, solTypeM, solIndexMat, solTypeMat;
    if (updateMat) {
      solIndexM = _sol->GetIndex ("M"); // nodes
      solTypeM = _sol->GetSolutionType (solIndexM);

      solIndexMat = _sol->GetIndex ("Mat"); // element
      solTypeMat = _sol->GetSolutionType (solIndexMat);

      _sol->_Sol[solIndexM]->zero();
      _sol->_Sol[solIndexMat]->zero();
    }

    std::map<unsigned, std::vector < std::vector < std::vector < std::vector < double > > > > > aX;
    // set all element with at least one marker to 3 and all nodes of the element to 1
    for (unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {

      unsigned iel = _particles[iMarker]->GetMarkerElement();
      unsigned ielType =  _mesh->GetElementType (iel);
      bool elementUpdate = (aX.find (iel) != aX.end()) ? false : true;    //update if iel was never updated

      std::vector <double> xi1 = _particles[iMarker]->GetMarkerLocalCoordinates();

      _particles[iMarker]->FindLocalCoordinates (2., aX[iel], elementUpdate, _sol, 0.);

      std::vector <double> xi = _particles[iMarker]->GetMarkerLocalCoordinates();

      if (updateMat) {
        for (unsigned j = 0; j < _mesh->GetElementDofNumber (iel, solTypeM); j++) {
          unsigned jdof = _mesh->GetSolutionDof (j, iel, solTypeM);
          _sol->_Sol[solIndexM]->set (jdof, 1.);
        }

        unsigned idofMat = _mesh->GetSolutionDof (0, iel, solTypeMat);
        _sol->_Sol[solIndexMat]->set (idofMat, 3.);
      }
    }
    if (updateMat) {
      _sol->_Sol[solIndexM]->close();
      _sol->_Sol[solIndexMat]->close();

      //set all nodes M with element Mat = 0 to 0
      for (int iel = _mesh->GetElementOffset(_iproc); iel < _mesh->GetElementOffset(_iproc + 1); iel++) {

        unsigned idofMat = _mesh->GetSolutionDof (0, iel, solTypeMat);
        unsigned  material = (*_sol->_Sol[solIndexMat]) (idofMat);

        unsigned nDofsM = _mesh->GetElementDofNumber (iel, solTypeM);  // number of mass dofs
        for (unsigned i = 0; i < nDofsM; i++) {
          unsigned idof = _mesh->GetSolutionDof (i, iel, solTypeM); // global to global mapping for mass solution
          if (material == 0 || (*_sol->_Bdc[solIndexM]) (idof) == 0) {
            _sol->_Sol[solIndexM]->set (idof, 0.);
          }
        }
      }
      _sol->_Sol[solIndexM]->close();

      for (int iel = _mesh->GetElementOffset(_iproc); iel < _mesh->GetElementOffset(_iproc + 1); iel++) {

        unsigned idofMat = _mesh->GetSolutionDof (0, iel, solTypeMat);
        unsigned  material = (*_sol->_Sol[solIndexMat]) (idofMat);
        if (material == 3) {
          unsigned nDofsM = _mesh->GetElementDofNumber (iel, solTypeM);  // number of mass dofs
          double counter = 0.;
          for (unsigned i = 0; i < nDofsM; i++) {
            unsigned idof = _mesh->GetSolutionDof (i, iel, solTypeM); // global to global mapping for mass solution
            double value = (*_sol->_Sol[solIndexM]) (idof);
            counter += value;
          }
          _sol->_Sol[solIndexMat]->set (idofMat, counter);
        }
      }
      _sol->_Sol[solIndexMat]->close();


      for (unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {

        unsigned iel = _particles[iMarker]->GetMarkerElement();
        unsigned ielType =  _mesh->GetElementType (iel);
        for (unsigned j = 0; j < _mesh->GetElementDofNumber (iel, solTypeM); j++) {

          unsigned jdof = _mesh->GetSolutionDof (j, iel, solTypeM);
          _sol->_Sol[solIndexM]->set (jdof, 1.);
        }
      }
      _sol->_Sol[solIndexM]->close();

    }
  }


  void Line::UpdateLineMPM() {


    for (unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {
      unsigned elem =  _particles[iMarker]->GetMarkerElement();
      _particles[iMarker]->GetElementSerial (elem, _sol, 0.);
      _particles[iMarker]->SetIprocMarkerPreviousElement (elem);
    }

    for (unsigned jproc = 0; jproc < _nprocs; jproc++) {
      for (unsigned iMarker = _markerOffset[jproc]; iMarker < _markerOffset[jproc + 1]; iMarker++) {

        unsigned elem =  _particles[iMarker]->GetMarkerElement();
        MPI_Bcast (& elem, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
        _particles[iMarker]->SetMarkerElement (elem);

        if (elem != UINT_MAX) {
          unsigned mproc = _particles[iMarker]->GetMarkerProc (_sol); //WARNING you don't know if this is your real process
          _particles[iMarker]->SetMarkerProc (mproc);

          if (mproc != jproc) { //this means, if we think the particle moved from jproc (which means element serial shouldn't have found the actual element)
            unsigned prevElem = _particles[iMarker]->GetIprocMarkerPreviousElement();
            _particles[iMarker]->GetElement (prevElem, jproc, _sol, 0.);
            _particles[iMarker]->SetIprocMarkerPreviousElement (prevElem);
          }

          elem = _particles[iMarker]->GetMarkerElement(); //TODO shouldn't this be inside the if ? If we don't go in the if I think we already know who elem is

          if (elem != UINT_MAX) {   // if it is not outside the domain
            unsigned mproc = _particles[iMarker]->GetMarkerProc (_sol); //actual mproc

            if (mproc != jproc) { //there is no need to send/receive if the particle didn't change process (which is when mproc == jproc)
              if (jproc == _iproc) {

                unsigned order = 0;
                std::vector <double> MPMQuantities = _particles[iMarker]->GetMPMQuantities();
                unsigned MPMsize =  _particles[iMarker]->GetMPMSize();
                MPI_Send (&MPMQuantities[0], MPMsize, MPI_DOUBLE, mproc, order, PETSC_COMM_WORLD);

                std::vector < std::vector < double > > Fp = _particles[iMarker]->GetDeformationGradient();
                for (unsigned i = 0; i < _dim; i++) {
                  MPI_Send (&Fp[i][0], _dim, MPI_DOUBLE, mproc, order + 1, PETSC_COMM_WORLD);
                }

                _particles[iMarker]->FreeVariables();

              }
              else if (mproc == _iproc) {

                unsigned order = 0;
                _particles[iMarker]->InitializeVariables (order);

                unsigned MPMsize =  _particles[iMarker]->GetMPMSize();
                std::vector <double> MPMQuantities (MPMsize);

                MPI_Recv (&MPMQuantities[0], MPMsize, MPI_DOUBLE, jproc, order, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                _particles[iMarker]->SetMPMQuantities (MPMQuantities);

                std::vector < std::vector < double > > Fp (_dim);
                for (unsigned i = 0; i < _dim; i++) {
                  Fp[i].resize (_dim);
                }

                for (unsigned i = 0; i < _dim; i++) {
                  MPI_Recv (&Fp[i][0], _dim, MPI_DOUBLE, jproc, order + 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                _particles[iMarker]->SetDeformationGradient (Fp);

              }
            }
          }
        }

        if (elem == UINT_MAX && jproc != 0) {   // elem = UINT_MAX, but not yet in jproc = 0
          if (jproc == _iproc) {
            std::vector<double> x (_dim);
            x = _particles[iMarker]->GetIprocMarkerCoordinates();
            MPI_Send (&x[0], _dim, MPI_DOUBLE, 0, 1, PETSC_COMM_WORLD);
            _particles[iMarker]->FreeVariables();
          }
          else if (_iproc == 0) {
            std::vector<double> x (_dim);
            _particles[iMarker]->InitializeX();
            MPI_Recv (&x[0], _dim, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
            _particles[iMarker]->SetIprocMarkerCoordinates (x);
          }
        }
      }
    }

    //END find new _elem and _mproc

    UpdateLine();

  }

  void Line::SetParticlesMass (const double& volume, const double& density) {
    double particlesMass = density * volume / _size;
    for (unsigned i = _markerOffset[_iproc]; i < _markerOffset[_iproc  + 1]; i++) {
      _particles[i]->SetMarkerMass (particlesMass);
    }
  }


  void Line::ScaleParticleMass (double scale (const std::vector <double>& x)) {
    for (unsigned i = _markerOffset[_iproc]; i < _markerOffset[_iproc  + 1]; i++) {
      std::vector<double> x (_dim);
      x = _particles[i]->GetIprocMarkerCoordinates();
      double mass = _particles[i]->GetMarkerMass();
      _particles[i]->SetMarkerMass (mass * scale (x));
    }

  }

  void Line::GetExtrema (std::vector <double>& xMin, std::vector <double>& xMax) {
    xMin.resize (_dim);
    xMax.resize (_dim);

    std::vector < double > xMinLocal (_dim, 1.0e100);
    std::vector < double > xMaxLocal (_dim, -1.0e100);
    std::vector< double > x;
    for (unsigned i = _markerOffset[_iproc]; i < _markerOffset[_iproc  + 1]; i++) {
      x = _particles[i]->GetIprocMarkerCoordinates();
      for (unsigned k = 0; k < _dim; k++) {
        xMinLocal[k] = (x[k] < xMinLocal[k]) ? x[k] : xMinLocal[k];
        xMaxLocal[k] = (x[k] > xMaxLocal[k]) ? x[k] : xMaxLocal[k];
      }
    }
    for (unsigned k = 0; k < _dim; k++) {
      MPI_Allreduce (&xMinLocal[0], &xMin[0], static_cast<int> (_dim), MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
      MPI_Allreduce (&xMaxLocal[0], &xMax[0], static_cast<int> (_dim), MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    }
  }


}







