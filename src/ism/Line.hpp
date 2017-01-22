/*=========================================================================

 Program: FEMuS
 Module: Line
 Authors: Eugenio Aulisa and Giacomo Capodaglio

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/


#ifndef __femus_ism_Line_hpp__
#define __femus_ism_Line_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MarkerTypeEnum.hpp"
#include "ParallelObject.hpp"
#include "Mesh.hpp"
#include "Marker.hpp"

#include "vector"
#include "map"
#include "MyVector.hpp"

namespace femus {

  class Line : public ParallelObject {
    public:
      Line(std::vector < std::vector < double > > x, const std::vector <MarkerType> &markerType,
           Mesh *mesh, const unsigned & solType, const unsigned &size, const bool &debug = false) {

        std::vector < Marker*> particles(size);

        _dim = mesh->GetDimension();
        _size = size;

        _markerOffset.resize(_nprocs + 1);
        _markerOffset[_nprocs] = size;

        _particles.resize(size);
        _printList.resize(size);

        for(unsigned j = 0; j < size; j++) {
          particles[j] = new Marker(x[j], markerType[j], mesh, solType, debug);
        }

//         std::cout << "FIRST OF ALL" << std::endl;
//         for(unsigned i = 0; i < size; i++) {
//           std::cout << "Particle: " << i << " , " << "Processor: " << particles[i]->GetMarkerProc() << " , "
//                     << "Element: " << particles[i]->GetMarkerElement() << " " << std::endl;
//         }


        //BEGIN reorder the markers by proc
        unsigned counter = 0;
        for(unsigned iproc = 0; iproc < _nprocs; iproc++) {
          bool markersInProc = false;
          unsigned  offsetCounter = counter;
          for(unsigned j = 0; j < size; j++) {
            unsigned markerProc = particles[j]->GetMarkerProc();
            if(markerProc == iproc) {
              _particles[counter] = particles[j];
              _printList[j] = counter;
              counter++;
              markersInProc = true;
            }
          }
          if(markersInProc == true) {
            if(iproc == 0) {
              _markerOffset[iproc] = 0;
            }
            else {
              unsigned firstCounter = 0;
              for(unsigned jproc = 0; jproc < iproc; jproc++) {
                if(_markerOffset[jproc] == UINT_MAX) firstCounter++;
              }
              if(firstCounter == iproc) _markerOffset[iproc] = 0;
              else {
                _markerOffset[iproc] = offsetCounter;
              }
            }
          }
          else {
            _markerOffset[iproc] = UINT_MAX;
          }
//           std::cout << "_markerOffset[" << iproc << "]= " << _markerOffset[iproc] << std::endl;
        }
        //END reorder the markers by proc


//         std::cout << "AFTER THE REORDERING BY PROC" << std::endl;
//         for(unsigned i = 0; i < size; i++) {
//           std::cout << "Particle: " << i << " , " << "Processor: " << _particles[i]->GetMarkerProc() << " , "
//                     << "Element: " << _particles[i]->GetMarkerElement() << " " << std::endl;
//         }


        //BEGIN reorder markers also by element

        for(unsigned j = 0; j < _size; j++) {
          particles[j] = _particles[j];
        }

        unsigned* elementList = new unsigned [mesh->GetNumberOfElements()];
        std::vector < unsigned> printList(size);

        //flags to see if we already ordered by that element, 0 = not considered yet
        for(unsigned iFlag = 0; iFlag < mesh->GetNumberOfElements(); iFlag++) {
          elementList[iFlag] = 0;
        }

        printList = _printList;

        for(unsigned iproc = 0; iproc < _nprocs; iproc++) {

          if(_markerOffset[iproc] != UINT_MAX) {

            counter = 0;

            unsigned upperBound;
            bool upperBoundFound = 0;
            bool finishedProcs = 0;
            unsigned lproc = iproc + 1;
            while(upperBoundFound + finishedProcs == 0) {
              if(_markerOffset[lproc] != UINT_MAX) {
                upperBound = _markerOffset[lproc];
                upperBoundFound = 1;
              }
              else {
                if(lproc < _nprocs) lproc ++;
                else finishedProcs = 1;
              }
            }

            
            std::cout << "upperBound =" << upperBound << std::endl;
                
            
            for(unsigned jp = _markerOffset[iproc]; jp < upperBound; jp++) {

              unsigned jel = particles[jp]->GetMarkerElement();

              if(elementList[jel] == 0) {

                elementList[jel] = 1;

                _particles[_markerOffset[iproc] + counter] = particles[jp];
                _printList[jp] = _markerOffset[iproc] + counter;
                counter++;


                for(unsigned ip = _markerOffset[iproc]; ip < upperBound; ip++) {
                  unsigned iel = particles[ip]->GetMarkerElement();
                  if(ip != jp && iel == jel) {
                    //std::cout << " jel =" << jel << " , " << "ip = " << ip << " , " << "jp = " << jp <<  std::endl;
                    elementList[iel] = 1;
                    _particles[_markerOffset[iproc] + counter] = particles[ip];
                    _printList[ip] = _markerOffset[iproc] + counter; //TODO FIX
                    counter++;
                  }
                }
              }
            }
          }
        }

        //END reorder markers also by element

//         for(unsigned i = 0; i < _size; i++) {
//           std::cout << "_printList[ " << i << "] = " << _printList[i] << std::endl;
//         }


//         std::cout << "AFTER THE REORDERING BY ELEM" << std::endl;
//         for(unsigned i = 0; i < size; i++) {
//           std::cout << "Particle: " << i << " , " << "Processor: " << _particles[i]->GetMarkerProc() << " , "
//                     << "Element: " << _particles[i]->GetMarkerElement() << " " << std::endl;
//         }


        _line.resize(1);
        _line[0].resize(size + 1);

        for(unsigned j = 0; j < size; j++) {
          _particles[_printList[j]]->GetMarkerCoordinates(_line[0][j]);
        }
        _particles[_printList[0]]->GetMarkerCoordinates(_line[0][size]);


        for(unsigned j = 0; j < size; j++) {
          delete particles[j];
        }
        delete[] elementList;
        std::vector < unsigned > ().swap(printList);

      };


      std::vector < std::vector < std::vector < double > > > GetLine() {
        return _line;
      }

      void UpdateLine(); //TODO update the function

      void AdvectionParallel(Solution* sol, const unsigned &n, const double& T, const unsigned &order);

    private:

      std::vector < std::vector < std::vector < double > > > _line;
      std::vector < Marker*> _particles;
      std::vector < unsigned > _markerOffset;
      std::vector < unsigned > _printList;
      unsigned _size;
      unsigned _dim;

      std::vector < std::vector < std::vector < double > > > _K;

      static const double _a[4][4][4];
      static const double _b[4][4];
      static const double _c[4][4];

  };
} //end namespace femus



#endif




