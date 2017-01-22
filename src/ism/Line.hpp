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
        _markerOffset.resize(_nprocs);
        _particles.resize(size);
        _printList.resize(size);

        for(unsigned j = 0; j < size; j++) {
          particles[j] = new Marker(x[j], markerType[j], mesh, solType, debug);
        }


        //BEGIN reorder the markers by proc
        unsigned counter = 0;
        for(unsigned iproc = 0; iproc < _nprocs; iproc++) {
          for(unsigned j = 0; j < size; j++) {
            unsigned markerProc = particles[j]->GetMarkerProc();
            if(markerProc == iproc) {
              _particles[counter] = particles[j];
              _printList[j] = counter;
              counter++;
            }
          }
          _markerOffset[iproc] = counter;
        }
        //END reorder the markers by proc

        //BEGIN reorder markers also by element

        for(unsigned j = 0; j < _size; j++) {
          particles[j] = _particles[j];
        }

        unsigned* elementList = new unsigned [mesh->GetNumberOfElements()];
        for(unsigned iFlag = 0; iFlag < mesh->GetNumberOfElements(); iFlag++) {
          elementList[iFlag] = 0;
        }

        for(unsigned iproc = 0; iproc < _nprocs; iproc++) {

          unsigned counter = 1;

          for(unsigned jp = _markerOffset[iproc]; jp < _markerOffset[iproc + 1]; jp++) {

            unsigned jel = particles[jp]->GetMarkerElement();

            if(elementList[jel] != 1) {

              elementList[jel] = 1;

              for(unsigned ip = _markerOffset[iproc]; ip < _markerOffset[iproc + 1]; ip++) {
                unsigned iel = _particles[ip]->GetMarkerElement();
                if(ip != jp && iel == jel) {
                  elementList[iel] = 1;
                  _particles[_markerOffset[iproc] + counter] = particles[ip];
                  _printList[ip] = _markerOffset[iproc] + counter; //TODO FIX
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


        //BEGIN TEST TO SEE IF IT's TRUE,  IT WORKS :)

        for(unsigned i = 0; i < size; i++) {
          std::cout << "Particle: " << i << " , " << "Processor: " << _particles[i]->GetMarkerProc() << " , "
                    << "Element: " << _particles[i]->GetMarkerElement() << " " << std::endl;
        }

        //END TESTS

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




