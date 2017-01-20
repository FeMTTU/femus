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
	
	unsigned* printList = new unsigned [_size];

        for(unsigned j = 0; j < size; j++) {
          particles[j] = new Marker(x[j], markerType[j], mesh, solType, debug);
	  printList[j] = j;
        }


        //BEGIN reorder the markers by proc
        unsigned counter = 0;
        for(unsigned iproc = 0; iproc < _nprocs; iproc++) {
          for(unsigned j = 0; j < size; j++) {
            unsigned markerProc = particles[j]->GetMarkerProc();
            if(markerProc == iproc) {

//               //BEGIN TESTS PART 1
// 	      std::cout << " PRIMA PARTE" << std::endl;
//               std::cout << "proc = " << particles[j]->GetMarkerProc() << std::endl;
//               std::cout << "elem = " << particles[j]->GetMarkerElement() << std::endl;
//               std::vector < double > trial(_dim, 0);
//               particles[j]->GetMarkerCoordinates(trial);
//               for(unsigned l = 0; l < _dim; l++) {
//                 std::cout << "trail[" << l << "]=" << trial[l] << std::endl;
//               }
//               trial = particles[j]->GetMarkerLocalCoordinates();
//               for(unsigned l = 0; l < _dim; l++) {
//                 std::cout << "xi trail[" << l << "]=" << trial[l] << std::endl;
//               }
//               //END TESTS PART 1

              _particles[counter] = particles[j];
	      _printList[j] = counter;

//               //BEGIN TESTS PART 2
//               std::cout << " SECONDA PARTE" << std::endl;
// 	      std::cout << "proc = " << _particles[counter]->GetMarkerProc() << std::endl;
//               std::cout << "elem = " << _particles[counter]->GetMarkerElement() << std::endl;
//               _particles[counter]->GetMarkerCoordinates(trial);
//               for(unsigned l = 0; l < _dim; l++) {
//                 std::cout << "trail[" << l << "]=" << trial[l] << std::endl;
//               }
//               trial = _particles[counter]->GetMarkerLocalCoordinates();
//               for(unsigned l = 0; l < _dim; l++) {
//                 std::cout << "xi trail[" << l << "]=" << trial[l] << std::endl;
//               }
//               //END TESTS PART 2


              counter++;
            }
          }
          _markerOffset[iproc] = counter;
        }
        //END reorder the markers by proc

        _line.resize(1);
        _line[0].resize(size + 1);

        for(unsigned j = 0; j < size; j++) {
          _particles[_printList[j]]->GetMarkerCoordinates(_line[0][j]);
        }
        _particles[_printList[0]]->GetMarkerCoordinates(_line[0][size]);


        for(unsigned j = 0; j < size; j++) {
          delete particles[j];
        }
        delete[] printList;

      };


      std::vector < std::vector < std::vector < double > > > GetLine() {
        return _line;
      }

      void UpdateParticles();

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



