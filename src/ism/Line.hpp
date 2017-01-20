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
              counter++;
            }
          }
          _markerOffset[iproc] = counter;
        }
        //END reorder the markers by proc

        _line.resize(1);
        _line[0].resize(size + 1);

        for(unsigned j = 0; j < size; j++) {
          _particles[j]->GetMarkerCoordinates(_line[0][j]);
        }
        _particles[0]->GetMarkerCoordinates(_line[0][size]);


        for(unsigned j = 0; j < size; j++) {
          delete particles[j];
        }

      };


     std::vector < std::vector < std::vector < double > > > GetLine() {
        return _line;
      }

      void UpdateParticles();
      
      void AdvectionParallel(Solution* sol, const unsigned &n, const double& T, const unsigned &order);
      
    private:

      std::vector < std::vector < std::vector < double > > > _line;
      std::vector < Marker*> _particles;
      std::vector < unsigned> _markerOffset;
      unsigned _size;
      unsigned _dim;
      
      std::vector < std::vector < std::vector < double > > > _K;
      
            static const double _a[4][4][4];
      static const double _b[4][4];
      static const double _c[4][4];

  };
} //end namespace femus



#endif



