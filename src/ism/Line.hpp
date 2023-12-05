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

#include <vector>
#include <map>
#include "MyVector.hpp"

namespace femus {

  class Line : public ParallelObject {
    public:

      Line (const std::vector < std::vector < double > > x,
            const std::vector <MarkerType>& markerType,
            Solution* sol, const unsigned& solType);

      Line (const std::vector < std::vector < double > > x, const std::vector < double > &mass,
            const std::vector <MarkerType>& markerType,
            Solution* sol, const unsigned& solType);

      Line (const std::vector < std::vector < double > > x,
            const std::vector < std::vector < std::vector < double > > > &tangent,
            const std::vector <MarkerType>& markerType,
            Solution* sol, const unsigned& solType);


      ~Line();

      typedef void (*ForceFunction) (const std::vector <double>& xMarker, std::vector <double>& Fm, const unsigned& material);

      void GetLine (std::vector < std::vector < double > >& line) {
        line = _line;
      }


      void GetStreamLine (std::vector < std::vector < std::vector < double > > >& line, const unsigned& step) {
        for (unsigned i = 0; i < _size; i++) {
          line[i].resize (step + 1);
          line[i][step] = _line[i];
        }
      }

      std::vector <unsigned> GetMarkerOffset() {
        return _markerOffset;
      }

      std::vector <Marker*> GetParticles() {
        return _particles;
      }

      void AdvectionParallel (const unsigned& n, const double& T, const unsigned& order, ForceFunction Force = NULL);

      void UpdateLine();

      unsigned NumberOfParticlesOutsideTheDomain();

      void GetParticlesToGridMaterial (const bool updateMat = true);

      void UpdateLineMPM();

      void SetParticlesMass (const double& volume, const double& density);

      void ScaleParticleMass (double scale (const std::vector <double>& x));

      void GetExtrema (std::vector <double>& xMin, std::vector <double>& xMax);

    private:
      std::vector < std::vector < double > > _line;
      std::vector < Marker*> _particles;
      std::vector < unsigned > _markerOffset;
      std::vector < unsigned > _printList;
      unsigned _size;
      unsigned _dim;

      void Reorder (std::vector < Marker*> &particles);

      static const double _a[4][4][4];
      static const double _b[4][4];
      static const double _c[4][4];

      std::vector< double > _time;
      const Mesh* _mesh;
      const Solution* _sol;

    protected:
      std::vector< unsigned int >::reference sol (const char* arg1);
  };
} //end namespace femus



#endif




