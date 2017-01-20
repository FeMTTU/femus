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

   void Line::UpdateParticles() {

    std::vector < Marker*> particlesOld(_size);
    particlesOld = _particles;

    unsigned counter = 0;
    for(unsigned iproc = 0; iproc < _nprocs; iproc++) {
      for(unsigned j = 0; j < _size; j++) {
        unsigned markerProc = particlesOld[j]->GetMarkerProc();
        if(markerProc == iproc) {
          _particles[counter] = particlesOld[j];
          counter++;
        }
      }
      _markerOffset[iproc] = counter;
    }

    std::vector < Marker* > ().swap(particlesOld);
  }

   
  
}
