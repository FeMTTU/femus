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




  void Line::AdvectionParallel(Solution* sol, const unsigned &n, const double& T, const unsigned &order) {

//     //BEGIN  Initialize the parameters for all processors, to be used when awake
// 
//     vector < unsigned > solVIndex(_dim);
//     solVIndex[0] = sol->GetIndex("U");    // get the position of "U" in the ml_sol object
//     solVIndex[1] = sol->GetIndex("V");    // get the position of "V" in the ml_sol object
//     if(_dim == 3) solVIndex[2] = sol->GetIndex("W");       // get the position of "V" in the ml_sol object
//     unsigned solVType = sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
// 
//     std::vector < double > phi;
//     std::vector < std::vector<double > > V(2);
//     std::vector < std::vector < std::vector < double > > > aV;
// 
//     std::vector < std::vector < double > > x;
//     std::vector < std::vector < double > > x0;
//     x.resize(_size);
//     x0.resize(_size);
//     for(unsigned j = 0; j < _size; j++) {
//       x[j].assign(_dim, 0.);
//       x0[j].assign(_dim, 0.);
//     }
// 
//     std::vector <unsigned> step;
//     step.assign(_size, 0);
//     //END
// 
//     //BEGIN Numerical integration scheme
//     // determine the step size
//     double h = T / n;
// 
//     for(unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {
// 
//       unsigned elem = _particles[iMarker]->GetMarkerElement();
//       _particles[iMarker]->GetMarkerCoordinates(x[iMarker]);
//       _particles[iMarker]->GetMarkerx0(x0[iMarker]);
// 
//       bool integrationIsOver = (elem != UINT_MAX) ? false : true;
// 
//       _K.resize(_size);
//       for(unsigned k = 0; k < _size; k++) {
//         _K[k].resize(order);
//         for(unsigned j = 0; j < order; j++) {
//           _K[k][j].resize(_dim);
//         }
//       }
// 
//       while(integrationIsOver == false) {
//         unsigned mprocOld = _iproc;
// 
//         bool pcElemUpdate = true ;
//         while(step[iMarker] < n * order) {
// 
//           unsigned tstep = step[iMarker] / order;
//           unsigned istep = step[iMarker] % order;
// 
//           if(istep == 0) {
//             x0[iMarker] = x[iMarker];
//             for(unsigned k = 0; k < order; k++) {
//               _K[iMarker][k].assign(_dim, 0.);
//             }
//           }
// 
//           //  std::cout << " -----------------------------" << "step = " <<  step << " tstep = " << tstep << " istep = " << istep << " -----------------------------" << std::endl;
//           //  std::cout << " _iproc = " << _iproc << std::endl;
// 
//           Marker::updateVelocity(V, sol, solVIndex, solVType, aV, phi, pcElemUpdate); //send xk
// 
//           double s = (tstep + _c[order - 1][istep]) / n;
// 
//           for(unsigned k = 0; k < _dim; k++) {
//             _K[iMarker][istep][k] = (s * V[0][k] + (1. - s) * V[1][k]) * h;
//           }
// 
//           step[iMarker] += 1 ;
//           istep++;
// 
//           if(istep < order) {
// 
//             for(unsigned k = 0; k < _dim; k++) {
// 	      x[iMarker][k] = x0[iMarker][k];
//               for(unsigned j = 0; j < order; j++) {
//                 x[iMarker][k] +=  _a[order - 1][istep][j] * _K[iMarker][j][k];
//               }
//             }
//           }
//           
//           //WARNING arrivato qui
//           
//           else if(istep == order) {
// 
//             for(unsigned i = 0; i < _dim; i++) {
//               _x[i] = _x0[i];
//               for(unsigned j = 0; j < order; j++) {
//                 _x[i] += _b[order - 1][j] * _K[j][i];
//               }
//             }
//           }
// 
// 
//           //BEGIN TO BE REMOVED
//           for(unsigned i = 0; i < _dim; i++) {
//             //    std::cout << "_x[" << i << "] = " << _x[i] ;
//             //    std::cout << " " ;
//           }
//           //  std::cout << std::endl;
//           //END TO BE REMOVED
// 
//           pcElemUpdate = false;
//           unsigned iel = _elem;
//           previousElem = _elem;
//           Marker::GetElementSerial(previousElem);
// 
//           if(_elem == UINT_MAX) { //out of the domain
//             //    std::cout << " the marker has been advected outside the domain " << std::endl;
//             break;
//           }
//           else if(iel != _elem && _iproc != _mproc) { //different element different process
//             break;
//           }
//           else if(iel != _elem) { //different element same process
//             pcElemUpdate = true;
//             Marker::FindLocalCoordinates(solVType, _aX, pcElemUpdate);
//           }
//           else { //same element same process
//             Marker::FindLocalCoordinates(solVType, _aX, pcElemUpdate);
//           }
//         }
// 
// 
//       }
// 
// 
//     }


  }






}
