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


  void Line::AdvectionParallel(Solution* sol, const unsigned &n, const double& T, const unsigned &order) {

    //BEGIN  Initialize the parameters for all processors, to be used when awake

    vector < unsigned > solVIndex(_dim);
    solVIndex[0] = sol->GetIndex("U");    // get the position of "U" in the ml_sol object
    solVIndex[1] = sol->GetIndex("V");    // get the position of "V" in the ml_sol object
    if(_dim == 3) solVIndex[2] = sol->GetIndex("W");       // get the position of "V" in the ml_sol object
    unsigned solVType = sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

    std::vector < double > phi;
    std::vector < std::vector<double > > V(2);
    std::vector < std::vector < std::vector < double > > > aV;


    std::vector <unsigned> step;
    step.assign(_size, 0);

    std::vector < double > x;
    std::vector < double > x0;
    //END

    //BEGIN Numerical integration scheme
    // determine the step size
    double h = T / n;

    unsigned previousElem = UINT_MAX;

    for(unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {

      _particles[iMarker]->GetMarkerCoordinates(x);
      x0 = _particles[iMarker]->GetMarkerx0();

      unsigned currentElem = _particles[iMarker]->GetMarkerElement();

      bool integrationIsOver = (currentElem != UINT_MAX) ? false : true;

      _K.resize(_size);
      for(unsigned k = 0; k < _size; k++) {
        _K[k].resize(order);
        for(unsigned j = 0; j < order; j++) {
          _K[k][j].resize(_dim);
        }
      }

      while(integrationIsOver == false) {
        unsigned mprocOld = _iproc; //TODO don't know if we need this

        bool pcElemUpdate = (previousElem == currentElem) ? false : true; //update only if the marker is in a different element

        while(step[iMarker] < n * order) {

          unsigned tstep = step[iMarker] / order;
          unsigned istep = step[iMarker] % order;

          if(istep == 0) {
            x0 = x;
            for(unsigned k = 0; k < order; k++) {
              _K[iMarker][k].assign(_dim, 0.);
            }
          }

          _particles[iMarker]->updateVelocity(V, sol, solVIndex, solVType, aV, phi, pcElemUpdate); //send xk

          double s = (tstep + _c[order - 1][istep]) / n;

          for(unsigned k = 0; k < _dim; k++) {
            _K[iMarker][istep][k] = (s * V[0][k] + (1. - s) * V[1][k]) * h;
          }

          step[iMarker] += 1 ;
          istep++;

          if(istep < order) {

            for(unsigned k = 0; k < _dim; k++) {
              x[k] = x0[k];
              for(unsigned j = 0; j < order; j++) {
                x[k] +=  _a[order - 1][istep][j] * _K[iMarker][j][k];
              }
            }
          }

          else if(istep == order) {

            for(unsigned i = 0; i < _dim; i++) {
              x[i] = x0[i];
              for(unsigned j = 0; j < order; j++) {
                x[i] += _b[order - 1][j] * _K[iMarker][j][i];
              }
            }
          }

          _particles[iMarker]->SetMarkerx0(x0);
          _particles[iMarker]->SetMarkerCoordinates(x);

          previousElem = currentElem;

          _particles[iMarker]->GetElementSerial(currentElem);

          currentElem = _particles[iMarker]->GetMarkerElement();

          if(currentElem == UINT_MAX) { //out of the domain
            //    std::cout << " the marker has been advected outside the domain " << std::endl;
            break;
          }
          else if(previousElem != currentElem && _iproc != _particles[iMarker]->GetMarkerProc()) { //different element different process
            break;
          }
          else if(previousElem != currentElem) { //different element same process
            pcElemUpdate = true;
            _particles[iMarker]->FindLocalCoordinates(solVType, _aX[iMarker], pcElemUpdate);
          }
          else { //same element same process
            _particles[iMarker]->FindLocalCoordinates(solVType, _aX[iMarker], pcElemUpdate);
          }
        }
      }
    }


    //TODO TO BE CONTINUED...


  }
}
