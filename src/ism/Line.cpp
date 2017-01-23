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


  void Line::UpdateLine(){
    
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
    double h = T / n;

    //END


    //BEGIN initialize marker instances
    unsigned step;
    std::vector < double > x;
    std::vector < double > x0;
    std::vector < std::vector < double > > K;
    std::vector < std::vector < std::vector < double > > > aX;
    //END

    unsigned integrationIsOverCounter = 0; // when integrationIsOverCounter = _size (which is the number of particles) it means all particles have been advected;
    // when the step of each marker is = n * order and/or marker element is UINT_MAX we increase by one such counter

    //BEGIN Numerical integration scheme

    while(integrationIsOverCounter != _size) {

      unsigned previousElem = UINT_MAX;

      for(unsigned iMarker = _markerOffset[_iproc]; iMarker < _markerOffset[_iproc + 1]; iMarker++) {


        //BEGIN extraction of the marker instances
        x.resize(_dim);
        x0.resize(_dim);
        _particles[iMarker]->GetMarkerCoordinates(x);
        x0 = _particles[iMarker]->GetMarker_x0();

        K.resize(order);
        for(unsigned j = 0; j < order; j++) {
          K[j].resize(_dim);
        }

        K = _particles[iMarker]->GetMarker_K();
        step = _particles[iMarker]->GetMarkerStep();
        aX = _particles[iMarker]->GetMarker_aX();


        unsigned currentElem = _particles[iMarker]->GetMarkerElement();

        //END

        bool markerOutsideDomain = (currentElem != UINT_MAX) ? false : true;

        if(markerOutsideDomain == false) {
          unsigned mprocOld = _iproc; //TODO don't know if we need this

          bool pcElemUpdate = (previousElem == currentElem) ? false : true; //update only if the marker is in a different element

          while(step < n * order) {

            unsigned tstep = step / order;
            unsigned istep = step % order;

            if(istep == 0) {
              x0 = x;
              for(unsigned k = 0; k < order; k++) {
                K[k].assign(_dim, 0.);
              }
            }

            _particles[iMarker]->updateVelocity(V, sol, solVIndex, solVType, aV, phi, pcElemUpdate); //send xk

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

            _particles[iMarker]->SetMarker_x0(x0);
            _particles[iMarker]->SetMarkerCoordinates(x);
            _particles[iMarker]->SetMarkerStep(step);
            _particles[iMarker]->SetMarker_K(K);

            previousElem = currentElem;

            _particles[iMarker]->GetElementSerial(currentElem);

            currentElem = _particles[iMarker]->GetMarkerElement();

            if(currentElem == UINT_MAX) { //out of the domain
              //    std::cout << " the marker has been advected outside the domain " << std::endl;
              break; //TODO double check where this is taking us
            }
            else if(previousElem != currentElem && _iproc != _particles[iMarker]->GetMarkerProc()) { //different element different process
              break; //TODO double check where this is taking us
            }
            else if(previousElem != currentElem) { //different element same process
              pcElemUpdate = true;
              _particles[iMarker]->FindLocalCoordinates(solVType, aX, pcElemUpdate);
            }
            else { //same element same process
              _particles[iMarker]->FindLocalCoordinates(solVType, aX, pcElemUpdate);
            }
          }

          if(step == n * order) {
            integrationIsOverCounter++;
            step = UINT_MAX;
            _particles[iMarker]->SetMarkerStep(step);
            //     std::cout << "Integration is over, point in proc " << _mproc << std::endl;
          }

          
          _particles[iMarker]->GetElement(previousElem, mprocOld);
	  
	  currentElem = _particles[iMarker]->GetMarkerElement();
	  
          
          
          //BEGIN exchange of information to reorder _particles

          for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
	    MPI_Send(&integrationIsOverCounter, 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD);
            MPI_Send(&x0, _dim, MPI_DOUBLE, jproc, 1 , PETSC_COMM_WORLD);
            MPI_Send(&x, _dim, MPI_DOUBLE, jproc, 1 , PETSC_COMM_WORLD);
            MPI_Send(&step, 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD);
            for(unsigned j = 0; j < order; j++) {
              MPI_Send(&K[j], _dim, MPI_DOUBLE, jproc, 1 , PETSC_COMM_WORLD);
            }
            MPI_Send(&currentElem, 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD);
            for(unsigned i = 0; i < solVType + 1 ; i++) {
              for(unsigned j = 0; j < _dim; j++) {
                MPI_Send(&aX[i][j], 27, MPI_DOUBLE, jproc, 1 , PETSC_COMM_WORLD); //TODO qua dovremmo mandare solo tanti quanti il numero di dofs dell'elemento, non 27 che e' il  max
              }
            }

            if(jproc != _iproc) {
              MPI_Recv(&x0, _dim, MPI_DOUBLE, _iproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
              MPI_Recv(&x, _dim, MPI_DOUBLE, _iproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
              MPI_Recv(&step, 1, MPI_UNSIGNED, _iproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
	      MPI_Recv(&integrationIsOverCounter, 1, MPI_UNSIGNED, _iproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
              for(unsigned j = 0; j < order; j++) {
                MPI_Recv(&K[j], _dim, MPI_DOUBLE, _iproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
              }
              for(unsigned i = 0; i < solVType + 1 ; i++) {
                for(unsigned j = 0; j < _dim; j++) {
                  MPI_Recv(&aX[i][j], 27, MPI_DOUBLE, _iproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                }
              }
              MPI_Recv(&currentElem, 1, MPI_UNSIGNED, _iproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

              _particles[iMarker]->SetMarker_x0(x0);
              _particles[iMarker]->SetMarkerCoordinates(x);
              _particles[iMarker]->SetMarkerStep(step);
              _particles[iMarker]->SetMarker_K(K);
              _particles[iMarker]->SetMarkerElement(currentElem);
              _particles[iMarker]->SetMarker_aX(aX);
            }
          }

          //END exchange of information


        }

        if(step != UINT_MAX) {
          integrationIsOverCounter++;
          step = UINT_MAX;
          _particles[iMarker]->SetMarkerStep(step);
          for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
            MPI_Send(&step, 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD);
	    MPI_Send(&integrationIsOverCounter, 1, MPI_UNSIGNED, jproc, 1 , PETSC_COMM_WORLD);
            if(jproc != _iproc) {
              MPI_Recv(&step, 1, MPI_UNSIGNED, _iproc, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
              _particles[iMarker]->SetMarkerStep(step);
            }
          }
        }
      }
      
      UpdateLine(); //TODO write this function, it should reorder _particles like the constructor so that all procs can restart advecting
      
    }

    //TODO TO BE CONTINUED...


  }
}

