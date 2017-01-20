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

    //BEGIN  Initialize the parameters for all processors, to be used when awake

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


    //BEGIN Numerical integration scheme

    //loop on the markers that belong to _iproc
    for(unsigned i = _markerOffset[_iproc]; i < _markerOffset[_iproc + 1]; i++) {

      std::vector< double > x;
      
      _particles[i]->GetMarkerCoordinates(x);
      unsigned elem = _particles[i]->GetMarkerElement();

      bool integrationIsOver = (elem != UINT_MAX) ? false : true;

      unsigned step = 0.;

      _K.resize(order);
      for(unsigned k = 0; k < order; k++) {
        _K[k].resize(_dim);
      }


      while(integrationIsOver == false) {
        unsigned mprocOld = _iproc;
        unsigned previousElem;


        bool pcElemUpdate = true ;
        while(step < n * order) {

          unsigned tstep = step / order;
          unsigned istep = step % order;

          if(istep == 0) {
            _particles[i]->SetMarkerx0(x);
            for(unsigned k = 0; k < order; k++) {
              _K[k].assign(_dim, 0.);
            }
          }

          //  std::cout << " -----------------------------" << "step = " <<  step << " tstep = " << tstep << " istep = " << istep << " -----------------------------" << std::endl;
          //  std::cout << " _iproc = " << _iproc << std::endl;

          updateVelocity(V, sol, solVIndex, solVType, aV, phi, pcElemUpdate); //send xk

          double s = (tstep + _c[order - 1][istep]) / n;

          for(unsigned i = 0; i < _dim; i++) {
            _K[istep][i] = (s * V[0][i] + (1. - s) * V[1][i]) * h;
          }

          step++;
          istep++;

          if(istep < order) {

            for(unsigned i = 0; i < _dim; i++) {
              _x[i] = _x0[i];
              for(unsigned j = 0; j < order; j++) {
                _x[i] +=  _a[order - 1][istep][j] * _K[j][i];
              }
            }
          }
          else if(istep == order) {

            for(unsigned i = 0; i < _dim; i++) {
              _x[i] = _x0[i];
              for(unsigned j = 0; j < order; j++) {
                _x[i] += _b[order - 1][j] * _K[j][i];
              }
            }
          }


          pcElemUpdate = false;
          unsigned iel = _elem;
          previousElem = _elem;
          GetElementSerial(previousElem);

          if(_elem == UINT_MAX) { //out of the domain
            //    std::cout << " the marker has been advected outside the domain " << std::endl;
            break;
          }
          else if(iel != _elem && _iproc != _mproc) { //different element different process
            break;
          }
          else if(iel != _elem) { //different element same process
            pcElemUpdate = true;
            FindLocalCoordinates(solVType, _aX, pcElemUpdate);
          }
          else { //same element same process
            FindLocalCoordinates(solVType, _aX, pcElemUpdate);
          }
        }



      }


    }


  }
  
  void Marker::GetElement(unsigned &previousElem, const unsigned &previousMproc) {

    unsigned mprocOld = previousMproc;

    if(mprocOld == _iproc) {
      MPI_Send(&previousElem, 1, MPI_UNSIGNED, _mproc, 1 , PETSC_COMM_WORLD);
      MPI_Send(&_x[0], _dim, MPI_DOUBLE, _mproc, 2 , PETSC_COMM_WORLD);
      std::vector < double > ().swap(_x);
    }
    else if(_mproc == _iproc) {
      MPI_Recv(&previousElem, 1, MPI_UNSIGNED, mprocOld, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      _x.resize(_dim);
      MPI_Recv(&_x[0], _dim, MPI_DOUBLE, mprocOld, 2 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    mprocOld = _mproc;
    while(true) {     //Corretto ora dovrebbe funzionare. WARNING credo debba essere while(!found) perche' senno' non fa mai GetElementSerial e non cambia mai 716 in 724, pero' se si mette !found non ranna piu niente
      if(_mproc == _iproc) {
        GetElementSerial(previousElem);
      }
      MPI_Bcast(& _elem, 1, MPI_UNSIGNED, mprocOld, PETSC_COMM_WORLD);
      if(_elem == UINT_MAX) {
        //   std::cout << " the marker has been advected outside the domain " << std::endl;
        break;
      }
      else {
        _mproc = _mesh->IsdomBisectionSearch(_elem, 3);
        if(_mproc == mprocOld) {
          break;
        }
        else {
          if(mprocOld == _iproc) {
            MPI_Send(&previousElem, 1, MPI_UNSIGNED, _mproc, 1 , PETSC_COMM_WORLD);
            MPI_Send(&_x[0], _dim, MPI_DOUBLE, _mproc, 2 , PETSC_COMM_WORLD);
            std::vector < double > ().swap(_x);
          }
          else if(_mproc == _iproc) {
            MPI_Recv(&previousElem, 1, MPI_UNSIGNED, mprocOld, 1 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
            _x.resize(_dim);
            MPI_Recv(&_x[0], _dim, MPI_DOUBLE, mprocOld, 2 , PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          }
          mprocOld = _mproc;
        }
      }
    }
  }



  void Marker::ProjectVelocityCoefficients(Solution * sol, const std::vector<unsigned> &solVIndex,
      const unsigned & solVType,  const unsigned & nDofsV,
      const unsigned & ielType, std::vector < std::vector < std::vector < double > > > &a) {

    bool timeDependent = true;
    for(unsigned  k = 0; k < _dim; k++) {
      if(sol->GetSolutionTimeOrder(solVIndex[k]) != 2) {
        timeDependent = false;
      }
    }

    vector < vector < double > >  solV(_dim);    // local solution
    vector < vector < double > >  solVold(_dim);    // local solution

    for(unsigned  k = 0; k < _dim; k++) {
      solV[k].resize(nDofsV);
      if(timeDependent) {
        solVold[k].resize(nDofsV);
      }
    }
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = _mesh->GetSolutionDof(i, _elem, solVType);    // global to global mapping between solution node and solution dof
      for(unsigned  k = 0; k < _dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        if(timeDependent) {
          solVold[k][i] = (*sol->_SolOld[solVIndex[k]])(solVDof);
        }
      }
    }

    a.resize(2);

    ProjectNodalToPolynomialCoefficients(a[0], solV, ielType, solVType);
    if(timeDependent) {
      ProjectNodalToPolynomialCoefficients(a[1], solVold, ielType, solVType);
    }
    else {
      a[1] = a[0];
    }
  }


  void Marker::updateVelocity(std::vector< std::vector <double> > & V, Solution * sol,
                              const vector < unsigned > &solVIndex, const unsigned & solVType,
                              std::vector < std::vector < std::vector < double > > > &a,  std::vector < double > &phi,
                              const bool & pcElemUpdate) {


    unsigned nDofsV = _mesh->GetElementDofNumber(_elem, solVType);
    short unsigned ielType = _mesh->GetElementType(_elem);

    if(pcElemUpdate) {
      ProjectVelocityCoefficients(sol, solVIndex, solVType, nDofsV, ielType, a);
    }

    GetPolynomialShapeFunction(phi, _xi, ielType, solVType);

    V.resize(2);
    V[0].assign(_dim, 0.);
    V[1].assign(_dim, 0.);
    for(unsigned i = 0; i < _dim; i++) {
      for(unsigned j = 0; j < nDofsV; j++) {
        V[0][i] += a[0][i][j] * phi[j];
        V[1][i] += a[1][i][j] * phi[j];
      }
    }


  }

  void Marker::FindLocalCoordinates(const unsigned & solType, std::vector < std::vector < std::vector < double > > > &aX, const bool & pcElemUpdate) {

    //BEGIN TO BE REMOVED
   // std::cout << "ENTRIAMO NELL' INVERSE MAPPING, _elem =" << _elem << std::endl;
    for(unsigned i = 0; i < 3; i++) {
    //  std::cout << "_x[" << i << "]= " << _x[i] << std::endl;
    }
    //END TO BE REMOVED


    short unsigned elemType = _mesh->GetElementType(_elem);
    if(pcElemUpdate) {
      unsigned nDofs = _mesh->GetElementDofNumber(_elem, solType);

      //BEGIN extraction nodal coordinate values
      std::vector< std::vector < double > > xv(_dim);

      for(unsigned k = 0; k < _dim; k++) {
        xv[k].resize(nDofs);
      }

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned iDof  = _mesh->GetSolutionDof(i, _elem, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < _dim; k++) {
          xv[k][i] = (*_mesh->_topology->_Sol[k])(iDof);     // global extraction and local storage for the element coordinates
        }
      }
      //END extraction

      //BEGIN projection nodal to polynomial coefficients
      aX.resize(solType + 1);
      for(unsigned j = 0; j < solType + 1; j++) {
        ProjectNodalToPolynomialCoefficients(aX[j], xv, elemType, j);
      }
      //END projection nodal to polynomial coefficients

      //BEGIN find initial guess

      GetClosestPointInReferenceElement(xv, _x, elemType, _xi);

      //BEGIN TO BE REMOVED
     // std::cout << "INITIAL GUESS" << std::endl;
      for(unsigned i = 0; i < 3; i++) {
      //  std::cout << "_xi[" << i << "]= " << _xi[i] << std::endl;
      }
      //END TO BE REMOVED


      /*
      _xi.resize(_dim);
      for(int k = 0; k < _mesh->GetDimension(); k++) {
        _xi[k] = _localCentralNode[elemType][k];
      }*/
      //END find initial guess
    }

    //BEGIN Inverse mapping loop
    for(unsigned j = 0; j < solType; j++) {
      std::vector < double > phi;
      std::vector < std::vector < double > > gradPhi;
      bool convergence = false;
      while(!convergence) {
        GetPolynomialShapeFunctionGradient(phi, gradPhi, _xi, elemType, solType);
        convergence = GetNewLocalCoordinates(_xi, _x, phi, gradPhi, aX[solType]);
      }
    }
    //END Inverse mapping loop


  //  std::cout << "ED USCIAMO" << std::endl;

  }
  
  
  
}
