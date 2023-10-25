/*=========================================================================

 Program: FEMuS
 Module: Marker
 Authors: Eugenio Aulisa and Giacomo Capodaglio

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_ism_Marker_hpp__
#define __femus_ism_Marker_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MarkerTypeEnum.hpp"
#include "ParallelObject.hpp"
#include "NumericVector.hpp"
#include "MyVector.hpp"
#include "Mesh.hpp"

#include <vector>
#include <map>


namespace femus {

  class Marker : public ParallelObject {
    public:
      Marker (std::vector < double > x, const double &mass, const MarkerType &markerType, Solution *sol, const unsigned & solType, const bool &debug = false, const double &s1 = 0.) {
        //double s1 = 0.;
        _x = x;
        _markerType = markerType;
        _solType = solType;
        _dim = sol->GetMesh()->GetDimension();
        _step = 0;

        _MPMSize = 3 * _dim + 1; //removed density
        GetElement (1, UINT_MAX, sol, s1);

        if (_iproc == _mproc) {
          std::vector < std::vector < std::vector < std::vector < double > > > >aX;
          FindLocalCoordinates (_solType, aX, true, sol, s1);

          _MPMQuantities.resize (_MPMSize);
          _MPMQuantities[3 * _dim ] = mass; /*11.133 for the disk */ /*0.217013888889 for the beam */ ;  //mass //now it is computed in the main, zero is a default value
//           _MPMQuantities[3 * _dim + 1] = 10000.;  //density
          for (unsigned d = 0; d < 3 * _dim; d++) {
            _MPMQuantities[d] = 0.;
          }

          //unitialization of the deformation gradient of the particle to the identity matrix
          _Fp.resize (_dim);
          for (unsigned i = 0; i < _dim; i++) {
            _Fp[i].resize (_dim);
          }

          for (unsigned i = 0; i < _dim; i++) {
            for (unsigned j = 0; j < _dim; j++) {
              if (i == j) {
                _Fp[i][i] = 1.;
              }
              else {
                _Fp[i][j] = 0.;
              }
            }
          }

        }
        else {
          std::vector < double > ().swap (_x);
        }
      };

      double GetCoordinates (Solution *sol, const unsigned &k, const unsigned &i , const double &s) {
        if (!sol->GetIfFSI()) {
          return (*sol->GetMesh()->_topology->_Sol[k]) (i);
        }
        else {
          const char varname[3][3] = {"DX", "DY", "DZ"};
          unsigned solIndex = sol->GetIndex (&varname[k][0]);
          return (*sol->GetMesh()->_topology->_Sol[k]) (i)
                 + (1. - s) * (*sol->_SolOld[solIndex]) (i)
                 + s * (*sol->_Sol[solIndex]) (i);
        }
      }


      void SetIprocMarkerOldCoordinates (const std::vector <double> &x0) {
        _x0 = x0;
      }

      void SetIprocMarkerCoordinates (const std::vector <double> &x) {
        _x = x;
      }

      void SetIprocMarkerK (const std::vector < std::vector < double > > &K) {
        _K = K;
      }

      void SetIprocMarkerStep (const unsigned &step) {
        _step = step;
      }

      void SetMarkerElement (const unsigned &elem) {
        _elem = elem;
      }

      void SetMarkerProc (const unsigned &mproc) {
        _mproc = mproc;
      }



      void SetIprocMarkerPreviousElement (const unsigned &previousElem) {
        _previousElem = previousElem;
      }

      void GetNumberOfMeshElements (unsigned &elements, Solution *sol) {
        elements = sol->GetMesh()->GetNumberOfElements();
      }


      unsigned GetMarkerProc (Solution *sol) {
        _mproc = (_elem == UINT_MAX) ? 0 : sol->GetMesh()->IsdomBisectionSearch (_elem , 3);
        return _mproc;
      }

      void GetMarkerLocalCoordinates (std::vector< double > &xi) {

        xi.resize (_dim);
        if (_mproc == _iproc) {
          xi = _xi;
        }
        MPI_Bcast (&xi[0], _dim, MPI_DOUBLE, _mproc, PETSC_COMM_WORLD);
      }

//       void GetMarker_x0Line(std::vector<double> &x0) {
//  x0.resize(_dim);
//         x0 = _x0;
//       }


//      void GetMarker_KLine(std::vector < std::vector < double > > &K){
//  K = _K;
//       }

//       void GetMarkerStepLine(unsigned &step){
//  step = _step;
//       }

      unsigned GetMarkerElement() {
        return _elem;
      }

      void GetMarkerElementLine (unsigned &elem) {
        elem = _elem;
      }

      void SetMarkerMass (const double &mass) {
        _MPMQuantities[3 * _dim] = mass;
      }

//       void SetMarkerDensity(const double &density) {
//         _MPMQuantities[3 * _dim + 1] = density;
//       }

      double GetMarkerMass() {
        return _MPMQuantities[3 * _dim];
      }

//       double GetMarkerDensity() {
//         return _MPMQuantities[3 * _dim + 1];
//       }

      void SetMarkerVelocity (const std::vector <double>  &velocity) {
        if (_markerType != INTERFACE) {
          for (unsigned d = 0; d < _dim; d++) {
            _MPMQuantities[_dim + d] = velocity[d];
          }
        }
        else {
          std::cout <<  "Wrong markerType, velocity is not available for INTERFACE markers!\n";
          abort();
        }
      }

      void SetMarkerTangentGlobal (const std::vector < std::vector <double> >  &tangent) {
        if(_iproc == _mproc){
          SetMarkerTangent (tangent);
        }  
      }
      
      void SetMarkerTangent (const std::vector < std::vector <double> >  &tangent) {
        if (_markerType == INTERFACE) {
          for (unsigned k = 0; k < _dim - 1; k++) {
            for (unsigned d = 0; d < _dim; d++) {
              _MPMQuantities[_dim * (k + 1) + d] = tangent[k][d];
            }
          }
        }
        else {
          std::cout <<  "Wrong markerType, tengent is available only for INTERFACE markers!\n";
          abort();
        }
      }

      void GetMarkerVelocity (std::vector <double> & velocity) {
        if (_markerType != INTERFACE) {
          velocity.resize (_dim);
          for (unsigned d = 0; d < _dim; d++) {
            velocity[d] = _MPMQuantities[_dim + d];
          }
        }
        else {
          std::cout <<  "Wrong markerType, velocity is not available for INTERFACE markers!\n";
          abort();
        }
      }


      void GetMarkerTangent (std::vector < std::vector <double> > & tangent) {
        if (_markerType == INTERFACE) {
          tangent.resize (_dim - 1);
          for (unsigned k = 0; k < _dim - 1; k++) {
            tangent[k].resize (_dim);
            for (unsigned d = 0; d < _dim; d++) {
              tangent[k][d] = _MPMQuantities[_dim * (k + 1) + d];
            }
          }
        }
        else {
          std::cout <<  "Wrong markerType, tangent is available only for INTERFACE markers!\n";
          abort();
        }
      }

      void UpdateParticleVelocities (const std::vector <double> newParticleAcceleration, const double dt) {
        if (_markerType != INTERFACE) {
          for (unsigned d = 0; d < _dim; d++) {
            _MPMQuantities[_dim + d] += 0.5 * (_MPMQuantities[2 * _dim + d] + newParticleAcceleration[d]) * dt ;
          }
        }
        else {
          std::cout <<  "Wrong MarkerType, velocity is not available for INTERFACE markers!\n";
          abort();
        }
      }

      void SetMarkerAcceleration (const std::vector <double>  &acceleration) {
        if (_markerType != INTERFACE) {
          for (unsigned d = 0; d < _dim; d++) {
            _MPMQuantities[2 * _dim + d] = acceleration[d];
          }
        }
        else {
          std::cout <<  "Wrong MarkerType, acceleration is not available for INTERFACE markers!\n";
          abort();
        }
      }

      void GetMarkerAcceleration (std::vector <double> & acceleration) {
        if (_markerType != INTERFACE) {
          acceleration.resize (_dim);
          for (unsigned d = 0; d < _dim; d++) {
            acceleration[d] = _MPMQuantities[2 * _dim + d];
          }
        }
        else {
          std::cout <<  "Wrong MarkerType, acceleration is not available for INTERFACE markers!\n";
          abort();
        }
      }

      void SetMarkerDisplacement (const std::vector <double>  &displacement) {
        for (unsigned d = 0; d < displacement.size(); d++) {
          _MPMQuantities[d] = displacement[d];
        }
      }

      void GetMarkerDisplacement (std::vector <double> & displacement) {
        for (unsigned d = 0; d < displacement.size(); d++) {
          displacement[d] = _MPMQuantities[d];
        }
      }

      void UpdateParticleCoordinates() {
        for (unsigned i = 0; i < _dim; i++) {
          _x[i] += _MPMQuantities[i];
        }
      }

      void SetMPMQuantities (const std::vector <double>  &MPMQuantities) {
        _MPMQuantities = MPMQuantities;
      }

      void SetDeformationGradient (const std::vector < std::vector < double > > Fp) {
        _Fp = Fp;
      }

      std::vector <double> GetMPMQuantities() {
        return _MPMQuantities;
      }

      unsigned GetMPMSize() {
        return _MPMSize;
      }

      std::vector < std::vector < double > > GetDeformationGradient() {
        return _Fp;
      }

      std::vector<double> GetMarkerLocalCoordinates() {
        return _xi;
      }

      void GetMarkerLocalCoordinatesLine (std::vector<double> &xi) {
        xi.resize (_dim);
        if (_mproc == _iproc) {
          xi = _xi;
        }
        MPI_Bcast (&xi[0], _dim, MPI_DOUBLE, _mproc, PETSC_COMM_WORLD);
      }

      std::vector< double > GetIprocMarkerCoordinates() {
        return _x;
      }

      std::vector< double > GetIprocMarkerOldCoordinates() {
        return _x0;
      }

      unsigned GetIprocMarkerStep() {
        return _step;
      }

      std::vector< std::vector < double> > GetIprocMarkerK() {
        return _K;
      }

      void GetMarkerCoordinates (std::vector< double > &xn) {
        xn.resize (_dim);
        if (_mproc == _iproc) {
          xn = _x;
        }
        MPI_Bcast (&xn[0], _dim, MPI_DOUBLE, _mproc, PETSC_COMM_WORLD);
      }

      void GetMarkerCoordinates (std::vector< MyVector <double > > &xn) {
        if (_mproc == _iproc) {
          for (unsigned d = 0; d < _dim; d++) {
            unsigned size = xn[d].size();
            xn[d].resize (size + 1);
            xn[d][size] = _x[d];
          }
        }
      }

      unsigned GetIprocMarkerPreviousElement() {
        return _previousElem;
      }

      void InitializeMarkerForAdvection (const unsigned & order) {
        _x0 = _x;
        _step = 0;
        _K.resize (order);
        for (unsigned j = 0; j < order; j++) {
          _K[j].assign (_dim, 0.);
        }
      }

      void InitializeVariables (const unsigned & order) {
        _xi.resize (_dim);
        _x0.resize (_dim);
        _K.resize (order);
        for (unsigned j = 0; j < order; j++) {
          _K[j].assign (_dim, 0.);
        }
        _MPMQuantities.resize (3 * _dim + 2);
        _Fp.resize (_dim);
        for (unsigned i = 0; i < _dim; i++) {
          _Fp[i].resize (_dim);
        }
      }

      void InitializeX() {
        _x.resize (_dim);
      }

      void FreeVariables() {
        std::vector < double > ().swap (_xi);
        std::vector < double > ().swap (_x0);
        std::vector < std::vector < double > > ().swap (_K);
        std::vector < double > ().swap (_MPMQuantities);
        std::vector <  std::vector < double > > ().swap (_Fp);
      }

      void GetElement (const bool & useInitialSearch, const unsigned & initialElem, Solution * sol, const double & s);
      void GetElementSerial (unsigned & initialElem, Solution * sol, const double & s);
      void GetElement (unsigned & previousElem, const unsigned & previousMproc, Solution * sol, const double & s);


      MarkerType GetMarkerType() {
        return _markerType;
      };

      void InverseMappingTEST (std::vector< double > &x, Solution * sol, const double & s);
      void Advection (const unsigned & n, const double & T, Solution * sol);

      void updateVelocity (std::vector< std::vector <double> > & V,
                           const std::vector < unsigned > &solVIndex, const unsigned & solVType,
                           std::vector < std::vector < std::vector < double > > > &a,  std::vector < double > &phi,
                           const bool & pcElemUpdate, Solution * sol);

      void FindLocalCoordinates (const unsigned & solVType, std::vector < std::vector < std::vector < std::vector < double > > > > &aX,
                                 const bool & pcElemUpdate, Solution * sol, const double & s);

      void ProjectVelocityCoefficients (const std::vector<unsigned> &solVIndex,
                                        const unsigned & solVType,  const unsigned & nDofsV,
                                        const unsigned & ielType, std::vector < std::vector < std::vector < double > > > &a, Solution * sol);

      void GetMarkerS (const unsigned & n, const unsigned & order, double & s) {
        unsigned step = (_step == UINT_MAX) ? n * order : _step;
        unsigned tstep = step / order;
        unsigned istep = step % order;
        s = (tstep + _c[order - 1][istep]) / n;
      }


    private:


      std::vector< double > InverseMapping (const unsigned & currentElem, const unsigned & solutionType, const std::vector< double > &x);
      void InverseMapping (const unsigned & iel, const unsigned & solType,
                           const std::vector< double > &x, std::vector< double > &xi, Solution * sol, const double & s);

      unsigned GetNextElement2D (const unsigned & iel, const unsigned & previousElem, Solution * sol, const double & s);
      unsigned GetNextElement3D (const unsigned & iel, const unsigned & previousElem, Solution * sol, const double & s);
      unsigned GetNextElement3D (const unsigned & iel, const std::vector< unsigned > &searchHistory, Solution * sol, const double & s);
      int FastForward (const unsigned & currentElem, const unsigned & previousElem, Solution * sol, const double & s);

      std::vector < double > _x;
      std::vector < double > _x0;
      std::vector < double > _xi;
      unsigned _solType;
      MarkerType _markerType;
      //const Mesh * mesh;
      unsigned _elem;
      unsigned _previousElem; //for advection reasons
      unsigned _dim;

      unsigned _mproc; //processor who has the marker
      std::vector < std::vector < double > > _K;
      unsigned _step; //added for line

      static const double _localCentralNode[6][3];
      static const double _a[4][4][4];
      static const double _b[4][4];
      static const double _c[4][4];

      std::vector <double> _MPMQuantities; // _displacement[_dim] + _velocity[_dim] + _acceleration[_dim] + mass /* + density */
      unsigned _MPMSize;

      std::vector < std::vector <double> > _Fp; //deformation gradient of the particle


  };
} //end namespace femus



#endif


