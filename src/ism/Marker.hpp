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
#include "Mesh.hpp"

#include "vector"
#include "map"
#include "MyVector.hpp"

namespace femus {

  class Marker : public ParallelObject {
    public:
      Marker(std::vector < double > x, const MarkerType &markerType, Mesh *mesh, const unsigned & solType, const bool &debug = false) {
        _x = x;
        _markerType = markerType;
        _mesh = mesh;
        _solType = solType;
        _dim = _mesh->GetDimension();
        _step = 0;

        GetElement(1, UINT_MAX);

        if(_iproc == _mproc) {
          std::vector < std::vector < std::vector < double > > > aX;
          FindLocalCoordinates(_solType, aX, true);
        }
        else {
          std::vector < double > ().swap(_x);
        }
      };


      void SetIprocMarkerOldCoordinates(const std::vector <double> &x0) {
        _x0 = x0;
      }

      void SetIprocMarkerCoordinates(const std::vector <double> &x) {
        _x = x;
      }

      void SetIprocMarkerK(const std::vector < std::vector < double > > &K) {
        _K = K;
      }

      void SetIprocMarkerStep(const unsigned &step) {
        _step = step;
      }

      void SetMarkerElement(const unsigned &elem) {
        _elem = elem;
      }

      void SetMarkerProc(const unsigned &mproc) {
        _mproc = mproc;
      }


      void GetNumberOfMeshElements(unsigned &elements) {
        elements = _mesh->GetNumberOfElements();
      }

      unsigned GetMarkerProc() {
	_mproc = (_elem == UINT_MAX)? 0 : _mesh->IsdomBisectionSearch(_elem , 3);
        return _mproc;
      }

//       void GetMarker_x0Line(std::vector<double> &x0) {
// 	x0.resize(_dim);
//         x0 = _x0;
//       }


//      void GetMarker_KLine(std::vector < std::vector < double > > &K){
// 	K = _K;
//       }

//       void GetMarkerStepLine(unsigned &step){
// 	step = _step;
//       }

      unsigned GetMarkerElement() {
        return _elem;
      }

      void GetMarkerElementLine(unsigned &elem) {
        elem = _elem;
      }

      std::vector<double> GetMarkerLocalCoordinates() {
        return _xi;
      }

      void GetMarkerLocalCoordinatesLine(std::vector<double> &xi) {
        xi.resize(_dim);
        if(_mproc == _iproc) {
          xi = _xi;
        }
        MPI_Bcast(&xi[0], _dim, MPI_DOUBLE, _mproc, PETSC_COMM_WORLD);
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

      void InitializeMarkerForAdvection(const unsigned &order) {
        _x0 = _x;
        _step = 0;
        _K.resize(order);
        for(unsigned j = 0; j < order; j++) {
          _K[j].assign(_dim, 0.);
        }
      }

      void InitializeX0andK(const unsigned &order) {
	_xi.resize(_dim);
        _x0.resize(_dim);
        _K.resize(order);
        for(unsigned j = 0; j < order; j++) {
          _K[j].assign(_dim, 0.);
        }
      }

      void FreeXiX0andK() {
        std::vector < double > ().swap(_xi);
        std::vector < double > ().swap(_x0);
        std::vector < std::vector < double > > ().swap(_K);
      }


      void GetMarkerCoordinates(std::vector< double > &xn) {
        xn.resize(_dim);
        if(_mproc == _iproc) {
          xn = _x;
        }
        MPI_Bcast(&xn[0], _dim, MPI_DOUBLE, _mproc, PETSC_COMM_WORLD);
      }

      void GetMarkerCoordinates(std::vector< MyVector <double > > &xn) {
        if(_mproc == _iproc) {
          for(unsigned d = 0; d < _dim; d++) {
            unsigned size = xn[d].size();
            xn[d].resize(size + 1);
            xn[d][size] = _x[d];
          }
        }
      }

      unsigned GetIprocMarkerPreviousElement(){
	return _previousElem;
      }
      
      void SetIprocMarkerPreviousElement( const unsigned &previousElem){
	_previousElem = previousElem;
      }
      
      void GetElement(const bool &useInitialSearch, const unsigned &initialElem);
      void GetElementSerial(unsigned &initialElem);
      void GetElement(unsigned &previousElem, const unsigned &previousMproc);


      MarkerType GetMarkerType() {
        return _markerType;
      };

      void InverseMappingTEST(std::vector< double > &x);
      void Advection(Solution* sol, const unsigned &n, const double& T);

      void updateVelocity(std::vector< std::vector <double> > & V, Solution * sol,
                          const vector < unsigned > &solVIndex, const unsigned & solVType,
                          std::vector < std::vector < std::vector < double > > > &a,  std::vector < double > &phi,
                          const bool & pcElemUpdate);

      void FindLocalCoordinates(const unsigned & solVType, std::vector < std::vector < std::vector < double > > > &aX,
                                const bool & pcElemUpdate);

      void ProjectVelocityCoefficients(Solution * sol, const std::vector<unsigned> &solVIndex,
                                       const unsigned &solVType,  const unsigned &nDofsV,
                                       const unsigned &ielType, std::vector < std::vector < std::vector < double > > > &a);



    private:


      std::vector< double > InverseMapping(const unsigned &currentElem, const unsigned &solutionType, const std::vector< double > &x);
      void InverseMapping(const unsigned &iel, const unsigned &solType,
                          const std::vector< double > &x, std::vector< double > &xi);

      unsigned GetNextElement2D(const unsigned &iel, const unsigned &previousElem);
      unsigned GetNextElement3D(const unsigned &iel, const unsigned &previousElem);
      int FastForward(const unsigned &currentElem, const unsigned &previousElem);

      std::vector < double > _x;
      std::vector < double > _x0;
      std::vector < double > _xi;
      unsigned _solType;
      MarkerType _markerType;
      const Mesh * _mesh;
      unsigned _elem;
      unsigned _previousElem; //for advection reasons
      unsigned _dim;

      unsigned _mproc; //processor who has the marker
      //std::vector < std::vector < std::vector < double > > > _aX;
      std::vector < std::vector < double > > _K;
      unsigned _step; //added for line

      static const double _localCentralNode[6][3];
      static const double _a[4][4][4];
      static const double _b[4][4];
      static const double _c[4][4];

  };
} //end namespace femus



#endif

