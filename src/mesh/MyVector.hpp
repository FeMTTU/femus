
#ifndef __femus_mesh_MyVector_hpp__
#define __femus_mesh_MyVector_hpp__


#include <iostream>
#include <vector>
#include <mpi.h>
#include <stdlib.h>

namespace femus {

  template <class Type> class MyVector {
    public:
      // ******************
      MyVector(const unsigned &globalSize) {
        _sVec.resize(globalSize);
        std::fill(_sVec.begin(), _sVec.end(), 0);
        _vecIsSerial = true;
        _vecIsParallel = false;
        _vecIsLocalized = false;
      }
      
      MyVector(const std::vector < unsigned > & vecOffset, const unsigned &iproc) {

        _vecOffset = vecOffset;
        _iproc = iproc;
        _nprocs = _vecOffset.size();

        _pVec.resize(_vecOffset[_iproc + 1] - vecOffset[_iproc]);
	std::fill(_pVec.begin(), _pVec.end(), 0);

        _vecIsSerial = false;
        _vecIsParallel = true;
        _vecIsLocalized = false;
      }

      // ******************
      ~MyVector() {
        std::vector<Type>().swap(_sVec);
        _vecIsSerial = false;
        std::vector<Type>().swap(_pVec);
        _vecIsParallel = false;
        std::vector<Type>().swap(_lVec);
        _vecIsLocalized = false;
      }

      // ******************
      void ScatterVector(const std::vector < unsigned > & vecOffset, const unsigned &iproc, const unsigned &nprocs) {
 
        _vecOffset = vecOffset;
        _iproc = iproc;
        _nprocs = nprocs;

        _pVec.resize(_vecOffset[_iproc + 1] - vecOffset[_iproc]);

        for (unsigned i = _vecOffset[_iproc]; i < _vecOffset[_iproc + 1]; i++) {
          _pVec[i - _vecOffset[_iproc] ] = _sVec[i];
        }

        std::vector<Type>().swap(_sVec);

        _vecIsSerial = false;
        _vecIsParallel = true;
        _vecIsLocalized = false;
      }

      // ******************
      void LocalizeVectorFromOneToAll(const unsigned &lproc) {

        if (!_vecIsParallel) {
          std::cout << "Error in MyVector.LocalizeVectorFromOneToAll(), vector is in a wrong status" << std::endl;
          abort();
        }

        _lproc = lproc;
        _lVec.resize(_vecOffset[_lproc + 1] - _vecOffset[_lproc]);

        if (_iproc == _lproc) {
          for (unsigned i = 0; i < _lVec.resize(); i++) {
            _lVec[i] = _pVec[i];
          }
        }

        MPI_Bcast(&_lVec[0], _lVec.size(), sizeof(Type), _lproc, MPI_COMM_WORLD);
        _vecIsLocalized = true;
      }

      // ******************
      void FreeLocalizedVector() {
        std::vector<Type>().swap(_lVec);
        _vecIsLocalized = false;
      }

      // ******************
      void Set(const unsigned &i, const Type &value) {
        if (_vecIsSerial) {
          _sVec[i] = value;
        }
        else if (_vecIsParallel && !_vecIsLocalized) {
          _pVec[i - _vecOffset[_iproc] ] = value;
        }
        else {
          std::cout << "Error in MyVector.Set(), vector is in a wrong status" << std::endl;
          abort();
        }
      }

      // ******************
      Type Get(const unsigned &i) {
        if (_vecIsSerial) {
          return _sVec[i];
        }
        else if (_vecIsParallel && !_vecIsLocalized) {
          return _pVec[i - _vecOffset[_iproc] ];
        }
        else if (_vecIsLocalized) {
          return _pVec[i - _vecOffset[_lproc] ];
        }
        else {
          std::cout << "Error in MyVector.Get(), vector is in a wrong status" << std::endl;
          abort();
        }
      }

    private:
      std::vector< Type > _sVec;
      bool _vecIsSerial;

      std::vector< Type > _pVec;
      bool _vecIsParallel;

      std::vector < unsigned > _vecOffset;
      unsigned _iproc;
      unsigned _nprocs;

      std::vector< Type > _lVec;
      bool _vecIsLocalized;
      unsigned _lproc;

  };


} //end namespace femus

#endif
