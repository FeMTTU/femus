
/*=========================================================================

 Program: FEMuS
 Module: MyVector
 Authors: Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/


#ifndef __femus_mesh_MyVector_hpp__
#define __femus_mesh_MyVector_hpp__


#include <iostream>
#include <vector>
#include <stdlib.h>

#include <mpi.h>
#include <boost/mpi/datatype.hpp>



namespace femus {

  template <class Type> class MyVector {

    public:
      // ******************
      MyVector() {
        _vecIsSerial = false;
        _vecIsParallel = false;
        _vecIsLocalized = false;

        int iproc, nprocs;

        MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

        _iproc = static_cast < unsigned >(iproc);
        _nprocs = static_cast < unsigned >(nprocs);

        Type a;
        _MPI_MYDATATYPE = boost::mpi::get_mpi_datatype(a);
      }


      // ******************
      MyVector(const unsigned &size, const Type value = 0) {
        MyVector();
        resize(size, value);
      }

      // ******************
      MyVector(const std::vector < unsigned > &offset, const Type value = 0) {
        MyVector();
        resize(offset, value);
      }

      // ******************
      void resize(const unsigned &globalSize, const Type value = 0) {
        clearParallel();
        clearLocalized();
        _sVec.resize(globalSize, value);
        _vecIsSerial = true;
      }

      // ******************
      void resize(const std::vector < unsigned > &offset, const Type value = 0) {
        if(_nprocs != _offset.size() - 1) {
          std::cout << "Error in MyVector.Resize(...), offset.size() != from nprocs" << std::endl;
          abort();
        }
        clearSerial();
        clearLocalized();
        _offset = offset;
        _pVec.resize(_offset[_iproc + 1] - offset[_iproc], value);
        _vecIsParallel = true;
      }

      // ******************
      ~MyVector() {
        clear();
      }

      // ******************
      void clearSerial() {
        std::vector<Type>().swap(_sVec);
        _vecIsSerial = false;
      }

      // ******************
      void clearParallel() {
        std::vector<Type>().swap(_pVec);
        _vecIsParallel = false;
      }

      // ******************
      void clearLocalized() {
        std::vector<Type>().swap(_lVec);
        _vecIsLocalized = false;
      }

      // ******************
      void clear() {
        clearSerial();
        clearParallel();
        clearLocalized();
      }


      // ******************
      unsigned size() {
        if(_vecIsSerial)
          return _sVec.size();
        else if(_vecIsParallel)
          return _offset[_nprocs];
        else {
          std::cout << "Error in MyVector.size(), vector is in a wrong status" << std::endl;
          abort();
        }
      }

      // ******************
      unsigned begin() {
        if(_vecIsSerial)
          return 0;
        else if(_vecIsParallel && !_vecIsLocalized)
          return _offset[_iproc];
        else if(_vecIsLocalized) {
          return _offset[_lproc];
        }
        else {
          std::cout << "Error in MyVector.begin(), vector is in a wrong status" << std::endl;
          abort();
        }
      }

      // ******************
      unsigned end() {
        if(_vecIsSerial)
          return _sVec.size();
        else if(_vecIsParallel && !_vecIsLocalized)
          return _offset[_iproc + 1];
        else if(_vecIsLocalized) {
          return _offset[_lproc + 1];
        }
        else {
          std::cout << "Error in MyVector.end(), vector is in a wrong status" << std::endl;
          abort();
        }

      }

      // ******************
      void scatter(const std::vector < unsigned > & offset) {

        _offset = offset;

        if(!_vecIsSerial) {
          std::cout << "Error in MyVector.scatter(), vector is in a wrong status" << std::endl;
          abort();
        }

        if(_nprocs != _offset.size() - 1) {
          std::cout << "Error in MyVector.resize(...), offset.size() != from nprocs" << std::endl;
          abort();
        }

        _pVec.resize(_offset[_iproc + 1] - offset[_iproc]);

        for(unsigned i = _offset[_iproc]; i < _offset[_iproc + 1]; i++) {
          _pVec[i - _offset[_iproc] ] = _sVec[i];
        }

        clearSerial();
        clearLocalized();
        _vecIsParallel = true;
      }

      // ******************
      void localizeToAll(const unsigned &lproc) {

        if(!_vecIsParallel) {
          std::cout << "Error in MyVector.LocalizeToAll(), vector is in a wrong status" << std::endl;
          abort();
        }

        _lproc = lproc;
        _lVec.resize(_offset[_lproc + 1] - _offset[_lproc]);

        if(_iproc == _lproc) {
          for(unsigned i = 0; i < _lVec.size(); i++) {
            _lVec[i] = _pVec[i];
          }
        }

        MPI_Bcast(&_lVec[0], _lVec.size(), _MPI_MYDATATYPE, _lproc, MPI_COMM_WORLD);
        _vecIsLocalized = true;
      }

      // ******************
      void set(const unsigned &i, const Type &value) {
        if(_vecIsSerial)
          _sVec[i] = value;
        else if(_vecIsParallel && !_vecIsLocalized)
          _pVec[i - _offset[_iproc] ] = value;
        else {
          std::cout << "Error in MyVector.set(), vector is in a wrong status" << std::endl;
          abort();
        }
      }

      // ******************
      Type get(const unsigned &i) {
        if(_vecIsSerial)
          return _sVec[i];
        else if(_vecIsParallel && !_vecIsLocalized)
          return _pVec[i - _offset[_iproc] ];
        else if(_vecIsLocalized)
          return _lVec[i - _offset[_lproc] ];
        else {
          std::cout << "Error in MyVector.get(), vector is in a wrong status" << std::endl;
          abort();
        }
      }

    private:
      std::vector< Type > _sVec;
      bool _vecIsSerial;

      std::vector< Type > _pVec;
      bool _vecIsParallel;

      std::vector < unsigned > _offset;
      unsigned _iproc;
      unsigned _nprocs;
      std::vector< Type > _lVec;
      bool _vecIsLocalized;
      unsigned _lproc;

      MPI_Datatype _MPI_MYDATATYPE;

  };


} //end namespace femus

#endif
