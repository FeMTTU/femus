
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
        _vecIsallocated = false;
        _serial = false;
	_parallel = false;
	
        int iproc, nprocs;

        MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

        _iproc = static_cast < unsigned >(iproc);
        _nprocs = static_cast < unsigned >(nprocs);

        _dummy = 0;
        _MPI_MYDATATYPE = boost::mpi::get_mpi_datatype(_dummy);
	
	_begin = 0;
	_end = 0;
	_size = 0;
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
        _vec.resize(globalSize, value);
        
	_begin = 0;
	_end = globalSize;
	_size = globalSize;
      }

      // ******************
      void resize(const std::vector < unsigned > &offset, const Type value = 0) {
        if(_nprocs != _offset.size() - 1) {
          std::cout << "Error in MyVector.Resize(...), offset.size() != from nprocs" << std::endl;
          abort();
        }
        clearSerial();
        _offset = offset;
        _vec.resize(_offset[_iproc + 1] - offset[_iproc], value);
        
	_vecIsallocated = true;
	_begin = _offset[_iproc];
	_end = _offset[_iproc+1];
	_size = _offset[_nprocs];
      }

      // ******************
      ~MyVector() {
        clear();
      }

      // ******************
      void clearSerial() {
        std::vector<Type>().swap(_sVec);
      }

      // ******************
      void clearParallel() {
        std::vector<Type>().swap(_pVec);
      }

      // ******************
      void clearLocalized() {
        std::vector<Type>().swap(_lVec);
      	
	_vec.swap(_pVec);
	std::vector<Type>().swap(_pVec);
	_begin = _offset[_iproc];
	_end = _offset[_iproc+1];
	_size = _offset[_nprocs];
      }

      // ******************
      void clear() {
        clearSerial();
        clearParallel();
        std::vector<Type>().swap(_vec);
      }


      // ******************
      unsigned size() {
	return _size;
      }

      // ******************
      unsigned begin() {
	return _begin;
      }

      // ******************
      unsigned end() {
	return _end;
      }

      // ******************
      void scatter(const std::vector < unsigned > &offset) {

        _offset = offset;

        if(!_vecIsallocated) {
          std::cout << "Error in MyVector.scatter(), vector is in a wrong status" << std::endl;
          abort();
        }

        if(_nprocs != _offset.size() - 1) {
          std::cout << "Error in MyVector.resize(...), offset.size() != from nprocs" << std::endl;
          abort();
        }

        _pVec.resize(_offset[_iproc + 1] - offset[_iproc]);

        for(unsigned i = _offset[_iproc]; i < _offset[_iproc + 1]; i++) {
          _pVec[i - _offset[_iproc] ] = _vec[i];
        }

        clearSerial();
	
        _vecIsallocated = true;
	_vec.swap(_pVec);
	std::vector<Type>().swap(_pVec);
	_begin = _offset[_iproc];
	_end = _offset[_iproc+1];
	_size = _offset[_nprocs];
      }

      // ******************
      void localizeToAll(const unsigned &lproc) {
	
        if(_serial) {
          std::cout << "Error in MyVector.LocalizeToAll(), vector is in a wrong status" << std::endl;
          abort();
        }

        _vec.swap(_pVec);
        
        _lproc = lproc;
        _lVec.resize(_offset[_lproc + 1] - _offset[_lproc]);

        if(_iproc == _lproc) {
          for(unsigned i = 0; i < _lVec.size(); i++) {
            _lVec[i] = _pVec[i];
          }
        }

        MPI_Bcast(&_lVec[0], _lVec.size(), _MPI_MYDATATYPE, _lproc, MPI_COMM_WORLD);
        
	_vecIsallocated = true;
	
	_vec.swap(_lVec);
	std::vector<Type>().swap(_lVec);
	_begin = _offset[_lproc];
	_end = _offset[_lproc+1];
	
      }

      // ******************
      void localizeToOne(const unsigned &lproc, const unsigned &kproc) {

        if(!_serial) {
          std::cout << "Error in MyVector.LocalizeToAll(), vector is in a wrong status" << std::endl;
          abort();
        }

        _vec.swap(_pVec);
        
        if(_iproc == kproc)
          _vecIsallocated = true;
        else
          _vecIsallocated = false;

        if(_iproc == lproc || _iproc == kproc) {
          _lproc = lproc;
          _lVec.resize(_offset[_lproc + 1] - _offset[_lproc]);

          if(_iproc == _lproc) {
            for(unsigned i = 0; i < _lVec.size(); i++) {
              _lVec[i] = _pVec[i];
            }
          }

          if(lproc != kproc) {
            if(_iproc == lproc)
              MPI_Send(&_lVec[0], _lVec.size(), _MPI_MYDATATYPE, kproc, 1, MPI_COMM_WORLD);
            else
              MPI_Recv(&_lVec[0], _lVec.size(), _MPI_MYDATATYPE, lproc, 1, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
          }
        }
        
	_vec.swap(_lVec);
	std::vector<Type>().swap(_lVec);
	_begin = _offset[_lproc];
	_end = _offset[_lproc+1];
      }

      // ******************
      Type& operator[](const unsigned &i) {
	if(_vecIsallocated){
	  return _vec[i-_begin];
	}
	else{  
	  return _dummy;
	}
	
      }

      // *****************
      friend std::ostream& operator<<(std::ostream& os, MyVector<Type>& vec) {
	
	if(vec._serial) {	  
          for(unsigned i = vec.begin(); i < vec.end(); i++) {
            os << i << " " << vec[i] << std::endl;
          }
        }
        else {
          for(int j = 0; j < vec._nprocs; j++) {
            vec.localizeToOne(j, 0);
            for(unsigned i = vec.begin(); i < vec.end(); i++) {
              os << i << " " << vec[i] << std::endl;
            }
            vec.clearLocalized();
          }
        }
        return os;
      }

    private:
      
      bool _serial;
      bool _parallel;
      
      unsigned _iproc;
      unsigned _nprocs;
      MPI_Datatype _MPI_MYDATATYPE;
      Type _dummy;
      
      unsigned _begin;
      unsigned _end;
      unsigned _size;
      
      std::vector< Type > _vec;
      bool _vecIsallocated;
      
      std::vector< Type > _sVec;
     
      std::vector< Type > _pVec;
      std::vector < unsigned > _offset;

      std::vector< Type > _lVec;
      unsigned _lproc;


  };


} //end namespace femus



#endif
