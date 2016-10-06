
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

#include <string>

#include <iostream>
#include <vector>
#include <stdlib.h>

#include <mpi.h>
#include <boost/mpi/datatype.hpp>

#include "MyVector.hpp"

namespace femus {

  // ******************
  template <class Type> MyVector<Type>::MyVector() {

    int iproc, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    _iproc = static_cast < unsigned >(iproc);
    _nprocs = static_cast < unsigned >(nprocs);

    _dummy = 0;
    _MY_MPI_DATATYPE = boost::mpi::get_mpi_datatype(_dummy);

    _vecIsAllocated = false;
    _serial = true;

    _begin = 0;
    _end = 0;
    _size = 0;
  }

  // ******************
  template <class Type> MyVector<Type>::MyVector(const unsigned &size, const Type value) {
    MyVector();
    resize(size, value);
  }

  // ******************
  template <class Type> MyVector<Type>::MyVector(const std::vector < unsigned > &offset, const Type value) {
    MyVector();
    resize(offset, value);
  }

  // ******************
  template <class Type> MyVector<Type>::~MyVector() {
    clear();
  }

  // ******************
  template <class Type> void MyVector<Type>::resize(const unsigned &size, const Type value) {

    _vec.resize(size, value);

    _vecIsAllocated = true;
    _serial = true;

    _begin = 0;
    _end = size;
    _size = size;

  }

  // ******************
  template <class Type> void MyVector<Type>::resize(const std::vector < unsigned > &offset, const Type value) {

    _offset = offset;
    if(_nprocs != _offset.size() - 1) {
      std::cout << "Error in MyVector.Resize(...), offset.size() != from nprocs" << std::endl;
      abort();
    }
       
    _vec.resize(_offset[_iproc + 1] - offset[_iproc], value);
   
    _vecIsAllocated = true;
    _serial = false;

    _begin = _offset[_iproc];
    _end = _offset[_iproc + 1];
    _size = _offset[_nprocs];

  }

  // ******************
  template <class Type> void MyVector<Type>::clear() {
    std::vector<Type>().swap(_vec);
    std::vector<Type>().swap(_vec2);
    _vecIsAllocated = false;
    _serial = true;
  }


  // ******************
  template <class Type> unsigned MyVector<Type>::size() {
    return _size;
  }

  // ******************
  template <class Type> unsigned MyVector<Type>::begin() {
    return _begin;
  }

  // ******************
  template <class Type> unsigned MyVector<Type>::end() {
    return _end;
  }

  // ******************
  template <class Type> void MyVector<Type>::scatter(const std::vector < unsigned > &offset) {

    _offset = offset;

    if(!_serial) {
      std::cout << "Error in MyVector.scatter(), vector is in " << status() << " status" << std::endl;
      abort();
    }

    if(_nprocs != _offset.size() - 1) {
      std::cout << "Error in MyVector.resize(...), offset.size() != from nprocs" << std::endl;
      abort();
    }

    _vec.swap(_vec2);

    _vec.resize(_offset[_iproc + 1] - _offset[_iproc]);

    for(unsigned i = _offset[_iproc]; i < _offset[_iproc + 1]; i++) {
      _vec[i - _offset[_iproc] ] = _vec2[i];
    }

    std::vector<Type>().swap(_vec2);

    _serial = false;

    _begin = _offset[_iproc];
    _end = _offset[_iproc + 1];
    _size = _offset[_nprocs];

  }

  // ******************
  template <class Type> void MyVector<Type>::scatter() {
    buildOffset();
    scatter(_offset);
  }

  // ******************
  template <class Type> void MyVector<Type>::buildOffset() {
    unsigned locsize = _size / _nprocs;
    unsigned reminder = _size % _nprocs;
    _offset.resize(_nprocs + 1);
    _offset[0] = 0;
    for(unsigned i = 1; i < _nprocs; i++) {
      _offset[i] = _offset[i - 1] + locsize;
      if(i <= reminder) _offset[i] += 1;
    }
    _offset[_nprocs] = _size;
  }

  // ******************
  template <class Type> std::vector<unsigned>  MyVector<Type>::getOffset() {
    return _offset;
  }

  // ******************
  template <class Type> void MyVector<Type>::localizeToAll(const unsigned &lproc) {

    if(_serial) {
      std::cout << "Error in MyVector.LocalizeToAll(), vector is in " << status() << " status" << std::endl;
      abort();
    }

    if(_iproc != lproc) {
      _vec.swap(_vec2);
      _vec.resize(_offset[lproc + 1] - _offset[lproc]);
    }

    MPI_Bcast(&_vec[0], _vec.size(), _MY_MPI_DATATYPE, lproc, MPI_COMM_WORLD);

    _begin = _offset[lproc];
    _end = _offset[lproc + 1];

    _lproc = lproc;

  }

  // ******************
  template <class Type> void MyVector<Type>::localizeToOne(const unsigned &lproc, const unsigned &kproc) {

    if(_serial) {
      std::cout << "Error in MyVector.LocalizeToAll(), vector is in " << status() << " status" << std::endl;
      abort();
    }

    if(_iproc == lproc) {
      if(_iproc != kproc){
	MPI_Send(&_vec[0], _vec.size(), _MY_MPI_DATATYPE, kproc, 1, MPI_COMM_WORLD);
	_vecIsAllocated = false;
      }
    }
    else if(_iproc == kproc) {
      _vec.swap(_vec2);
      _vec.resize(_offset[lproc + 1] - _offset[lproc]);
      MPI_Recv(&_vec[0], _vec.size(), _MY_MPI_DATATYPE, lproc, 1, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
      _begin = _offset[lproc];
      _end = _offset[lproc + 1];
    }
    else {
      _vecIsAllocated = false;
    }

    _lproc = lproc;
  }

  // ******************
  template <class Type> void MyVector<Type>::clearLocalized() {

    if(_lproc != _iproc && _vecIsAllocated) {
      _vec.swap(_vec2);
      std::vector<Type>().swap(_vec2);
      _begin = _offset[_iproc];
      _end = _offset[_iproc + 1];
      _size = _offset[_nprocs];
    }
    else {
      _vecIsAllocated = true;
    }

  }

  // ****************
  template <class Type> const std::string & MyVector<Type>::status() {

    if(!_vecIsAllocated)
      _status = "UNINITIALIZED";
    else if(_serial)
      _status = "SERIAL";
    else
      _status = "PARALLEL";

    return _status;
  }

  // ******************
  template <class Type> Type& MyVector<Type>::operator[](const unsigned &i) {
    return (_vecIsAllocated) ? _vec[i - _begin] : _dummy;
  }

  // Explicit template instantiation
  template class MyVector<float>;
  template class MyVector<double>;
  template class MyVector<long double>;
  template class MyVector<int>;
  template class MyVector<short int>;
  template class MyVector<long int>;
  template class MyVector<short unsigned int>;
  template class MyVector<unsigned int>;
  template class MyVector<long unsigned int>;
  template class MyVector<char>;

} //end namespace femus






