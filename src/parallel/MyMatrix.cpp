
/*=========================================================================

 Program: FEMuS
 Module: MyMatrix
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

#include "MyMatrix.hpp"

namespace femus {

  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix() {
    init();
  }

  template <class Type> void MyMatrix<Type>::init() {
    int iproc, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    _iproc = static_cast < unsigned >(iproc);
    _nprocs = static_cast < unsigned >(nprocs);

    Type dummy = 0;
    _MY_MPI_DATATYPE = boost::mpi::get_mpi_datatype(dummy);

    _matIsAllocated = false;
    _serial = true;

    _begin = 0;
    _end = 0;
    _size = 0;
  }
  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix(const unsigned &rsize, const unsigned &csize, const Type value) {
    init();
    resize(rsize, csize, value);
  }
  
  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix(const std::vector < unsigned > &offset, const unsigned &csize, const Type value) {
    init();
    resize(offset, csize, value);
  }

  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix(const MyVector < unsigned > &rowSize, const Type value) {
    init();
    
    _rowSize = rowSize;
    _begin = _rowSize.begin();
    _end = _rowSize.end();
    _size = _rowSize.size();
    _matSize.resize(_nprocs);
    
    if( (_rowSize.status()).compare("UNINITIALIZED") == 0 ){
      std::cout << "In function MyMatrix(const MyVector < unsigned > &rowSize, const Type value), rowSize is UNINITIALIZED"<<std::endl;
      abort();
    }
    else  if( (_rowSize.status()).compare("PARALLEL") == 0 ){
      _offset = _rowSize.getOffset();
      _rowOffset.resize(_offset);
      _matSize.scatter();
      _serial = false;
    }
    else{
      _rowOffset.resize(_size);
    }
       
    _rowOffset[_rowOffset.begin()] = 0;
    unsigned matsize = _rowSize[_rowOffset.begin()];
    for(unsigned i = _rowOffset.begin() + 1; i < _rowOffset.end(); i++) {
      _rowOffset[i] = _rowOffset[i - 1] + _rowSize[i - 1] ;
      matsize += _rowSize[i - 1];
    }
    _matSize[_iproc] = matsize;
    _mat.resize(_matSize[_iproc], value);

    _matIsAllocated = true;
  }
  
  // ******************
  template <class Type> MyMatrix<Type>::~MyMatrix() {
    clear();
  }

  // ******************
  template <class Type> void MyMatrix<Type>::resize(const unsigned &rsize, const unsigned &csize, const Type value) {

    _begin = 0;
    _end = rsize;
    _size = _end - _begin;

    _rowSize.resize(_size, csize);
    _rowOffset.resize(_size);
    _rowOffset[_rowOffset.begin()] = 0;
    for(unsigned i = _rowOffset.begin() + 1; i < _rowOffset.end(); i++) {
      _rowOffset[i] = _rowOffset[i - 1] + _rowSize[i - 1] ;
    }

    _matSize.resize(_nprocs);
    _matSize[_iproc] = _size * csize;
    _mat.resize(_matSize[_iproc], value);

    _matIsAllocated = true;
    _serial = true;

  }

  // ******************
  template <class Type> void MyMatrix<Type>::resize(const std::vector < unsigned > &offset, const unsigned &csize, const Type value) {

    _offset = offset;
    if(_nprocs != _offset.size() - 1) {
      std::cout << "Error in MyMatrix.Resize(...), offset.size()" << _offset.size()
                << "!= from nprocs+1=" << _nprocs + 1 << std::endl;
      abort();
    }

    _begin = _offset[_iproc];
    _end = _offset[_iproc + 1];
    _size = _end - _begin;

    _rowSize.resize(offset, csize);
    _rowOffset.resize(offset);
    _rowOffset[_rowOffset.begin()] = 0;
    for(unsigned i = _rowOffset.begin() + 1; i < _rowOffset.end(); i++) {
      _rowOffset[i] = _rowOffset[i - 1] + _rowSize[i - 1] ;
    }

    _matSize.resize(_nprocs);
    _matSize.scatter();
    _matSize[_iproc] = _size * csize;
    _mat.resize(_matSize[_iproc], value);

    _matIsAllocated = true;
    _serial = false;

  }

  // ******************
  template <class Type> void MyMatrix<Type>::clear() {
    std::vector<Type>().swap(_mat);
    std::vector<Type>().swap(_mat2);
    _rowOffset.clear();
    _matIsAllocated = false;
    _serial = true;
  }

  template <class Type> Type* MyMatrix<Type>::operator[](const unsigned &i) {
    return &_mat[ _rowOffset[i]];
  }

  template <class Type> Type& MyMatrix<Type>::operator()(const unsigned &i, const unsigned &j) {
    return _mat[ _rowOffset[i] + j];
  }

  // ******************
  template <class Type> unsigned MyMatrix<Type>::size() {
    return _size;
  }

  template <class Type> unsigned MyMatrix<Type>::size(const unsigned &i) {
    return (_matIsAllocated) ? _rowSize[i] : 0;
  }

  // ******************
  template <class Type> unsigned MyMatrix<Type>::begin() {
    return _begin;
  }

  template <class Type> unsigned MyMatrix<Type>::begin(const unsigned &i) {
    return 0;
  }

  // ******************
  template <class Type> unsigned MyMatrix<Type>::end() {
    return _end;
  }

  template <class Type> unsigned MyMatrix<Type>::end(const unsigned &i) {
    return (_matIsAllocated) ? _rowSize[i] : 0;
  }

  // ******************
  template <class Type> void MyMatrix<Type>::scatter(const std::vector < unsigned > &offset) {

    _offset = offset;

    if(!_serial) {
      std::cout << "Error in MyMatrix.scatter(), matrix is in " << status() << " status" << std::endl;
      abort();
    }

    if(_nprocs != _offset.size() - 1) {
      std::cout << "Error in MyMatrix.resize(...), offset.size() != from nprocs" << std::endl;
      abort();
    }

    _begin = _offset[_iproc];
    _end = _offset[_iproc + 1];
    _size = _end - _begin;

    MyVector < unsigned > rowOffset2 = _rowOffset;

    _rowSize.scatter(_offset);
    _rowOffset.scatter(_offset);
    _rowOffset[_rowOffset.begin()] = 0;
    unsigned matsize = _rowSize[_rowOffset.begin()];
    for(unsigned i = _rowOffset.begin() + 1; i < _rowOffset.end(); i++) {
      _rowOffset[i] = _rowOffset[i - 1] + _rowSize[i - 1] ;
      matsize += _rowSize[i - 1];
    }

    _matSize.scatter();
    _matSize[_iproc] = matsize;

    _mat.swap(_mat2);

    _mat.resize(_matSize[_iproc]);

    for(unsigned i = _offset[_iproc]; i < _offset[_iproc + 1]; i++) {
      for(unsigned j = begin(i); j < end(i); j++) {
        _mat[_rowOffset[i] + j] = _mat2[rowOffset2[i] + j];
      }
    }

    std::vector<Type>().swap(_mat2);

    _serial = false;

  }

  // ******************
  template <class Type> void MyMatrix<Type>::scatter() {
    buildOffset();
    scatter(_offset);
  }

  // ******************
  template <class Type> void MyMatrix<Type>::buildOffset() {
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
  template <class Type> std::vector<unsigned>  MyMatrix<Type>::getOffset() {
    return _offset;
  }

  // ******************
  template <class Type> void MyMatrix<Type>::localize(const unsigned &lproc) {

    if(_serial) {
      std::cout << "Error in MyMatrix.LocalizeToAll(), matrix is in " << status() << " status" << std::endl;
      abort();
    }

    _matSize.localize(lproc);
    _rowSize.localize(lproc);
    _rowOffset.localize(lproc);

    if(_iproc != lproc) {
      _mat.swap(_mat2);
      _mat.resize(_matSize[lproc]);
    }

    MPI_Bcast(&_mat[0], _matSize[lproc], _MY_MPI_DATATYPE, lproc, MPI_COMM_WORLD);

    _begin = _offset[lproc];
    _end = _offset[lproc + 1];
    _size = _end - _begin;

    _lproc = lproc;
  }

  // ******************
  template <class Type> void MyMatrix<Type>::clearLocalized() {

    _matSize.clearLocalized();
    _rowSize.clearLocalized();
    _rowOffset.clearLocalized();

    if(_lproc != _iproc) {
      _mat.swap(_mat2);
      std::vector<Type>().swap(_mat2);
      _begin = _offset[_iproc];
      _end = _offset[_iproc + 1];
      _size = _end - _begin;
    }
  }

  // ****************
  template <class Type> const std::string & MyMatrix<Type>::status() {

    if(!_matIsAllocated)
      _status = "UNINITIALIZED";
    else if(_serial)
      _status = "SERIAL";
    else
      _status = "PARALLEL";

    return _status;
  }

  // ******************

  // Explicit template instantiation
  template class MyMatrix<float>;
  template class MyMatrix<double>;
  template class MyMatrix<long double>;
  template class MyMatrix<int>;
  template class MyMatrix<short int>;
  template class MyMatrix<long int>;
  template class MyMatrix<short unsigned int>;
  template class MyMatrix<unsigned int>;
  template class MyMatrix<long unsigned int>;
  template class MyMatrix<char>;

} //end namespace femus






