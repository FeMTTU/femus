
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

    int iproc, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    _iproc = static_cast < unsigned >(iproc);
    _nprocs = static_cast < unsigned >(nprocs);

    _dummy = 0;
    _MY_MPI_DATATYPE = boost::mpi::get_mpi_datatype(_dummy);

    _matIsAllocated = false;
    _serial = true;

    _begin = 0;
    _end = 0;
    _size = 0;

  }

  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix(const unsigned &rsize, const unsigned &csize, const Type value) {
    MyMatrix();
    resize(rsize, csize, value);
  }

  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix(const std::vector < unsigned > &offset, const unsigned &csize, const Type value) {
    MyMatrix();
    resize(offset, csize, value);
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
      std::cout << "Error in MyMatrix.Resize(...), offset.size() != from nprocs" << std::endl;
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


  //This function should never be used (as mat[i][j])
  //when the matrix is in localizedToOne status.
  //Use mat(i,j) instead
  template <class Type> Type* MyMatrix<Type>::operator[](const unsigned &i) {
    return &_mat[ _rowOffset[i]];
  }


  template <class Type> Type& MyMatrix<Type>::operator()(const unsigned &i, const unsigned &j) {
    return (_matIsAllocated) ? _mat[ _rowOffset[i] + j] : _dummy;
  }


  // ******************

  template <class Type> unsigned MyMatrix<Type>::size() {
    return _size;
  }

  template <class Type> unsigned MyMatrix<Type>::size(const unsigned & i) {
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

    //std::cout<< _offset[_iproc] <<" "<< _offset[_iproc + 1] << " "<<_matSize[_iproc] << std::endl;
    
    for(unsigned i = _offset[_iproc]; i < _offset[_iproc + 1]; i++) {
      //std::cout<< _rowOffset[i] << " " << end(i) << std::endl;
      for(unsigned j = begin(i); j< end(i); j++){
	
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
  template <class Type> void MyMatrix<Type>::localizeToAll(const unsigned &lproc) {

    if(_serial) {
      std::cout << "Error in MyMatrix.LocalizeToAll(), matrix is in " << status() << " status" << std::endl;
      abort();
    }

    _matSize.localizeToAll(lproc);
    _rowSize.localizeToAll(lproc);
    _rowOffset.localizeToAll(lproc);

    
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
//
//   // ******************
//   template <class Type> void MyMatrix<Type>::localizeToOne(const unsigned &lproc, const unsigned &kproc) {
//
//     if(_serial) {
//       std::cout << "Error in MyMatrix.LocalizeToAll(), vector is in " << status() << " status" << std::endl;
//       abort();
//     }
//
//     if(_iproc == lproc) {
//       if(_iproc != kproc){
// 	MPI_Send(&_vec[0], _vec.size(), _MY_MPI_DATATYPE, kproc, 1, MPI_COMM_WORLD);
// 	_matIsAllocated = false;
//       }
//     }
//     else if(_iproc == kproc) {
//       _vec.swap(_vec2);
//       _vec.resize(_offset[lproc + 1] - _offset[lproc]);
//       MPI_Recv(&_vec[0], _vec.size(), _MY_MPI_DATATYPE, lproc, 1, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
//       _begin = _offset[lproc];
//       _end = _offset[lproc + 1];
//     }
//     else {
//       _matIsAllocated = false;
//     }
//
//     _lproc = lproc;
//   }
//
//   // ******************
  template <class Type> void MyMatrix<Type>::clearLocalized() {

    _matSize.clearLocalized();
    _rowOffset.clearLocalized();
    _rowSize.clearLocalized();
    
    if(_lproc != _iproc && _matIsAllocated) {
      _mat.swap(_mat2);
      std::vector<Type>().swap(_mat2);
      _begin = _offset[_iproc];
      _end = _offset[_iproc + 1];
      _size = _end - _begin;
    }
    else {
      _matIsAllocated = true;
    }

  }
//
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






