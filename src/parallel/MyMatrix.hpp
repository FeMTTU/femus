
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

#ifndef __femus_mesh_MyMatrix_hpp__
#define __femus_mesh_MyMatrix_hpp__


#include <iostream>
#include <vector>
#include <stdlib.h>

#include <mpi.h>
#include <boost/mpi/datatype.hpp>

#include "MyVector.hpp"

namespace femus {

  template <class Type> class MyMatrix {

    public:
      // ******************
      MyMatrix();

      // ******************
      MyMatrix(const unsigned &rsize, const unsigned &csize, const Type value = 0);

      // ******************
      MyMatrix(const std::vector < unsigned > &offset, const unsigned &csize, const Type value = 0);

      // ******************
      ~MyMatrix();

      // ******************
      void resize(const unsigned &rsize, const unsigned &csize, const Type value = 0);

      // ******************
      void resize(const std::vector < unsigned > &offset, const unsigned &csize, const Type value = 0);

      // ******************
      void clear();

      // ******************
      unsigned size();

      // ******************
      unsigned begin();

      // ******************
      unsigned end();

      // ******************
      void scatter(const std::vector < unsigned > &offset);

      // ******************
      void scatter();

      // ******************
      void buildOffset();

      // ******************
      std::vector<unsigned> getOffset();

      // ******************
      void localizeToAll(const unsigned &lproc);

      // ******************
      void localizeToOne(const unsigned &lproc, const unsigned &kproc);

      // ******************
      void clearLocalized();

      // ****************
      const std::string &status();

      // ******************
      Type& operator()(const unsigned &i, const unsigned &j);

      MyMatrix<Type>& operator()(const unsigned &i);

      // *****************
      friend std::ostream& operator<<(std::ostream& os, MyMatrix<Type>& mat) {

        os << mat.status() << std::endl;

        if(mat._matIsAllocated) {
          if(mat._serial) {
            for(unsigned i = mat.begin(); i < mat.end(); i++) {
	      for(unsigned j = mat(i).begin(); j < mat(i).end(); j++) {
                os << mat(i, j) << " ";
              }
              os << std::endl;
            }
          }
          else {
//             for(int j = 0; j < vec._nprocs; j++) {
//               vec.localizeToAll(j);
//               for(unsigned i = vec.begin(); i < vec.end(); i++) {
//                 os << i << " " << vec[i] << std::endl;
//               }
//               vec.clearLocalized();
//             }
          }
        }
        return os;
      }

    private:

      std::string _status;
      bool _serial;
      bool _matIsAllocated;

      unsigned _iproc;
      unsigned _nprocs;
      MPI_Datatype _MY_MPI_DATATYPE;
      Type _dummy;

      unsigned _begin;
      unsigned _end;
      unsigned _size;

      MyVector < unsigned > _rowOffset;
      MyVector < unsigned > _rowSize;

      std::vector< Type > _mat;
      std::vector< Type > _mat2;
      std::vector < unsigned > _offset;

      unsigned _lproc;
      bool _irowIsSet;
      unsigned _irow;
  };





} //end namespace femus



#endif
