
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


#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>

#include <mpi.h>



namespace femus {

  template <class Type> class MyVector {

    public:
      // ******************
      MyVector();

      // ******************
      MyVector(const unsigned &size, const Type value = 0);

      // ******************
      MyVector(const std::vector < unsigned > &offset, const Type value = 0);

      // ******************
      ~MyVector();

      // ******************
      void init();

      //*******************
      void resize(const unsigned &size, const Type value = 0);

      // ******************
      void resize(const std::vector < unsigned > &offset, const Type value = 0);

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
      void stack();

      // ******************
      void buildOffset();

      // ******************
      std::vector<unsigned> getOffset();

      // ******************
      void broadcast(const unsigned &lproc);

      // ******************
      void clearBroadcast();

      // ****************
      const std::string &status();

      // ******************
      Type& operator[](const unsigned &i);

      // ******************
      const Type& operator[](const unsigned &i) const;

      // *****************
      friend std::ostream& operator<<(std::ostream& os, MyVector<Type>& vec) {

        os << vec.status() << std::endl;

        if(vec._vecIsAllocated) {
          if(vec._serial) {
            for(unsigned i = vec.begin(); i < vec.end(); i++) {
              os << i << " " << vec[i] << std::endl;
            }
          }
          else {
            for(int j = 0; j < vec._nprocs; j++) {
              vec.broadcast(j);
              for(unsigned i = vec.begin(); i < vec.end(); i++) {
                os << i << " " << vec[i] << std::endl;
              }
              os << std::endl;
              vec.clearBroadcast();
            }
          }
        }
        return os;
      }

    private:

      std::string _status;
      bool _serial;
      bool _vecIsAllocated;

      unsigned _iproc;
      unsigned _nprocs;
      MPI_Datatype _MY_MPI_DATATYPE;

      unsigned _begin;
      unsigned _end;
      unsigned _size;

      std::vector< Type > _vec;
      std::vector< Type > _vec2;
      std::vector < unsigned > _offset;

      unsigned _lproc;
  };





} //end namespace femus



#endif
