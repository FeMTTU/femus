 /*=========================================================================

 Program: FEMUS
 Module: MultiLevelProblem
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __MultiLevelMesh_hpp__
#define __MultiLevelMesh_hpp__
#include <vector>
//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem_type;
class mesh;


/**
* This class is a black box container to handle multilevel mesh.
*/

class MultiLevelMesh {

private:
 
  unsigned short _gridn, _gridr;
  /** Array of mesh */
  std::vector <mesh*> _level;
  std::vector <bool> _type_elem_flag;
protected:
  
  
 public:
     const elem_type *_type_elem[6][5]; 
  
 
  /** Constructor */
  MultiLevelMesh(const unsigned short &igridn,const unsigned short &igridr,
		 const char mesh_file[], const char GaussOrder[], const double Lref, 
		 bool (* SetRefinementFlag)(const double &x, const double &y, const double &z, 
					    const int &ElemGroupNumber,const int &level));

  /** Destructor */
  ~MultiLevelMesh();

  mesh* GetLevel(const unsigned i) {return _level[i];};
  
  unsigned GetNumberOfGrid();
  unsigned GetNumberOfGridTotallyRefined();
};

#endif