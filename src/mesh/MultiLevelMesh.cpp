
//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelMesh.hpp"
#include "ExplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "System.hpp"
#include "ElemType.hpp"
#include "Elem.hpp"
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"
#include "FEMTTUConfig.h"
#include "Parameter.hpp"
#include "hdf5.h"


//C++ include
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>   
#include <string.h>
using std::cout;
using std::endl;
using std::min;
using std::string;

MultiLevelMesh::~MultiLevelMesh() {

  for (unsigned i=0; i<_gridn; i++) {
    delete _level[i];
  }

  for (unsigned i=0; i<3; i++)
    for (unsigned j=0; j<3; j++)
      if (1!=i || 2!=j) delete type_elem[i][j];

  for (unsigned i=3; i<5; i++)
    for (unsigned j=0; j<3; j++)
      if (4!=i || 2!=j) delete type_elem[i][j];

  delete type_elem[5][0];
  delete type_elem[5][1];

  delete type_elem[0][3];
  delete type_elem[0][4];

  delete type_elem[3][3];
  delete type_elem[3][4];
};

// //---------------------------------------------------------------------------------------------------
MultiLevelMesh::MultiLevelMesh(const unsigned short &igridn,const unsigned short &igridr, const char mesh_file[], const char GaussOrder[],
						       const double Lref, bool (* SetRefinementFlag)(const double &x, const double &y, const double &z, 
												     const int &ElemGroupNumber,const int &level)):
			       _gridn(igridn), 
			       _gridr(igridr){
		       
  type_elem[0][0]=new const elem_type("hex","linear",GaussOrder);
  type_elem[0][1]=new const elem_type("hex","quadratic",GaussOrder);
  type_elem[0][2]=new const elem_type("hex","biquadratic",GaussOrder);
  type_elem[0][3]=new const elem_type("hex","constant",GaussOrder);
  type_elem[0][4]=new const elem_type("hex","disc_linear",GaussOrder);


  type_elem[1][0]=new const elem_type("tet","linear",GaussOrder);
  type_elem[1][1]=new const elem_type("tet","biquadratic",GaussOrder);
  type_elem[1][2]=type_elem[1][1];

  type_elem[2][0]=new const elem_type("wedge","linear",GaussOrder);
  type_elem[2][1]=new const elem_type("wedge","quadratic",GaussOrder);
  type_elem[2][2]=new const elem_type("wedge","biquadratic",GaussOrder);

  type_elem[3][0]=new const elem_type("quad","linear",GaussOrder);
  type_elem[3][1]=new const elem_type("quad","quadratic",GaussOrder);
  type_elem[3][2]=new const elem_type("quad","biquadratic",GaussOrder);
  type_elem[3][3]=new const elem_type("quad","constant",GaussOrder);
  type_elem[3][4]=new const elem_type("quad","disc_linear",GaussOrder);


  type_elem[4][0]=new const elem_type("tri","linear",GaussOrder);
  type_elem[4][1]=new const elem_type("tri","biquadratic",GaussOrder);
  type_elem[4][2]=type_elem[4][1];

  type_elem[5][0]=new const elem_type("line","linear",GaussOrder);
  type_elem[5][1]=new const elem_type("line","biquadratic",GaussOrder);
  type_elem[5][2]=type_elem[5][1];

  cout << "MESH DATA: " << endl;

  _level.resize(_gridn);
    
  //coarse mesh
  
  _level[0]=new mesh(mesh_file, Lref);
 
    
  //totally refined meshes
  for (unsigned i=1; i<_gridr; i++) {
    _level[i-1u]->_coordinate->SetElementRefiniement(1);
    _level[i] = new mesh(i,_level[i-1],type_elem); 
  }
  
  //partially refined meshes
  for (unsigned i=_gridr; i<_gridn; i++) {
    if(SetRefinementFlag==NULL) {
      cout << "Set Refinement Region flag is not defined! " << endl;
      exit(1);
    }
    else {
      mesh::_SetRefinementFlag = SetRefinementFlag;
      _level[i-1u]->_coordinate->SetElementRefiniement(2);
    }
    _level[i] = new mesh(i,_level[i-1],type_elem); 
  }
  _level[_gridn-1u]->_coordinate->SetElementRefiniement(0);
      
  unsigned refindex = _level[0]->GetRefIndex();
  elem_type::_refindex=refindex;
    
}