/*=========================================================================

 Program: FEMUS
 Module: MultiLevelMesh
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelMesh.hpp"
#include "ElemType.hpp"
#include "Elem.hpp"
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "FemusConfig.hpp"
#include "MeshRefinement.hpp"
#include "Domain.hpp"


//C++ include
#include <iostream>


namespace femus {

using std::cout;
using std::endl;

MultiLevelMesh::~MultiLevelMesh() {

    for (unsigned i=0; i<_level0.size(); i++) {
        delete _level0[i];
    }

    for(unsigned i=0;i<6;i++){
      if( _finiteElementGeometryFlag[i])
      for(unsigned j=0;j<5;j++){
	delete _finiteElement[i][j];
      }
    }
    
    if(_writer != NULL) delete _writer;

}

//---------------------------------------------------------------------------------------------------
MultiLevelMesh::MultiLevelMesh(): _gridn0(0)
  {

  _finiteElementGeometryFlag.resize(6,false);

  for(int i=0; i<6; i++) {
    for(int j=0; j<5; j++) {
      _finiteElement[i][j] = NULL;
    }
  }
  _writer = NULL;

}


  void MultiLevelMesh::BuildElemType(const char GaussOrder[]){
    if(_finiteElementGeometryFlag[0]) {
      _finiteElement[0][0]=new const elem_type_3D("hex","linear",GaussOrder);
      _finiteElement[0][1]=new const elem_type_3D("hex","quadratic",GaussOrder);
      _finiteElement[0][2]=new const elem_type_3D("hex","biquadratic",GaussOrder);
      _finiteElement[0][3]=new const elem_type_3D("hex","constant",GaussOrder);
      _finiteElement[0][4]=new const elem_type_3D("hex","disc_linear",GaussOrder);
    }
    if(_finiteElementGeometryFlag[1]) {
      _finiteElement[1][0]=new const elem_type_3D("tet","linear",GaussOrder);
      _finiteElement[1][1]=new const elem_type_3D("tet","quadratic",GaussOrder);
      _finiteElement[1][2]=new const elem_type_3D("tet","biquadratic",GaussOrder);
      _finiteElement[1][3]=new const elem_type_3D("tet","constant",GaussOrder);
      _finiteElement[1][4]=new const elem_type_3D("tet","disc_linear",GaussOrder);
    }
    if(_finiteElementGeometryFlag[2]) {
      _finiteElement[2][0]=new const elem_type_3D("wedge","linear",GaussOrder);
      _finiteElement[2][1]=new const elem_type_3D("wedge","quadratic",GaussOrder);
      _finiteElement[2][2]=new const elem_type_3D("wedge","biquadratic",GaussOrder);
      _finiteElement[2][3]=new const elem_type_3D("wedge","constant",GaussOrder);
      _finiteElement[2][4]=new const elem_type_3D("wedge","disc_linear",GaussOrder);
    }
    if(_finiteElementGeometryFlag[3]) {
      _finiteElement[3][0]=new const elem_type_2D("quad","linear",GaussOrder);
      _finiteElement[3][1]=new const elem_type_2D("quad","quadratic",GaussOrder);
      _finiteElement[3][2]=new const elem_type_2D("quad","biquadratic",GaussOrder);
      _finiteElement[3][3]=new const elem_type_2D("quad","constant",GaussOrder);
      _finiteElement[3][4]=new const elem_type_2D("quad","disc_linear",GaussOrder);
    }
    if(_finiteElementGeometryFlag[4]) {
      _finiteElement[4][0]=new const elem_type_2D("tri","linear",GaussOrder);
      _finiteElement[4][1]=new const elem_type_2D("tri","quadratic",GaussOrder);
      _finiteElement[4][2]=new const elem_type_2D("tri","biquadratic",GaussOrder);
      _finiteElement[4][3]=new const elem_type_2D("tri","constant",GaussOrder);
      _finiteElement[4][4]=new const elem_type_2D("tri","disc_linear",GaussOrder);
    }
    _finiteElementGeometryFlag[5]=1;
    _finiteElement[5][0]=new const elem_type_1D("line","linear",GaussOrder);
    _finiteElement[5][1]=new const elem_type_1D("line","quadratic",GaussOrder);
    _finiteElement[5][2]=new const elem_type_1D("line","biquadratic",GaussOrder);
    _finiteElement[5][3]=new const elem_type_1D("line","constant",GaussOrder);
    _finiteElement[5][4]=new const elem_type_1D("line","disc_linear",GaussOrder);
    _level0[0]->SetFiniteElementPtr(_finiteElement);
  }


//---------------------------------------------------------------------------------------------------
MultiLevelMesh::MultiLevelMesh(const unsigned short &igridn,const unsigned short &igridr,
			       const char mesh_file[], const char GaussOrder[], const double Lref,
			       bool (* SetRefinementFlag)(const std::vector < double > &x,
							  const int &ElemGroupNumber,const int &level) ):
    _gridn0(igridn)
    {


    _level0.resize(_gridn0);
    _finiteElementGeometryFlag.resize(5,false);

    //coarse mesh
    _level0[0] = new Mesh();
    std::cout << " Reading corse mesh from file: " << mesh_file << std::endl;
    _level0[0]->ReadCoarseMesh(mesh_file, Lref,_finiteElementGeometryFlag);

    BuildElemType(GaussOrder);

    //totally refined meshes
    for (unsigned i=1; i<igridr; i++) {
        MeshRefinement meshcoarser(*_level0[i-1u]);
        meshcoarser.FlagAllElementsToBeRefined();
	_level0[i] = new Mesh();
	MeshRefinement meshfiner(*_level0[i]);
        meshfiner.RefineMesh(i,_level0[i-1],_finiteElement);
    }

    if(SetRefinementFlag==NULL){
    }
    else{
      Mesh::_SetRefinementFlag = SetRefinementFlag;
      Mesh::_IsUserRefinementFunctionDefined = true;
    }


    //partially refined meshes
    for (unsigned i=igridr; i<_gridn0; i++) {
      if(!Mesh::_IsUserRefinementFunctionDefined) {
        cout << "Set Refinement Region flag is not defined! " << endl;
        exit(1);
      }
      else {
	MeshRefinement meshcoarser(*_level0[i-1u]);
        meshcoarser.FlagElementsToBeRefined();
      }
      _level0[i] = new Mesh();
      MeshRefinement meshfiner(*_level0[i]);
      meshfiner.RefineMesh(i,_level0[i-1],_finiteElement);
      //_level0[i]->RefineMesh(i,_level0[i-1],_finiteElement);
    }

    unsigned refindex = _level0[0]->GetRefIndex();
    elem_type::_refindex=refindex;

    _gridn=_gridn0;
    _level.resize(_gridn);
    for(int i=0; i<_gridn; i++)
        _level[i]=_level0[i];
    
    _writer = NULL;

}

void MultiLevelMesh::ReadCoarseMesh(const char mesh_file[], const char GaussOrder[], const double Lref)
{
    
  const bool read_groups = true; //by default groups are read
  const bool read_boundary_groups = true; //by default boundary groups are read
  
    ReadCoarseMesh(mesh_file, GaussOrder, Lref, read_groups, read_boundary_groups);
   
}


void MultiLevelMesh::ReadCoarseMesh(const char mesh_file[], const char GaussOrder[], const double Lref, const bool read_groups, const bool read_boundary_groups)
{
    _gridn0 = 1;

    _level0.resize(_gridn0);
    _finiteElementGeometryFlag.resize(5,false);

    //coarse mesh
    _level0[0] = new Mesh();
    std::cout << " Reading corse mesh from file: " << mesh_file << std::endl;
    _level0[0]->ReadCoarseMesh(mesh_file, Lref,_finiteElementGeometryFlag, read_groups, read_boundary_groups);

    BuildElemType(GaussOrder);

    _gridn=_gridn0;
    _level.resize(_gridn);
    _level[0] = _level0[0];

}



void MultiLevelMesh::GenerateCoarseBoxMesh(
        const unsigned int nx, const unsigned int ny, const unsigned int nz,
        const double xmin, const double xmax,
        const double ymin, const double ymax,
        const double zmin, const double zmax,
        const ElemType type, const char GaussOrder[])
{
    _gridn0 = 1;

    _level0.resize(_gridn0);
    _finiteElementGeometryFlag.resize(5,false);

    //coarse mesh
    _level0[0] = new Mesh();
    std::cout << " Building brick mesh using the built-in mesh generator" << std::endl;

    _level0[0]->GenerateCoarseBoxMesh(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,type,_finiteElementGeometryFlag);

    BuildElemType(GaussOrder);

    _gridn=_gridn0;
    _level.resize(_gridn);
    _level[0] = _level0[0];

}


void MultiLevelMesh::RefineMesh( const unsigned short &igridn, const unsigned short &igridr,
                                 bool (* SetRefinementFlag)(const std::vector < double >& x,
                                         const int &ElemGroupNumber,const int &level))
{

    _gridn0 = igridn;

    _level0.resize(_gridn0);

    //totally refined meshes
    for (unsigned i=1; i<igridr; i++) {
      MeshRefinement meshcoarser(*_level0[i-1u]);
      meshcoarser.FlagAllElementsToBeRefined();

      _level0[i] = new Mesh();
      MeshRefinement meshfiner(*_level0[i]);
      meshfiner.RefineMesh(i,_level0[i-1],_finiteElement);
    }

    //partially refined meshes

    if(SetRefinementFlag==NULL){

    }
    else{
      Mesh::_SetRefinementFlag = SetRefinementFlag;
      Mesh::_IsUserRefinementFunctionDefined = true;
    }

    for (unsigned i=igridr; i<_gridn0; i++) {
      if(Mesh::_IsUserRefinementFunctionDefined == false) {
        cout << "Set Refinement Region flag is not defined! " << endl;
        exit(1);
      }
      else {
	MeshRefinement meshcoarser(*_level0[i-1u]);
        meshcoarser.FlagElementsToBeRefined();
	//meshcoarser.FlagOnlyEvenElementsToBeRefined();
      }
      _level0[i] = new Mesh();
      MeshRefinement meshfiner(*_level0[i]);
      meshfiner.RefineMesh(i,_level0[i-1],_finiteElement);
    }

    unsigned refindex = _level0[0]->GetRefIndex();
    elem_type::_refindex=refindex;


    _gridn=_gridn0;
    _level.resize(_gridn);
    for(int i=0; i<_gridn; i++)
        _level[i]=_level0[i];

}

void MultiLevelMesh::AddAMRMeshLevel()
{

  //AMR refine mesh
   _level0.resize(_gridn0+1u);

  MeshRefinement meshcoarser(*_level0[_gridn0-1u]);
  meshcoarser.FlagElementsToBeRefined();

  _level0[_gridn0] = new Mesh();
  MeshRefinement meshfiner(*_level0[_gridn0]);
  meshfiner.RefineMesh(_gridn0,_level0[_gridn0-1u],_finiteElement);

  _level.resize(_gridn+1u);
  _level[_gridn]=_level0[_gridn0];

  _gridn0++;
  _gridn++;
}


//---------------------------------------------------------------------------------------------

void MultiLevelMesh::EraseCoarseLevels(unsigned levels_to_be_erased) {
    _gridn -= levels_to_be_erased;
    for(int i=0; i<_gridn; i++) {
      _level[i]=_level0[i+levels_to_be_erased];
      _level[i]->SetLevel(i);
    }
}

//---------------------------------------------------------------------------------------------

void MultiLevelMesh::PrintInfo() {
    std::cout << " Number of uniform mesh refinement: " << _gridn << std::endl;
    for(int i=0; i<_gridn; i++) {
      _level[i]->PrintInfo();
    }
}

    /** Get the dimension of the problem (1D, 2D, 3D) from one Mesh (level 0 always exists, after initialization) */
    const unsigned MultiLevelMesh::GetDimension() const {
      return _level0[LEV_PICK]->GetDimension();
    }

// ========================================================
  void MultiLevelMesh::SetDomain(Domain* domain_in)  {

    _domain = domain_in;

   return;
  }

 // ========================================================
  Domain* MultiLevelMesh::GetDomain() const {

   return _domain;

  }

} //end namespace femus


