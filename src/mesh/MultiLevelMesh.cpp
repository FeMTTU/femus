
//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelMesh.hpp"
#include "ElemType.hpp"
#include "Elem.hpp"
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "FEMTTUConfig.h"


//C++ include
#include <iostream>


namespace femus {


using std::cout;
using std::endl;

MultiLevelMesh::~MultiLevelMesh() {

    for (unsigned i=0; i<_gridn0; i++) {
        delete _level0[i];
    }

    for (unsigned i=0; i<3; i++) {
        if(_type_elem_flag[i]) {
            for (unsigned j=0; j<3; j++)
                if (1!=i || 2!=j) delete _type_elem[i][j];
        }
    }

    for (unsigned i=3; i<5; i++) {
        if(_type_elem_flag[i]) {
            for (unsigned j=0; j<3; j++)
                if (4!=i || 2!=j) delete _type_elem[i][j];
        }
    }

    if(_type_elem_flag[0]) {
        delete _type_elem[0][3];
        delete _type_elem[0][4];
    }
    if(_type_elem_flag[3]) {
        delete _type_elem[3][3];
        delete _type_elem[3][4];
    }

    delete _type_elem[5][0];
    delete _type_elem[5][1];

};

//---------------------------------------------------------------------------------------------------
MultiLevelMesh::MultiLevelMesh(const unsigned short &igridn,const unsigned short &igridr, const char mesh_file[], const char GaussOrder[],
                               const double Lref, bool (* SetRefinementFlag)(const double &x, const double &y, const double &z,
                                       const int &ElemGroupNumber,const int &level)):
    _gridn0(igridn),
    _gridr0(igridr) {

    _level0.resize(_gridn0);
    _type_elem_flag.resize(5,false);
    
    //coarse mesh
    _level0[0] = new mesh();
    std::cout << " Reading corse mesh from file: " << mesh_file << std::endl;
    _level0[0]->ReadCoarseMesh(mesh_file, Lref,_type_elem_flag);

    if(_type_elem_flag[0]) {
        _type_elem[0][0]=new const elem_type("hex","linear",GaussOrder);
        _type_elem[0][1]=new const elem_type("hex","quadratic",GaussOrder);
        _type_elem[0][2]=new const elem_type("hex","biquadratic",GaussOrder);
        _type_elem[0][3]=new const elem_type("hex","constant",GaussOrder);
        _type_elem[0][4]=new const elem_type("hex","disc_linear",GaussOrder);
    }
    if(_type_elem_flag[1]) {
        _type_elem[1][0]=new const elem_type("tet","linear",GaussOrder);
        _type_elem[1][1]=new const elem_type("tet","biquadratic",GaussOrder);
        _type_elem[1][2]=_type_elem[1][1];
    }
    if(_type_elem_flag[2]) {
        _type_elem[2][0]=new const elem_type("wedge","linear",GaussOrder);
        _type_elem[2][1]=new const elem_type("wedge","quadratic",GaussOrder);
        _type_elem[2][2]=new const elem_type("wedge","biquadratic",GaussOrder);
    }
    if(_type_elem_flag[3]) {
        _type_elem[3][0]=new const elem_type("quad","linear",GaussOrder);
        _type_elem[3][1]=new const elem_type("quad","quadratic",GaussOrder);
        _type_elem[3][2]=new const elem_type("quad","biquadratic",GaussOrder);
        _type_elem[3][3]=new const elem_type("quad","constant",GaussOrder);
        _type_elem[3][4]=new const elem_type("quad","disc_linear",GaussOrder);
    }
    if(_type_elem_flag[4]) {
        _type_elem[4][0]=new const elem_type("tri","linear",GaussOrder);
        _type_elem[4][1]=new const elem_type("tri","biquadratic",GaussOrder);
        _type_elem[4][2]=_type_elem[4][1];
    }
    _type_elem[5][0]=new const elem_type("line","linear",GaussOrder);
    _type_elem[5][1]=new const elem_type("line","biquadratic",GaussOrder);
    _type_elem[5][2]=_type_elem[5][1];

    //totally refined meshes
    for (unsigned i=1; i<_gridr0; i++) {
        _level0[i-1u]->FlagAllElementsToBeRefined();
        _level0[i] = new mesh();
        _level0[i]->RefineMesh(i,_level0[i-1],_type_elem);
    }

    //partially refined meshes
    for (unsigned i=_gridr0; i<_gridn0; i++) {
        if(SetRefinementFlag==NULL) {
            cout << "Set Refinement Region flag is not defined! " << endl;
            exit(1);
        }
        else {
            mesh::_SetRefinementFlag = SetRefinementFlag;
            _level0[i-1u]->FlagElementsToBeRefinedByUserDefinedFunction();
        }
        _level0[i] = new mesh();
        _level0[i]->RefineMesh(i,_level0[i-1],_type_elem);

    }

    unsigned refindex = _level0[0]->GetRefIndex();
    elem_type::_refindex=refindex;


    _gridn=_gridn0;
    _gridr=_gridr0;
    _level.resize(_gridn);
    for(int i=0; i<_gridn; i++)
        _level[i]=_level0[i];

}

void MultiLevelMesh::ReadCoarseMesh(const char mesh_file[], const char GaussOrder[], const double Lref)
{
    _gridn0 = 1;
    _gridr0 = 1;

    _level0.resize(_gridn0);
    _type_elem_flag.resize(5,false);
    
    //coarse mesh
    _level0[0] = new mesh();
    std::cout << " Reading corse mesh from file: " << mesh_file << std::endl;
    _level0[0]->ReadCoarseMesh(mesh_file, Lref,_type_elem_flag);

    if(_type_elem_flag[0]) {
        _type_elem[0][0]=new const elem_type("hex","linear",GaussOrder);
        _type_elem[0][1]=new const elem_type("hex","quadratic",GaussOrder);
        _type_elem[0][2]=new const elem_type("hex","biquadratic",GaussOrder);
        _type_elem[0][3]=new const elem_type("hex","constant",GaussOrder);
        _type_elem[0][4]=new const elem_type("hex","disc_linear",GaussOrder);
    }
    if(_type_elem_flag[1]) {
        _type_elem[1][0]=new const elem_type("tet","linear",GaussOrder);
        _type_elem[1][1]=new const elem_type("tet","biquadratic",GaussOrder);
        _type_elem[1][2]=_type_elem[1][1];
    }
    if(_type_elem_flag[2]) {
        _type_elem[2][0]=new const elem_type("wedge","linear",GaussOrder);
        _type_elem[2][1]=new const elem_type("wedge","quadratic",GaussOrder);
        _type_elem[2][2]=new const elem_type("wedge","biquadratic",GaussOrder);
    }
    if(_type_elem_flag[3]) {
        _type_elem[3][0]=new const elem_type("quad","linear",GaussOrder);
        _type_elem[3][1]=new const elem_type("quad","quadratic",GaussOrder);
        _type_elem[3][2]=new const elem_type("quad","biquadratic",GaussOrder);
        _type_elem[3][3]=new const elem_type("quad","constant",GaussOrder);
        _type_elem[3][4]=new const elem_type("quad","disc_linear",GaussOrder);
    }
    if(_type_elem_flag[4]) {
        _type_elem[4][0]=new const elem_type("tri","linear",GaussOrder);
        _type_elem[4][1]=new const elem_type("tri","biquadratic",GaussOrder);
        _type_elem[4][2]=_type_elem[4][1];
    }
    _type_elem[5][0]=new const elem_type("line","linear",GaussOrder);
    _type_elem[5][1]=new const elem_type("line","biquadratic",GaussOrder);
    _type_elem[5][2]=_type_elem[5][1];

    _gridn=_gridn0;
    _gridr=_gridr0;
    _level.resize(_gridn);
    _level[0] = _level0[0];

}

void MultiLevelMesh::BuildBrickCoarseMesh( const unsigned int nx,
        const unsigned int ny,
        const unsigned int nz,
        const double xmin, const double xmax,
        const double ymin, const double ymax,
        const double zmin, const double zmax,
        const ElemType type,
        const char GaussOrder[]
                                         )
{
    _gridn0 = 1;
    _gridr0 = 1;

    _level0.resize(_gridn0);
    _type_elem_flag.resize(5,false);
    
    //coarse mesh
    _level0[0] = new mesh();
    std::cout << " Building brick mesh using the built-in mesh generator" << std::endl;
    _level0[0]->BuildBrick(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,type,_type_elem_flag);

    if(_type_elem_flag[0]) {
        _type_elem[0][0]=new const elem_type("hex","linear",GaussOrder);
        _type_elem[0][1]=new const elem_type("hex","quadratic",GaussOrder);
        _type_elem[0][2]=new const elem_type("hex","biquadratic",GaussOrder);
        _type_elem[0][3]=new const elem_type("hex","constant",GaussOrder);
        _type_elem[0][4]=new const elem_type("hex","disc_linear",GaussOrder);
    }
    if(_type_elem_flag[1]) {
        _type_elem[1][0]=new const elem_type("tet","linear",GaussOrder);
        _type_elem[1][1]=new const elem_type("tet","biquadratic",GaussOrder);
        _type_elem[1][2]=_type_elem[1][1];
    }
    if(_type_elem_flag[2]) {
        _type_elem[2][0]=new const elem_type("wedge","linear",GaussOrder);
        _type_elem[2][1]=new const elem_type("wedge","quadratic",GaussOrder);
        _type_elem[2][2]=new const elem_type("wedge","biquadratic",GaussOrder);
    }
    if(_type_elem_flag[3]) {
        _type_elem[3][0]=new const elem_type("quad","linear",GaussOrder);
        _type_elem[3][1]=new const elem_type("quad","quadratic",GaussOrder);
        _type_elem[3][2]=new const elem_type("quad","biquadratic",GaussOrder);
        _type_elem[3][3]=new const elem_type("quad","constant",GaussOrder);
        _type_elem[3][4]=new const elem_type("quad","disc_linear",GaussOrder);
    }
    if(_type_elem_flag[4]) {
        _type_elem[4][0]=new const elem_type("tri","linear",GaussOrder);
        _type_elem[4][1]=new const elem_type("tri","biquadratic",GaussOrder);
        _type_elem[4][2]=_type_elem[4][1];
    }
    _type_elem[5][0]=new const elem_type("line","linear",GaussOrder);
    _type_elem[5][1]=new const elem_type("line","biquadratic",GaussOrder);
    _type_elem[5][2]=_type_elem[5][1];

    _gridn=_gridn0;
    _gridr=_gridr0;
    _level.resize(_gridn);
    _level[0] = _level0[0];

}


void MultiLevelMesh::BuildRectangleCoarseMesh( const unsigned int nx,
        const unsigned int ny,
        const double xmin, const double xmax,
        const double ymin, const double ymax,
        const ElemType type,
        const char GaussOrder[]
                                             )
{

    BuildBrickCoarseMesh(nx,ny,0,xmin,xmax,ymin,ymax,0.,0.,type,GaussOrder);

}


void MultiLevelMesh::RefineMesh( const unsigned short &igridn, const unsigned short &igridr,
                                 bool (* SetRefinementFlag)(const double &x, const double &y, const double &z,
                                         const int &ElemGroupNumber,const int &level))
{

    _gridn0 = igridn;
    _gridr0 = igridr;

    _level0.resize(_gridn0);

    //totally refined meshes
    for (unsigned i=1; i<_gridr0; i++) {
        _level0[i-1u]->FlagAllElementsToBeRefined();
        _level0[i] = new mesh();
        _level0[i]->RefineMesh(i,_level0[i-1],_type_elem);
    }

    //partially refined meshes
    for (unsigned i=_gridr0; i<_gridn0; i++) {
        if(SetRefinementFlag==NULL) {
            cout << "Set Refinement Region flag is not defined! " << endl;
            exit(1);
        }
        else {
            mesh::_SetRefinementFlag = SetRefinementFlag;
            _level0[i-1u]->FlagElementsToBeRefinedByUserDefinedFunction();
        }
        _level0[i] = new mesh();
        _level0[i]->RefineMesh(i,_level0[i-1],_type_elem);

    }

    unsigned refindex = _level0[0]->GetRefIndex();
    elem_type::_refindex=refindex;


    _gridn=_gridn0;
    _gridr=_gridr0;
    _level.resize(_gridn);
    for(int i=0; i<_gridn; i++)
        _level[i]=_level0[i];

}


void MultiLevelMesh::AddMeshLevel( bool (* SetRefinementFlag)(const double &x, const double &y, const double &z,
                                  const int &ElemGroupNumber,const int &level))
{
 
  _level0.resize(_gridn0+1u);

    
  //AMR refine mesh
            
  if(SetRefinementFlag==NULL) {
    cout << "Set Refinement Region flag is not defined! " << endl;
    exit(1);
  }
  
  mesh::_SetRefinementFlag = SetRefinementFlag;
  _level0[_gridn0-1u]->FlagElementsToBeRefinedByUserDefinedFunction();
  
  _level0[_gridn0] = new mesh();
  _level0[_gridn0]->RefineMesh(_gridn0,_level0[_gridn0-1u],_type_elem);
    
  _level.resize(_gridn+1u);
  _level[_gridn]=_level0[_gridn0];

  _gridn0++;
  _gridn++;
}




//---------------------------------------------------------------------------------------------

void MultiLevelMesh::EraseCoarseLevels(unsigned levels_to_be_erased) {
    if(levels_to_be_erased >= _gridr) {
        levels_to_be_erased = _gridr-1;
        cout<<"Warning the number of levels to be erased has been reduced to"<<levels_to_be_erased;
    }
    _gridr -= levels_to_be_erased;
    _gridn -= levels_to_be_erased;
    for(int i=0; i<_gridn; i++) {
        _level[i]=_level0[i+levels_to_be_erased];
        _level[i]->SetGridNumber(i);
    }
}

//---------------------------------------------------------------------------------------------

void MultiLevelMesh::MarkStructureNode() {
    for (unsigned i=0; i<_gridn0; i++) _level0[i]->AllocateAndMarkStructureNode();
}


//---------------------------------------------------------------------------------------------

void MultiLevelMesh::print_info() {
    std::cout << " Number of uniform mesh refinement: " << _gridn << std::endl;
    for(int i=0; i<_gridn; i++) {
        _level[i]->print_info();
    }
}


} //end namespace femus


