/*=========================================================================

Program: FEMUS
Module: MultiLevelMesh
Authors: Eugenio Aulisa, Simone Bnà

Copyright (c) FEMUS
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __MultiLevelMesh_hpp__
#define __MultiLevelMesh_hpp__
#include <vector>

#include "ElemTypeEnum.hpp"


namespace femus {


//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem_type;
class mesh;


/**
* This class is a black box container to handle multilevel mesh.
*/

class MultiLevelMesh {

public:

    /** Constructors */
    MultiLevelMesh() {};

    MultiLevelMesh(const unsigned short &igridn,const unsigned short &igridr,
                   const char mesh_file[], const char GaussOrder[], const double Lref,
                   bool (* SetRefinementFlag)(const double &x, const double &y, const double &z,
                           const int &ElemGroupNumber,const int &level));
    
    /** Destructor */
    ~MultiLevelMesh();

    /** Read the coarse-mesh from an input file (call the right reader from the extension) */
    void ReadCoarseMesh(const char mesh_file[], const char GaussOrder[], const double Lref);

    /** Built-in cube-structured mesh generator */
    void BuildBrickCoarseMesh( const unsigned int nx,
                               const unsigned int ny,
                               const unsigned int nz,
                               const double xmin, const double xmax,
                               const double ymin, const double ymax,
                               const double zmin, const double zmax,
                               const ElemType type,
                               const char GaussOrder[]
                             );

    /** Built-in square-structured mesh generator */
    void BuildRectangleCoarseMesh( const unsigned int nx,
                                   const unsigned int ny,
                                   const double xmin, const double xmax,
                                   const double ymin, const double ymax,
                                   const ElemType type,
                                   const char GaussOrder[]
                                 );

    /** Refine the coarse mesh (totally or selectively (in according to the SetRefinementFlag user-function) ) */
    void RefineMesh(const unsigned short &igridn, const unsigned short &igridr,
                    bool (* SetRefinementFlag)(const double &x, const double &y, const double &z,
                            const int &ElemGroupNumber,const int &level));

    /** Add a partially refined mesh level in the AMR alghorithm **/
    void AddMeshLevel( bool (* SetRefinementFlag)(const double &x, const double &y, const double &z,
                      const int &ElemGroupNumber,const int &level));
    
    
    /** Get the mesh pointer to level i */
    mesh* GetLevel(const unsigned i) {
        return _level[i];
    };

    /** Get the number of grid */
    unsigned GetNumberOfGrid() {
        return _gridn;
    };

    /** Get the Number of grid totally uniformly refined */
    unsigned GetNumberOfGridTotallyRefined() {
        return _gridr;
    };

    /** Erase levels_to_be_erased levels from the mesh array */
    void EraseCoarseLevels(unsigned levels_to_be_erased);

    /** Mark the node as solid or fluid */
    void MarkStructureNode();

    /** Print the mesh info for each level */
    void print_info();

    // data
    const elem_type *_type_elem[6][5];

protected:

private:

    // data
    unsigned short _gridn0, _gridr0;
    unsigned short _gridn, _gridr;
    /** Array of mesh */
    std::vector <mesh*> _level0;
    std::vector <mesh*> _level;
    std::vector <bool> _type_elem_flag;

};


} //end namespace femus



#endif
