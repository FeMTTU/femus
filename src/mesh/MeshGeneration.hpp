/*=========================================================================

 Program: FEMuS
 Module:  MeshGeneration
 Authors: Simone Bnà

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __meshgeneration_hpp__
#define __meshgeneration_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------


namespace femus {

//using std::vector;
class mesh;

// ------------------------------------------------------------
// MeshTools::Generation namespace
namespace MeshTools
{
  /**
   * Tools for \p Mesh generation.
   */
  
  namespace Generation
  {

    /** Built-in cube-unstructured mesh generator */
    void BuildBrick ( mesh& mesh, 
                      const unsigned int nx,
                      const unsigned int ny,
                      const unsigned int nz,
                      const double xmin, const double xmax,
                      const double ymin, const double ymax,
                      const double zmin, const double zmax,
                      const ElemType type,
                      std::vector<bool> &type_elem_flag );


    /** 2D local-to-global map for built-in mesh generator */
    unsigned int idx(const ElemType type, const unsigned int nx, const unsigned int i, const unsigned int j);

    /** 3D local-to-global map for built-in mesh generator */
    unsigned int idx(const ElemType type, const unsigned int nx, const unsigned int ny, const unsigned int i,
                     const unsigned int j, const unsigned int k);

  }
  
}

} //end namespace femus



#endif