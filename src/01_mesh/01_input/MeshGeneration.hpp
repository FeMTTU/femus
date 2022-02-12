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

#ifndef __femus_mesh_MeshGeneration_hpp__
#define __femus_mesh_MeshGeneration_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------


namespace femus {


class Mesh;

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
    void BuildBox ( Mesh& mesh,
		      vector < vector < double> > &coords,
                      const unsigned int nx,
                      const unsigned int ny,
                      const unsigned int nz,
                      const double xmin, const double xmax,
                      const double ymin, const double ymax,
                      const double zmin, const double zmax,
                      const ElemType type,
                      std::vector<bool> &type_elem_flag );

  }
  
}

} //end namespace femus



#endif