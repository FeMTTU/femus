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

#include "ElemTypeEnum.hpp"


#include <vector>
#include <iostream>

namespace femus {


class Mesh;

// ------------------------------------------------------------
// MeshTools::Generation namespace
namespace MeshTools
{
  /**
   * Tools for \p Mesh generation.
   */
  
  class Generation
  {

  public:

    /** Built-in cube-unstructured mesh generator */
    static void BuildBox ( Mesh& mesh,
		      std::vector < std::vector < double> > &coords,
                      const unsigned int nx,
                      const unsigned int ny,
                      const unsigned int nz,
                      const double xmin, const double xmax,
                      const double ymin, const double ymax,
                      const double zmin, const double zmax,
                      const ElemType type,
                      std::vector<bool> &type_elem_flag );

  private:

    static  inline unsigned int idx(const ElemType type, const unsigned int nx, const unsigned int i, const unsigned int j);

    static  inline unsigned int idx(const ElemType type,
                         const unsigned int nx,
                         const unsigned int ny,
                         const unsigned int i,
                         const unsigned int j,
                         const unsigned int k);

  };



        /**
         * A useful inline function which replaces the #defines
         * used previously.  Not private since this is a namespace,
         * but would be if this were a class.  The first one returns
         * the proper node number for 2D elements while the second
         * one returns the node number for 3D elements.
         */
        inline
        unsigned int Generation::idx(const ElemType type, const unsigned int nx, const unsigned int i, const unsigned int j)
        {

          switch(type) {
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	  case TRI3:
// 	    {
// 	      return i + j*(nx+1);
// 	      break;
// 	    }
//
// 	  case QUAD8:
            case QUAD9:
            case TRI6: {
                return i + j * (2 * nx + 1);
                break;
              }

            default: {
                std::cout << "ERROR: Unrecognized or Not Supported 2D element type." << std::endl;
                exit(1);
              }
          }

          return -1; // invalid_uint
        }


        // Same as the function above, but for 3D elements
        inline
        unsigned int Generation::idx(const ElemType type,
                         const unsigned int nx,
                         const unsigned int ny,
                         const unsigned int i,
                         const unsigned int j,
                         const unsigned int k)
        {
          switch(type) {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	  case PRISM6:
// 	    {
// 	      return i + (nx+1)*(j + k*(ny+1));
// 	      break;
// 	    }

// 	  case HEX20:
            case HEX27:
// 	  case TET4:  // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
// 	  case PRISM15:
// 	  case PRISM18:
              {
                return i + (2 * nx + 1) * (j + k * (2 * ny + 1));
                break;
              }

            default: {
                std::cout << "ERROR: Unrecognized element type." << std::endl;
                exit(1);
              }
          }

          return -1;
        }




  
}

} //end namespace femus



#endif
