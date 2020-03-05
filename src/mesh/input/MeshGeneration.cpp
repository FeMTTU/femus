/*=========================================================================

 Program: FEMUS
 Module: MeshGeneration
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

// #include <iostream>
// #include <fstream>
// #include <cstdlib>
// #include <cmath>
// #include <cstring>
#include <assert.h>
#include "mpi.h"
// #include <algorithm>
#include "Mesh.hpp"
#include "MeshGeneration.hpp"



namespace femus
{


  namespace MeshTools
  {
    namespace Generation
    {
      namespace Private
      {

        /**
         * A useful inline function which replaces the #defines
         * used previously.  Not private since this is a namespace,
         * but would be if this were a class.  The first one returns
         * the proper node number for 2D elements while the second
         * one returns the node number for 3D elements.
         */
        inline
        unsigned int idx(const ElemType type, const unsigned int nx, const unsigned int i, const unsigned int j)
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
        unsigned int idx(const ElemType type,
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

// ------------------------------------------------------------
// MeshTools::Generation function for mesh generation
      void BuildBox(Mesh& mesh,
                    vector < vector < double> >& vt,
                    const unsigned int nx,
                    const unsigned int ny,
                    const unsigned int nz,
                    const double xmin, const double xmax,
                    const double ymin, const double ymax,
                    const double zmin, const double zmax,
                    const ElemType type,
                    std::vector<bool>& type_elem_flag)
      {

        using namespace MeshTools::Generation::Private;

        // Clear the mesh and start from scratch
        //mesh.clear(); // to be added

//   vector <vector <double> > vt;
//   vt.resize(3);

//   mesh.SetGridNumber(0);

        if(nz != 0)
          mesh.SetDimension(3);
        else if(ny != 0)
          mesh.SetDimension(2);
        else if(nx != 0)
          mesh.SetDimension(1);
        else
          mesh.SetDimension(0);

        unsigned ngroup;
        unsigned nbcd;

        switch(mesh.GetDimension()) {

            //---------------------------------------------------------------------
            // Build a 1D line
          case 1: {
              assert(nx != 0);
              assert(ny == 0);
              assert(nz == 0);
              assert(xmin < xmax);

              // Reserve elements
              switch(type) {
//           case INVALID_ELEM:
//           case EDGE2:
                case EDGE3:
//           case EDGE4:
                  {
                    mesh.SetNumberOfElements(nx);

                    ngroup = 1;
                    nbcd   = 2;
                    break;
                  }
//
                default: {
                    std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
                  }
              }

              // Reserve nodes
              switch(type) {
//           case INVALID_ELEM:
//           case EDGE2:
//             {
//               mesh.reserve_nodes(nx+1);
//               break;
//             }
//
                case EDGE3: {
                    mesh.SetNumberOfNodes((2 * nx + 1));
                    break;
                  }
//
//           case EDGE4:
//             {
//               mesh.reserve_nodes(3*nx+1);
//               break;
//             }
//
                default: {
                    std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
                  }
              }


              // Build the nodes, depends on whether we're using linears,
              // quadratics or cubics and whether using uniform grid or Gauss-Lobatto
              unsigned int node_id = 0;
              switch(type) {
//           case INVALID_ELEM:
//           case EDGE2:
//             {
//               for (unsigned int i=0; i<=nx; i++)
//               {
//                 if (gauss_lobatto_grid)
//                   mesh.add_point (Point(0.5*(std::cos(libMesh::pi*static_cast<double>(nx-i)/static_cast<double>(nx))+1.0),
//                         0,
//                         0), node_id++);
//                 else
//                   mesh.add_point (Point(static_cast<double>(i)/static_cast<double>(nx),
//                         0,
//                         0), node_id++);
//               }
//               break;
//             }
//
                case EDGE3: {
                    vt[0].resize(mesh.GetNumberOfNodes());
                    vt[1].resize(mesh.GetNumberOfNodes());
                    vt[2].resize(mesh.GetNumberOfNodes());

                    for(unsigned int i = 0; i <= 2 * nx; i++) {
//                 if (gauss_lobatto_grid)
//                 {
//                   // The x location of the point.
//                   double x=0.;
//
//                   // Shortcut quantities (do not depend on i)
//                   const double c = std::cos( libMesh::pi*i / static_cast<double>(2*nx) );
//
//                   // If i is even, compute a normal Gauss-Lobatto point
//                   if (i%2 == 0)
//                     x = 0.5*(1.0 - c);
//
//                   // Otherwise, it is the average of the previous and next points
//                   else
//                   {
//                     double cmin = std::cos( libMesh::pi*(i-1) / static_cast<double>(2*nx) );
//                     double cmax = std::cos( libMesh::pi*(i+1) / static_cast<double>(2*nx) );
//
//                     double gl_xmin = 0.5*(1.0 - cmin);
//                     double gl_xmax = 0.5*(1.0 - cmax);
//                     x = 0.5*(gl_xmin + gl_xmax);
//                   }
//
//                   mesh.add_point (Point(x,0.,0.), node_id++);
//                 }
//                 else
                      vt[0][node_id] = (static_cast<double>(i) / static_cast<double>(2 * nx)) * (xmax - xmin) + xmin;
                      vt[1][node_id] = 0.;
                      vt[2][node_id] = 0.;

                      node_id++;

                    }
                    break;
                  }

                default: {
                    std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
                  }

              }

              // Build the elements of the mesh
              unsigned iel = 0;
              mesh.el = new elem(mesh.GetNumberOfElements());
              mesh.el->SetElementGroupNumber(1);
              // Build the elements.  Each one is a bit different.
              switch(type) {
//             case INVALID_ELEM:
//             case EDGE2:
//               {
//                 for (unsigned int i=0; i<nx; i++)
//                 {
//                   Elem* elem = mesh.add_elem (new Edge2);
//                   elem->set_node(0) = mesh.node_ptr(i);
//                   elem->set_node(1) = mesh.node_ptr(i+1);
//
//                   if (i == 0)
//                     mesh.boundary_info->add_side(elem, 0, 0);
//
//                   if (i == (nx-1))
//                     mesh.boundary_info->add_side(elem, 1, 1);
//                 }
//               break;
//               }
//
                case EDGE3: {


                    unsigned LocalToGlobalNodePerElement[3];

                    for(unsigned int i = 0; i < nx; i++) {

                      mesh.el->SetElementGroup(iel, 1);
                      mesh.el->SetElementMaterial(iel, 2);
                      type_elem_flag[5] = true;
                      mesh.el->AddToElementNumber(1, "Line");
                      mesh.el->SetElementType(iel, 5);

                      LocalToGlobalNodePerElement[0] = 2 * i   + 1;
                      LocalToGlobalNodePerElement[1] = 2 * i + 2 + 1;
                      LocalToGlobalNodePerElement[2] = 2 * i + 1 + 1;

                      // connectivity
                      for(unsigned iloc = 0; iloc < 3; iloc++) {
                        mesh.el->SetElementDofIndex(iel, iloc, LocalToGlobalNodePerElement[iloc] - 1u);
                      }

                      if(i == 0) {
                        mesh.el->SetFaceElementIndex(iel, 0, -2);
                        mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(0, "left"));
                      }

                      if(i == (nx - 1)) {
                        mesh.el->SetFaceElementIndex(iel, 1, -3);
                        mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(1, "right"));
                      }


                      iel++;
                    }

                    break;
                  }

                default: {
                    std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
                  }
              }

// 	// Scale the nodal positions
// 	for (unsigned int p=0; p<mesh.n_nodes(); p++)
// 	  mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
//
//         // Add sideset names to boundary info
//         mesh.boundary_info->sideset_name(0) = "left";
//         mesh.boundary_info->sideset_name(1) = "right";
//
//         // Add nodeset names to boundary info
//         mesh.boundary_info->nodeset_name(0) = "left";
//         mesh.boundary_info->nodeset_name(1) = "right";

              std::cout << "Error: NotImplemented " << std::endl;
              break;
            }

            //---------------------------------------------------------------------
            // Build a 2D quadrilateral
          case 2: {
              assert(nx != 0);
              assert(ny != 0);
              assert(nz == 0);
              assert(xmin < xmax);
              assert(ymin < ymax);

              switch(type) {
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	  case QUAD8:
                case QUAD9: {
                    mesh.SetNumberOfElements(nx * ny);

                    ngroup = 1;
                    nbcd   = 4;
                    break;
                  }
//
// 	  case TRI3:
                case TRI6: {
                    mesh.SetNumberOfElements(2 * nx * ny);

                    ngroup = 1;
                    nbcd   = 4;
                    break;
                  }
//
                default: {
                    std::cout << "ERROR: Unrecognized or NotSupported 2D element type." << std::endl;
                    exit(1);
                  }
              }

// 	// Reserve nodes.  The quadratic element types
// 	// need to reserve more nodes than the linear types.
              switch(type) {
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	  case TRI3:
// 	    {
// 	      mesh.reserve_nodes( (nx+1)*(ny+1) );
// 	      break;
// 	    }
//
// 	  case QUAD8:
                case QUAD9:
                case TRI6: {
                    mesh.SetNumberOfNodes((2 * nx + 1) * (2 * ny + 1));
                    break;
                  }


                default: {
                    std::cout << "ERROR: Unrecognized or not Supported 2D element type." << std::endl;
                    exit(1);
                  }
              }



              // Build the nodes. Depends on whether you are using a linear
              // or quadratic element, and whether you are using a uniform
              // grid [ or the Gauss-Lobatto grid points - NOTIMPLENTED ].
              unsigned int node_id = 0;
              switch(type) {
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	  case TRI3:
// 	    {
// 	      for (unsigned int j=0; j<=ny; j++)
// 		for (unsigned int i=0; i<=nx; i++)
// 		  {
// 		    if (gauss_lobatto_grid)
// 		      {
// 			mesh.add_point (Point(0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(i)/static_cast<double>(nx))),
// 					      0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(j)/static_cast<double>(ny))),
// 					      0.), node_id++);
// 		      }
//
// 		    else
// 		      mesh.add_point (Point(static_cast<double>(i)/static_cast<double>(nx),
// 					    static_cast<double>(j)/static_cast<double>(ny),
// 					    0.), node_id++);
// 		  }
//
// 	      break;
// 	    }
//
// 	  case QUAD8:
                case QUAD9:
                case TRI6: {
                    vt[0].resize(mesh.GetNumberOfNodes());
                    vt[1].resize(mesh.GetNumberOfNodes());
                    vt[2].resize(mesh.GetNumberOfNodes());


                    for(unsigned int j = 0; j <= (2 * ny); j++)
                      for(unsigned int i = 0; i <= (2 * nx); i++) {
// 		    if (gauss_lobatto_grid) // NOT YET IMPLENTED
// 		      {
// 			// The x,y locations of the point.
// 			double x=0., y=0.;
//
// 			// Shortcut quantities (do not depend on i,j)
// 			const double a = std::cos( libMesh::pi / static_cast<double>(2*nx) );
// 			const double b = std::cos( libMesh::pi / static_cast<double>(2*ny) );
//
// 			// Shortcut quantities (depend on i,j)
// 			const double c = std::cos( libMesh::pi*i / static_cast<double>(2*nx) );
// 			const double d = std::cos( libMesh::pi*j / static_cast<double>(2*ny) );
//
// 			// If i is even, compute a normal Gauss-Lobatto point
// 			if (i%2 == 0)
// 			  x = 0.5*(1.0 - c);
//
// 			// Otherwise, it is the average of the previous and next points
// 			else
// 			  x = 0.5*(1.0 - a*c);
//
// 			// If j is even, compute a normal Gauss-Lobatto point
// 			if (j%2 == 0)
// 			  y = 0.5*(1.0 - d);
//
// 			// Otherwise, it is the average of the previous and next points
// 			else
// 			  y = 0.5*(1.0 - b*d);
//
//
// 			mesh.add_point (Point(x,y,0.), node_id++);
// 		      }


// 		    else
                        vt[0][node_id] = (static_cast<double>(i) / static_cast<double>(2 * nx)) * (xmax - xmin) + xmin;
                        vt[1][node_id] = (static_cast<double>(j) / static_cast<double>(2 * ny)) * (ymax - ymin) + ymin;
                        vt[2][node_id] = 0.;

//			  std::cout << "inode: " << node_id <<  " x: " << vt[0][node_id] << "   y: " << vt[1][node_id] << std::endl;

                        node_id++;

                      }

                    break;
                  }

                default: {
                    std::cout << "ERROR: Unrecognized or NotSupported 2D element type." << std::endl;
                    exit(1);
                  }
              }


              unsigned iel = 0;
              mesh.el = new elem(mesh.GetNumberOfElements());
              mesh.el->SetElementGroupNumber(1);
              // Build the elements.  Each one is a bit different.
              switch(type) {
//
// 	  case INVALID_ELEM:
// 	  case QUAD4:
// 	    {
// 	      for (unsigned int j=0; j<ny; j++)
// 		for (unsigned int i=0; i<nx; i++)
// 		  {
// 		    Elem* elem = mesh.add_elem(new Quad4);
//
// 		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
// 		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j)  );
// 		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+1,j+1));
// 		    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i,j+1)  );
//
// 		    if (j == 0)
// 		      mesh.boundary_info->add_side(elem, 0, 0);
//
// 		    if (j == (ny-1))
// 		      mesh.boundary_info->add_side(elem, 2, 2);
//
// 		    if (i == 0)
// 		      mesh.boundary_info->add_side(elem, 3, 3);
//
// 		    if (i == (nx-1))
// 		      mesh.boundary_info->add_side(elem, 1, 1);
// 		  }
// 	      break;
// 	    }
//
//
// 	  case TRI3:
// 	    {
// 	      for (unsigned int j=0; j<ny; j++)
// 		for (unsigned int i=0; i<nx; i++)
// 		  {
// 		    Elem* elem = NULL;
//
// 		    // Add first Tri3
// 		    elem = mesh.add_elem(new Tri3);
//
// 		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
// 		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j)  );
// 		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+1,j+1));
//
// 		    if (j == 0)
// 		      mesh.boundary_info->add_side(elem, 0, 0);
//
// 		    if (i == (nx-1))
// 		      mesh.boundary_info->add_side(elem, 1, 1);
//
// 		    // Add second Tri3
// 		    elem = mesh.add_elem(new Tri3);
//
// 		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
// 		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j+1));
// 		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i,j+1)  );
//
// 		    if (j == (ny-1))
// 		      mesh.boundary_info->add_side(elem, 1, 2);
//
// 		    if (i == 0)
// 		      mesh.boundary_info->add_side(elem, 2, 3);
// 		  }
// 	      break;
// 	    }
//
//
//
// 	  case QUAD8:
                case QUAD9: {

                    unsigned LocalToGlobalNodePerElement[9];

                    for(unsigned int j = 0; j < (2 * ny); j += 2)
                      for(unsigned int i = 0; i < (2 * nx); i += 2) {

                        mesh.el->SetElementGroup(iel, 1);
                        mesh.el->SetElementMaterial(iel, 2);
                        type_elem_flag[3] = true;
                        mesh.el->AddToElementNumber(1, "Quad");
                        mesh.el->SetElementType(iel, 3);

                        LocalToGlobalNodePerElement[0] = idx(type, nx, i, j) + 1;
                        LocalToGlobalNodePerElement[1] = idx(type, nx, i + 2, j) + 1;
                        LocalToGlobalNodePerElement[2] = idx(type, nx, i + 2, j + 2) + 1;
                        LocalToGlobalNodePerElement[3] = idx(type, nx, i, j + 2) + 1;
                        LocalToGlobalNodePerElement[4] = idx(type, nx, i + 1, j) + 1;
                        LocalToGlobalNodePerElement[5] = idx(type, nx, i + 2, j + 1) + 1;
                        LocalToGlobalNodePerElement[6] = idx(type, nx, i + 1, j + 2) + 1;
                        LocalToGlobalNodePerElement[7] = idx(type, nx, i, j + 1) + 1;
                        LocalToGlobalNodePerElement[8] = idx(type, nx, i + 1, j + 1) + 1;

                        // connectivity
                        for(unsigned iloc = 0; iloc < 9; iloc++) {
                          mesh.el->SetElementDofIndex(iel, iloc, LocalToGlobalNodePerElement[iloc] - 1u);
                        }

                        if(j == 0) {
                          mesh.el->SetFaceElementIndex(iel, 0, -2);
                          mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(0, "bottom"));
                        }

                        if(j == 2 * (ny - 1)) {
                          mesh.el->SetFaceElementIndex(iel, 2, -4);
                          mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(2, "top"));
                        }

                        if(i == 0) {
                          mesh.el->SetFaceElementIndex(iel, 3, -5);
                          mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(3, "left"));
                        }

                        if(i == 2 * (nx - 1)) {
                          mesh.el->SetFaceElementIndex(iel, 1, -3);
                          mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(1, "right"));
                        }

                        iel++;
                      }

                    break;
                  }

                case TRI6: {

                    unsigned LocalToGlobalNodePerElement[6];

                    for(unsigned int j = 0; j < (2 * ny); j += 2)
                      for(unsigned int i = 0; i < (2 * nx); i += 2) {

                        // Add first Tri6
                        mesh.el->SetElementGroup(iel, 1);
                        mesh.el->SetElementMaterial(iel, 2); // 2 == fluid
                        type_elem_flag[4] = true;
                        mesh.el->AddToElementNumber(1, "Triangle");
                        mesh.el->SetElementType(iel, 4);


                        LocalToGlobalNodePerElement[0] = idx(type, nx, i, j) + 1;
                        LocalToGlobalNodePerElement[1] = idx(type, nx, i + 2, j) + 1;
                        LocalToGlobalNodePerElement[2] = idx(type, nx, i + 2, j + 2) + 1;
                        LocalToGlobalNodePerElement[3] = idx(type, nx, i + 1, j) + 1;
                        LocalToGlobalNodePerElement[4] = idx(type, nx, i + 2, j + 1) + 1;
                        LocalToGlobalNodePerElement[5] = idx(type, nx, i + 1, j + 1) + 1;

                        // connectivity
                        for(unsigned iloc = 0; iloc < 6; iloc++) {
                          mesh.el->SetElementDofIndex(iel, iloc, LocalToGlobalNodePerElement[iloc] - 1u);
                        }

                        if(j == 0) {
                          mesh.el->SetFaceElementIndex(iel, 0, -2);
                          mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(0, "bottom"));
                        }

                        if(i == 2 * (nx - 1)) {
                          mesh.el->SetFaceElementIndex(iel, 1, -3);
                          mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(1, "right"));
                        }

                        iel++;

                        // Add second Tri6

                        mesh.el->SetElementGroup(iel, 1);
                        mesh.el->SetElementMaterial(iel, 2); // 2 == fluid
                        type_elem_flag[4] = true;
                        mesh.el->AddToElementNumber(1, "Triangle");
                        mesh.el->SetElementType(iel, 4);

                        LocalToGlobalNodePerElement[0] = idx(type, nx, i, j) + 1;
                        LocalToGlobalNodePerElement[1] = idx(type, nx, i + 2, j + 2) + 1;
                        LocalToGlobalNodePerElement[2] = idx(type, nx, i, j + 2) + 1;
                        LocalToGlobalNodePerElement[3] = idx(type, nx, i + 1, j + 1) + 1;
                        LocalToGlobalNodePerElement[4] = idx(type, nx, i + 1, j + 2) + 1;
                        LocalToGlobalNodePerElement[5] = idx(type, nx, i, j + 1) + 1;

                        // connectivity
                        for(unsigned iloc = 0; iloc < 6; iloc++) {
                          mesh.el->SetElementDofIndex(iel, iloc, LocalToGlobalNodePerElement[iloc] - 1u);
                        }

                        if(j == 2 * (ny - 1)) {
                          mesh.el->SetFaceElementIndex(iel, 1, -4);
                          mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(2, "top"));
                        }

                        if(i == 0) {
                          mesh.el->SetFaceElementIndex(iel, 2, -5);
                          mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(3, "left"));
                        }

                        iel++;

                      }

                    break;
                  }


                default: {
                    std::cout << "ERROR: Unrecognized or Not Supported 2D element type." << std::endl;
                    exit(1);
                  }
              }


//         // Add sideset names to boundary info
//         mesh.boundary_info->sideset_name(0) = "bottom";
//         mesh.boundary_info->sideset_name(1) = "right";
//         mesh.boundary_info->sideset_name(2) = "top";
//         mesh.boundary_info->sideset_name(3) = "left";
//
//         // Add nodeset names to boundary info
//         mesh.boundary_info->nodeset_name(0) = "bottom";
// 	mesh.boundary_info->nodeset_name(1) = "right";
// 	mesh.boundary_info->nodeset_name(2) = "top";
// 	mesh.boundary_info->nodeset_name(3) = "left";
//
              break;
            }

            //---------------------------------------------------------------------
            // Build a 3D mesh using hexahedral or prismatic elements.
          case 3: {
// 	libmesh_assert_not_equal_to (nx, 0);
// 	libmesh_assert_not_equal_to (ny, 0);
// 	libmesh_assert_not_equal_to (nz, 0);
// 	libmesh_assert_less (xmin, xmax);
// 	libmesh_assert_less (ymin, ymax);
// 	libmesh_assert_less (zmin, zmax);

              assert(nx != 0);
              assert(ny != 0);
              assert(nz != 0);
              assert(xmin < xmax);
              assert(ymin < ymax);
              assert(zmin < zmax);

              // Reserve elements.  Meshes with prismatic elements require
              // twice as many elements.
              switch(type) {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	  case HEX20:
                case HEX27:
// 	  case TET4:  // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
                  {
                    mesh.SetNumberOfElements(nx * ny * nz);
                    ngroup = 1;
                    nbcd   = 6;
                    break;
                  }
//
// 	  case PRISM6:
// 	  case PRISM15:
// 	  case PRISM18:
// 	    {
// 	      mesh.reserve_elem(2*nx*ny*nz);
// 	      break;
// 	    }
//
                default: {
                    std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
                    exit(1);
                  }
              }


              // Reserve nodes.  Quadratic elements need twice as many nodes as linear elements.
              switch(type) {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	  case PRISM6:
// 	    {
// 	      mesh.reserve_nodes( (nx+1)*(ny+1)*(nz+1) );
// 	      break;
// 	    }
//
// 	  case HEX20:
                case HEX27:
// 	  case TET4: // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
// 	  case PRISM15:
// 	  case PRISM18:
                  {
// 	      // FYI: The resulting TET4 mesh will have exactly
// 	      // 5*(nx*ny*nz) + 2*(nx*ny + nx*nz + ny*nz) + (nx+ny+nz) + 1
// 	      // nodes once the additional mid-edge nodes for the HEX27 discretization
// 	      // have been deleted.

                    //mesh.reserve_nodes( (2*nx+1)*(2*ny+1)*(2*nz+1) );
                    mesh.SetNumberOfNodes((2 * nx + 1) * (2 * ny + 1) * (2 * nz + 1));
                    break;
                  }

                default: {
                    std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
                    exit(1);
                  }
              }


              // Build the nodes.
              unsigned int node_id = 0;

              vt[0].resize(mesh.GetNumberOfNodes());
              vt[1].resize(mesh.GetNumberOfNodes());
              vt[2].resize(mesh.GetNumberOfNodes());

              switch(type) {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	  case PRISM6:
// 	    {
// 	      for (unsigned int k=0; k<=nz; k++)
// 		for (unsigned int j=0; j<=ny; j++)
// 		  for (unsigned int i=0; i<=nx; i++)
// 		    {
// 		      if (gauss_lobatto_grid)
// 			{
// 			  mesh.add_point (Point(0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(i)/static_cast<double>(nx))),
// 						0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(j)/static_cast<double>(ny))),
// 						0.5*(1.0 - std::cos(libMesh::pi*static_cast<double>(k)/static_cast<double>(nz)))), node_id++);
// 			}
//
// 		      else
// 			mesh.add_point(Point(static_cast<double>(i)/static_cast<double>(nx),
// 					     static_cast<double>(j)/static_cast<double>(ny),
// 					     static_cast<double>(k)/static_cast<double>(nz)), node_id++);
// 		    }
// 	      break;
// 	    }
//
// 	  case HEX20:
                case HEX27:
// 	  case TET4: // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
// 	  case PRISM15:
// 	  case PRISM18:
                  {
                    for(unsigned int k = 0; k <= (2 * nz); k++)
                      for(unsigned int j = 0; j <= (2 * ny); j++)
                        for(unsigned int i = 0; i <= (2 * nx); i++) {
// 		      if (gauss_lobatto_grid)
// 			{
// 			  // The x,y locations of the point.
// 			  double x=0., y=0., z=0.;
//
// 			  // Shortcut quantities (do not depend on i,j)
// 			  const double a = std::cos( libMesh::pi / static_cast<double>(2*nx) );
// 			  const double b = std::cos( libMesh::pi / static_cast<double>(2*ny) );
//
// 			  // Shortcut quantities (depend on i,j)
// 			  const double c = std::cos( libMesh::pi*i / static_cast<double>(2*nx) );
// 			  const double d = std::cos( libMesh::pi*j / static_cast<double>(2*ny) );
//
// 			  // Additional shortcut quantities (for 3D)
// 			  const double e = std::cos( libMesh::pi / static_cast<double>(2*nz) );
// 			  const double f = std::cos( libMesh::pi*k / static_cast<double>(2*nz) );
//
// 			  // If i is even, compute a normal Gauss-Lobatto point
// 			  if (i%2 == 0)
// 			    x = 0.5*(1.0 - c);
//
// 			  // Otherwise, it is the average of the previous and next points
// 			  else
// 			    x = 0.5*(1.0 - a*c);
//
// 			  // If j is even, compute a normal Gauss-Lobatto point
// 			  if (j%2 == 0)
// 			    y = 0.5*(1.0 - d);
//
// 			  // Otherwise, it is the average of the previous and next points
// 			  else
// 			    y = 0.5*(1.0 - b*d);
//
// 			  // If k is even, compute a normal Gauss-Lobatto point
// 			  if (k%2 == 0)
// 			    z = 0.5*(1.0 - f);
//
// 			  // Otherwise, it is the average of the previous and next points
// 			  else
// 			    z = 0.5*(1.0 - e*f);
//
//
// 			  mesh.add_point (Point(x,y,z), node_id++);
// 			}
//
// 		      else
// 			mesh.add_point(Point(static_cast<double>(i)/static_cast<double>(2*nx),
// 					     static_cast<double>(j)/static_cast<double>(2*ny),
// 					     static_cast<double>(k)/static_cast<double>(2*nz)), node_id++);
                          vt[0][node_id] = (static_cast<double>(i) / static_cast<double>(2 * nx)) * (xmax - xmin) + xmin;
                          vt[1][node_id] = (static_cast<double>(j) / static_cast<double>(2 * ny)) * (ymax - ymin) + ymin;
                          vt[2][node_id] = (static_cast<double>(k) / static_cast<double>(2 * nz)) * (zmax - zmin) + zmin;;

//			  std::cout << "inode: " << node_id <<  " x: " << vt[0][node_id] << "   y: " << vt[1][node_id] << std::endl;

                          node_id++;
                        }

                    break;
                  }


                default: {
                    std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
                    exit(1);
                  }
              }



              // Build the elements.
              unsigned iel = 0;
              mesh.el = new elem(mesh.GetNumberOfElements());
              mesh.el->SetElementGroupNumber(1);
              switch(type) {
// 	  case INVALID_ELEM:
// 	  case HEX8:
// 	    {
// 	      for (unsigned int k=0; k<nz; k++)
// 		for (unsigned int j=0; j<ny; j++)
// 		  for (unsigned int i=0; i<nx; i++)
// 		    {
// 		      Elem* elem = mesh.add_elem(new Hex8);
//
// 		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i,j,k)      );
// 		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
// 		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
// 		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i,j,k+1)    );
// 		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
// 		      elem->set_node(6) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 		      elem->set_node(7) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );
//
// 		      if (k == 0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
//
// 		      if (k == (nz-1))
// 			mesh.boundary_info->add_side(elem, 5, 5);
//
// 		      if (j == 0)
// 			mesh.boundary_info->add_side(elem, 1, 1);
//
// 		      if (j == (ny-1))
// 			mesh.boundary_info->add_side(elem, 3, 3);
//
// 		      if (i == 0)
// 			mesh.boundary_info->add_side(elem, 4, 4);
//
// 		      if (i == (nx-1))
// 			mesh.boundary_info->add_side(elem, 2, 2);
// 		    }
// 	      break;
// 	    }
//
//
//
//
// 	  case PRISM6:
// 	    {
// 	      for (unsigned int k=0; k<nz; k++)
// 		for (unsigned int j=0; j<ny; j++)
// 		  for (unsigned int i=0; i<nx; i++)
// 		    {
// 		      // First Prism
// 		      Elem* elem = NULL;
// 		      elem = mesh.add_elem(new Prism6);
//
// 		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i,j,k)      );
// 		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
// 		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
// 		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i,j,k+1)    );
// 		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
// 		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );
//
// 		      // Add sides for first prism to boundary info object
// 		      if (i==0)
// 			mesh.boundary_info->add_side(elem, 3, 4);
//
// 		      if (j==0)
// 			mesh.boundary_info->add_side(elem, 1, 1);
//
// 		      if (k==0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
//
// 		      if (k == (nz-1))
// 			mesh.boundary_info->add_side(elem, 4, 5);
//
// 		      // Second Prism
// 		      elem = mesh.add_elem(new Prism6);
//
// 		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
// 		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
// 		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
// 		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );
//
// 		      // Add sides for second prism to boundary info object
// 		      if (i == (nx-1))
// 			mesh.boundary_info->add_side(elem, 1, 2);
//
// 		      if (j == (ny-1))
// 			mesh.boundary_info->add_side(elem, 2, 3);
//
// 		      if (k==0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
//
// 		      if (k == (nz-1))
// 			mesh.boundary_info->add_side(elem, 4, 5);
// 		    }
// 	      break;
// 	    }
//
//
//
//
//
//
// 	  case HEX20:
                case HEX27:
// 	  case TET4: // TET4's are created from an initial HEX27 discretization
// 	  case TET10: // TET10's are created from an initial HEX27 discretization
// 	  case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
                  {
                    double LocalToGlobalNodePerElement[27];

                    for(unsigned int k = 0; k < (2 * nz); k += 2)
                      for(unsigned int j = 0; j < (2 * ny); j += 2)
                        for(unsigned int i = 0; i < (2 * nx); i += 2) {
// 		      Elem* elem = (type == HEX20) ?
// 			mesh.add_elem(new Hex20) :
// 			mesh.add_elem(new Hex27);

                          mesh.el->SetElementGroup(iel, 1);
                          mesh.el->SetElementMaterial(iel, 2);
                          type_elem_flag[0] = true; // hex
                          type_elem_flag[3] = true; // quad face
                          mesh.el->AddToElementNumber(1, "Hex");
                          mesh.el->SetElementType(iel, 0);

                          LocalToGlobalNodePerElement[0] = idx(type, nx, ny, i,  j,  k) + 1;
                          LocalToGlobalNodePerElement[1] = idx(type, nx, ny, i + 2, j,  k) + 1;
                          LocalToGlobalNodePerElement[2] = idx(type, nx, ny, i + 2, j + 2, k) + 1;
                          LocalToGlobalNodePerElement[3] = idx(type, nx, ny, i,  j + 2, k) + 1;
                          LocalToGlobalNodePerElement[4] = idx(type, nx, ny, i,  j,  k + 2) + 1;
                          LocalToGlobalNodePerElement[5] = idx(type, nx, ny, i + 2, j,  k + 2) + 1;
                          LocalToGlobalNodePerElement[6] = idx(type, nx, ny, i + 2, j + 2, k + 2) + 1;
                          LocalToGlobalNodePerElement[7] = idx(type, nx, ny, i,  j + 2, k + 2) + 1;
                          LocalToGlobalNodePerElement[8] = idx(type, nx, ny, i + 1, j,  k) + 1;
                          LocalToGlobalNodePerElement[9] = idx(type, nx, ny, i + 2, j + 1, k) + 1;
                          LocalToGlobalNodePerElement[10] = idx(type, nx, ny, i + 1, j + 2, k) + 1;
                          LocalToGlobalNodePerElement[11] = idx(type, nx, ny, i,  j + 1, k) + 1;

                          // mid-point - up
                          LocalToGlobalNodePerElement[12] = idx(type, nx, ny, i + 1, j,  k + 2) + 1;
                          LocalToGlobalNodePerElement[13] = idx(type, nx, ny, i + 2, j + 1, k + 2) + 1;
                          LocalToGlobalNodePerElement[14] = idx(type, nx, ny, i + 1, j + 2, k + 2) + 1;
                          LocalToGlobalNodePerElement[15] = idx(type, nx, ny, i,  j + 1, k + 2) + 1 ;

                          // vertices - middle
                          LocalToGlobalNodePerElement[16] = idx(type, nx, ny, i,  j,  k + 1) + 1;
                          LocalToGlobalNodePerElement[17] = idx(type, nx, ny, i + 2, j,  k + 1) + 1;
                          LocalToGlobalNodePerElement[18] = idx(type, nx, ny, i + 2, j + 2, k + 1) + 1;
                          LocalToGlobalNodePerElement[19] = idx(type, nx, ny, i,  j + 2, k + 1) + 1;

                          // mid-point - middle
                          LocalToGlobalNodePerElement[20] = idx(type, nx, ny, i + 1, j,  k + 1) + 1;
                          LocalToGlobalNodePerElement[21] = idx(type, nx, ny, i + 2, j + 1, k + 1) + 1;
                          LocalToGlobalNodePerElement[22] = idx(type, nx, ny, i + 1, j + 2, k + 1) + 1;
                          LocalToGlobalNodePerElement[23] = idx(type, nx, ny, i,  j + 1, k + 1) + 1;

                          // center - bottom
                          LocalToGlobalNodePerElement[24] = idx(type, nx, ny, i + 1, j + 1, k) + 1;

                          // center - top
                          LocalToGlobalNodePerElement[25] = idx(type, nx, ny, i + 1, j + 1, k + 2) + 1;

                          // center - middle
                          LocalToGlobalNodePerElement[26] = idx(type, nx, ny, i + 1, j + 1, k + 1) + 1;

                          // connectivity
                          for(unsigned iloc = 0; iloc < 27; iloc++) {
                            mesh.el->SetElementDofIndex(iel, iloc, LocalToGlobalNodePerElement[iloc] - 1u);
                          }

                          if(k == 0) {
                            mesh.el->SetFaceElementIndex(iel, 4, -2);
                            mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(0, "bottom"));
                          }


                          if(k == 2 * (nz - 1)) {
                            mesh.el->SetFaceElementIndex(iel, 5, -7);
                            mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(5, "top"));
                          }


                          if(j == 0) {
                            mesh.el->SetFaceElementIndex(iel, 0, -3);
                            mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(1, "front"));
                          }


                          if(j == 2 * (ny - 1)) {
                            mesh.el->SetFaceElementIndex(iel, 2, -5);
                            mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(3, "behind"));
                          }


                          if(i == 0) {
                            mesh.el->SetFaceElementIndex(iel, 3, -6);
                            mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(4, "left"));
                          }


                          if(i == 2 * (nx - 1)) {
                            mesh.el->SetFaceElementIndex(iel, 1, -4);
                            mesh._boundaryinfo.insert(std::pair<unsigned int, std::string>(2, "right"));
                          }


                          iel++;

// 		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k)  );
// 		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k)  );
// 		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k)  );
// 		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k)  );

// 		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+2));
// 		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+2));
// 		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+2));
// 		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+2));
// 		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k)  );

// 		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k)  );
// 		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k)  );
// 		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k)  );
// 		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+1));

// 		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+1));
// 		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+1));
// 		      elem->set_node(15) = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+1));
// 		      elem->set_node(16) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+2));
// 		      elem->set_node(17) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+2));
// 		      elem->set_node(18) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+2));
// 		      elem->set_node(19) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+2));
// 		      if ((type == HEX27) || (type == TET4) || (type == TET10) || (type == PYRAMID5))
// 			{
// 			  elem->set_node(20) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 			  elem->set_node(21) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+1));
// 			  elem->set_node(22) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+1));
// 			  elem->set_node(23) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+1));
// 			  elem->set_node(24) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+1));
// 			  elem->set_node(25) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
// 			  elem->set_node(26) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 			}
//
//
// 		      if (k == 0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
//
// 		      if (k == 2*(nz-1))
// 			mesh.boundary_info->add_side(elem, 5, 5);
//
// 		      if (j == 0)
// 			mesh.boundary_info->add_side(elem, 1, 1);
//
// 		      if (j == 2*(ny-1))
// 			mesh.boundary_info->add_side(elem, 3, 3);
//
// 		      if (i == 0)
// 			mesh.boundary_info->add_side(elem, 4, 4);
//
// 		      if (i == 2*(nx-1))
// 			mesh.boundary_info->add_side(elem, 2, 2);
                        }
                    break;
                  }
//
//
//
//
// 	  case PRISM15:
// 	  case PRISM18:
// 	    {
// 	      for (unsigned int k=0; k<(2*nz); k += 2)
// 		for (unsigned int j=0; j<(2*ny); j += 2)
// 		  for (unsigned int i=0; i<(2*nx); i += 2)
// 		    {
// 		      // First Prism
// 		      Elem* elem = NULL;
// 		      elem = ((type == PRISM15) ?
// 			      mesh.add_elem(new Prism15) :
// 			      mesh.add_elem(new Prism18));
//
// 		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k)  );
// 		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k)  );
// 		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k)  );
// 		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+2));
// 		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+2));
// 		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+2));
// 		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k)  );
// 		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k)  );
// 		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+1));
// 		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+1));
// 		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+1));
// 		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+2));
// 		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
// 		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+2));
// 		      if (type == PRISM18)
// 			{
// 			  elem->set_node(15) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+1));
// 			  elem->set_node(16) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 			  elem->set_node(17) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+1));
// 			}
//
// 		      // Add sides for first prism to boundary info object
// 		      if (i==0)
// 			mesh.boundary_info->add_side(elem, 3, 4);
//
// 		      if (j==0)
// 			mesh.boundary_info->add_side(elem, 1, 1);
//
// 		      if (k==0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
//
// 		      if (k == 2*(nz-1))
// 			mesh.boundary_info->add_side(elem, 4, 5);
//
//
// 		      // Second Prism
// 		      elem = ((type == PRISM15) ?
// 			      mesh.add_elem(new Prism15) :
// 			      mesh.add_elem(new Prism18));
//
// 		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k)     );
// 		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k)   );
// 		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i,j+2,k)     );
// 		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k+2)   );
// 		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+2) );
// 		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i,j+2,k+2)   );
// 		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k)  );
// 		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k)  );
// 		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
// 		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k+1)  );
// 		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+1));
// 		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,j+2,k+1)  );
// 		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+2));
// 		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+2));
// 		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
// 		      if (type == PRISM18)
// 			{
// 			  elem->set_node(15)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+1));
// 			  elem->set_node(16)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+1));
// 			  elem->set_node(17)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
// 			}
//
// 		      // Add sides for second prism to boundary info object
// 		      if (i == 2*(nx-1))
// 			mesh.boundary_info->add_side(elem, 1, 2);
//
// 		      if (j == 2*(ny-1))
// 			mesh.boundary_info->add_side(elem, 2, 3);
//
// 		      if (k==0)
// 			mesh.boundary_info->add_side(elem, 0, 0);
//
// 		      if (k == 2*(nz-1))
// 			mesh.boundary_info->add_side(elem, 4, 5);
//
//  		    }
//  	      break;
//  	    }
//

                default: {
                    std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
                    exit(1);
                  }
              }

// 	//.......................................
// 	// Scale the nodal positions
// 	for (unsigned int p=0; p<mesh.n_nodes(); p++)
// 	  {
// 	    mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
// 	    mesh.node(p)(1) = (mesh.node(p)(1))*(ymax-ymin) + ymin;
// 	    mesh.node(p)(2) = (mesh.node(p)(2))*(zmax-zmin) + zmin;
// 	  }
//
//
//
//
// 	// Additional work for tets and pyramids: we take the existing
// 	// HEX27 discretization and split each element into 24
// 	// sub-tets or 6 sub-pyramids.
// 	//
// 	// 24 isn't the minimum-possible number of tets, but it
// 	// obviates any concerns about the edge orientations between
// 	// the various elements.
// 	if ((type == TET4) ||
// 	    (type == TET10) ||
// 	    (type == PYRAMID5))
// 	  {
// 	    // Temporary storage for new elements. (24 tets per hex, 6 pyramids)
// 	    std::vector<Elem*> new_elements;
//
// 	    if ((type == TET4) || (type == TET10))
// 	      new_elements.reserve(24*mesh.n_elem());
// 	    else
// 	      new_elements.reserve(6*mesh.n_elem());
//
// 	    // Create tetrahedra or pyramids
// 	    {
// 	      MeshBase::element_iterator       el     = mesh.elements_begin();
// 	      const MeshBase::element_iterator end_el = mesh.elements_end();
//
// 	      for ( ; el != end_el;  ++el)
// 		{
// 		  // Get a pointer to the HEX27 element.
// 		  Elem* base_hex = *el;
//
// 		  // Get a pointer to the node located at the HEX27 centroid
// 		  Node* apex_node = base_hex->get_node(26);
//
// 		  for (unsigned int s=0; s<base_hex->n_sides(); ++s)
// 		    {
// 		      // Get the boundary ID for this side
// 		      boundary_id_type b_id = mesh.boundary_info->boundary_id(*el, s);
//
// 		      // Need to build the full-ordered side!
// 		      AutoPtr<Elem> side = base_hex->build_side(s);
//
// 		      if ((type == TET4) || (type == TET10))
// 			{
// 			  // Build 4 sub-tets per side
// 			  for (unsigned int sub_tet=0; sub_tet<4; ++sub_tet)
// 			    {
// 			      new_elements.push_back( new Tet4 );
// 			      Elem* sub_elem = new_elements.back();
// 			      sub_elem->set_node(0) = side->get_node(sub_tet);
// 			      sub_elem->set_node(1) = side->get_node(8);                           // centroid of the face
// 			      sub_elem->set_node(2) = side->get_node(sub_tet==3 ? 0 : sub_tet+1 ); // wrap-around
// 			      sub_elem->set_node(3) = apex_node;                                   // apex node always used!
//
// 			      // If the original hex was a boundary hex, add the new sub_tet's side
// 			      // 0 with the same b_id.  Note: the tets are all aligned so that their
// 			      // side 0 is on the boundary.
// 			      if (b_id != BoundaryInfo::invalid_id)
// 				mesh.boundary_info->add_side(sub_elem, 0, b_id);
// 			    }
// 			} // end if ((type == TET4) || (type == TET10))
//
// 		      else // type==PYRAMID5
// 			{
// 			  // Build 1 sub-pyramid per side.
// 			  new_elements.push_back(new Pyramid5);
// 			  Elem* sub_elem = new_elements.back();
//
// 			  // Set the base.  Note that since the apex is *inside* the base_hex,
// 			  // and the pyramid uses a counter-clockwise base numbering, we need to
// 			  // reverse the [1] and [3] node indices.
// 			  sub_elem->set_node(0) = side->get_node(0);
// 			  sub_elem->set_node(1) = side->get_node(3);
// 			  sub_elem->set_node(2) = side->get_node(2);
// 			  sub_elem->set_node(3) = side->get_node(1);
//
// 			  // Set the apex
// 			  sub_elem->set_node(4) = apex_node;
//
// 			  // If the original hex was a boundary hex, add the new sub_pyr's side
// 			  // 4 (the square base) with the same b_id.
// 			  if (b_id != BoundaryInfo::invalid_id)
// 			    mesh.boundary_info->add_side(sub_elem, 4, b_id);
// 			} // end else type==PYRAMID5
// 		    }
// 		}
// 	    }
//
//
// 	    // Delete the original HEX27 elements from the mesh, and the boundary info structure.
// 	    {
// 	      MeshBase::element_iterator       el     = mesh.elements_begin();
// 	      const MeshBase::element_iterator end_el = mesh.elements_end();
//
// 	      for ( ; el != end_el;  ++el)
// 		{
// 		  mesh.boundary_info->remove(*el); // Safe even if *el has no boundary info.
// 		  mesh.delete_elem(*el);
// 		}
// 	    }
//
// 	    // Add the new elements
// 	    for (unsigned int i=0; i<new_elements.size(); ++i)
// 	      mesh.add_elem(new_elements[i]);
//
// 	  } // end if (type == TET4,TET10,PYRAMID5
//
//
// 	// Use all_second_order to convert the TET4's to TET10's
// 	if (type == TET10)
// 	  {
// 	    mesh.all_second_order();
// 	  }
//
//         // Add sideset names to boundary info (Z axis out of the screen)
//         mesh.boundary_info->sideset_name(0) = "back";
//         mesh.boundary_info->sideset_name(1) = "bottom";
//         mesh.boundary_info->sideset_name(2) = "right";
//         mesh.boundary_info->sideset_name(3) = "top";
//         mesh.boundary_info->sideset_name(4) = "left";
//         mesh.boundary_info->sideset_name(5) = "front";
//
//         // Add nodeset names to boundary info
//         mesh.boundary_info->nodeset_name(0) = "back";
// 	mesh.boundary_info->nodeset_name(1) = "bottom";
// 	mesh.boundary_info->nodeset_name(2) = "right";
// 	mesh.boundary_info->nodeset_name(3) = "top";
// 	mesh.boundary_info->nodeset_name(4) = "left";
// 	mesh.boundary_info->nodeset_name(5) = "front";

              break;
            } // end case dim==3

          default: {
              std::cout << " Error! " << std::endl;
            }
        }

      }


    }

  }

}
