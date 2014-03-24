// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// This file was massively overhauled
#include "gambit_io.hpp"


//C++ includes
#include <fstream>
#include <set>
#include <cstring>

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"





// anonymous namespace to hold local data
namespace
{

/**
 * Defines a structure to hold boundary element information.
 *
 * We use a set because it keeps the nodes unique and ordered, and can be
 * easily compared to another set of nodes (the ones on the element side)
 */
struct boundaryElementInfo {
  std::set<unsigned int> nodes;
  unsigned int id;
};

/**
 * Defines mapping from libMesh element types to Gmsh element types.
 */
struct elementDefinition {
  std::string label;
  std::vector<unsigned int> nodes;
  ElemType type;
  unsigned int exptype;
  unsigned int dim;
  unsigned int nnodes;
};


// maps from a libMesh element type to the proper
// Gmsh elementDefinition.  Placing the data structure
// here in this anonymous namespace gives us the
// benefits of a global variable without the nasty
// side-effects
std::map<ElemType, elementDefinition> eletypes_exp;
std::map<unsigned int, elementDefinition> eletypes_imp;



// ------------------------------------------------------------
// helper function to initialize the eletypes map
void init_eletypes()
{
  if ( eletypes_exp.empty() && eletypes_imp.empty() ) {
    // This should happen only once.  The first time this method
    // is called the eletypes data struture will be empty, and
    // we will fill it.  Any subsequent calls will find an initialized
    // eletypes map and will do nothing.

    //==============================
    // setup the element definitions
    elementDefinition eledef;

    // use "swap trick" from Scott Meyer's "Effective STL" to initialize
    // eledef.nodes vector

    // POINT (only Gmsh)
    {
      eledef.exptype = 15;
      eledef.dim     = 0;
      eledef.nnodes  = 1;
      eledef.nodes.clear();

      // import only
      eletypes_imp[15] = eledef;
    }

    // EDGE2
    {
      eledef.type    = EDGE2;
      eledef.dim     = 1;
      eledef.nnodes  = 2;
      eledef.exptype = 1;
      eledef.nodes.clear();

      eletypes_exp[EDGE2] = eledef;
      eletypes_imp[1]     = eledef;
    }

    // EDGE3
    {
      eledef.type    = EDGE3;
      eledef.dim     = 1;
      eledef.nnodes  = 3;
      eledef.exptype = 8;
      eledef.nodes.clear();

      eletypes_exp[EDGE3] = eledef;
      eletypes_imp[8]     = eledef;
    }

    // TRI3 Gambit
    {
      eledef.type    = TRI3;
      eledef.dim     = 2;
      eledef.nnodes  = 3;
      eledef.exptype = 2;
      const unsigned int nodes[] = {0,1,2};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[TRI3] = eledef;
      eletypes_imp[2] = eledef;
    }

    // TRI6 Gambit
    {
      eledef.type    = TRI6;
      eledef.dim     = 2;
      eledef.nnodes  = 6;
      eledef.exptype = 9;
      const unsigned int nodes[] = {0,3,1,4,2,5};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[TRI6] = eledef;
      eletypes_imp[9]    = eledef;
    }

    // QUAD4 Gambit
    {
      eledef.type    = QUAD4;
      eledef.dim     = 2;
      eledef.nnodes  = 4;
      eledef.exptype = 3;
      const unsigned int nodes[] = {0,1,2,3};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[QUAD4] = eledef;
      eletypes_imp[3]     = eledef;
    }

    // QUAD8 Gambit
    // TODO: what should be done with this on writing?
    {
      eledef.type    = QUAD8;
      eledef.dim     = 2;
      eledef.nnodes  = 8;
      eledef.exptype = 100;
      const unsigned int nodes[] = {0,4,1,5,2,6,3,7};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[QUAD8] = eledef;
      eletypes_imp[19]    = eledef;
    }

    // QUAD9 Gambit
    {
      eledef.type    = QUAD9;
      eledef.dim     = 2;
      eledef.nnodes  = 9;
      eledef.exptype = 10;
       const unsigned int nodes[] =  {0,4,1,5,2,6,3,7,8};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[QUAD9] = eledef;
      eletypes_imp[10]    = eledef;
    }

    // HEX8 Gambit
    {
      eledef.type    = HEX8;
      eledef.dim     = 3;
      eledef.nnodes  = 8;
      eledef.exptype = 5;
      // eledef.nodes.clear();
      const unsigned int nodes[] = {4,5,0,1,7,6,3,2};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[HEX8] = eledef;
      eletypes_imp[5]    = eledef;
    }

    // HEX20 Gambit
    // TODO: what should be done with this on writing?
    {
      eledef.type    = HEX20;
      eledef.dim     = 3;
      eledef.nnodes  = 20;
      eledef.exptype = 101;
      //	const unsigned int nodes[] = {0,8,1,11,9,3,10,2,12,13,15,14,4,16,5,19,17,7,18,6};

      const unsigned int nodes[] = {4,16,5,12,13,0,8,1,19,17,11,9,7,18,6,15,14,3,10,2};
      // !! negative jacobian ???????????????????
//       const unsigned int nnodes = sizeof ( nodes ) /sizeof ( nodes[0] );
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );

      eletypes_exp[HEX20] = eledef;
      eletypes_imp[20]    = eledef;
    }

    // HEX27 Gambit
    {
      eledef.type    = HEX27;
      eledef.dim     = 3;
      eledef.nnodes  = 27;
      eledef.exptype = 12;
    //const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,15,16,19,17,18,20,21,24,22,23,25,26};
      //const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
      //const unsigned int nodes[] = {4,16,5,12,21,13,0,8,1,19,25,17,24,26,22,11,20,9,7,18,6,15,23,14,3,10,2};
      //const unsigned int nodes[] = {5,17,6,13,22,14,1,9,2,16,25,18,21,26,23,8,20,10,4,19,7,12,24,15,0,11,3};
/*      const unsigned int nodes[] = {0,2,8,6,18,20,26,24,1,5,7,3,9,11,17,15,19,23,25,21,4,10,14,16,12,22,13};*/


 const unsigned int nodes[] = {0,8,1,11,20,9,3,10,2,12,21,13,24,26,22,15,23,14,4,16,5,19,25,17,7,18,6};
 
//         const unsigned int nodes[] = {4,16,5,19,25,17,7,18,6,12,21,13,24,26,22,15,23,14,0,8,1,11,20,9,3,10,2};
    //  const unsigned int nodes[] = {1,8,0,9,20,11,2,10,3,13,21,12,22,26,24,14,23,15,5,16,4,17,25,19,6,18,7};
//        const unsigned int nodes[] =  {0,2,8,6,18,20,26,24,1,5,7,3,9,11,17,15,19,23,25,21,4,10,14,16,12,22,13};
       std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );

      eletypes_exp[HEX27] = eledef;
      eletypes_imp[12]    = eledef;
    }

    // TET4 Gambit
    {
      eledef.type    = TET4;
      eledef.dim     = 3;
      eledef.nnodes  = 4;
      eledef.exptype = 4;
      //eledef.nodes.clear();
      const unsigned int nodes[] = {0,1,2,3};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[TET4] = eledef;
      eletypes_imp[4]    = eledef;
    }

    // TET10 Gambit
    {
      eledef.type    = TET10;
      eledef.dim     = 3;
      eledef.nnodes  = 10;
      eledef.exptype = 11;
      const unsigned int nodes[] = {0,6,2,7,9,3,4,5,8,1};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[TET10] = eledef;
      eletypes_imp[11]    = eledef;
    }

    // PRISM6 Gambit
    {
      eledef.type    = PRISM6;
      eledef.dim     = 3;
      eledef.nnodes  = 6;
      eledef.exptype = 6;
      //eledef.nodes.clear();
      const unsigned int nodes[] = {2,0,1,5,3,4};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[PRISM6] = eledef;
      eletypes_imp[6]      = eledef;
    }

    // PRISM15 Gambit
    // TODO: what should be done with this on writing?
    {
      eledef.type    = PRISM15;
      eledef.dim     = 3;
      eledef.nnodes  = 15;
      eledef.exptype = 103;
      // eledef.nodes.clear();
      const unsigned int nodes[] = {4,12,3,13,14,5,10,9,11,1,6,0,7,8,2};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      eletypes_exp[PRISM15] = eledef;
      eletypes_imp[16] = eledef;
    }

    // PRISM18 Gambit
    {
      eledef.type    = PRISM18;
      eledef.dim     = 3;
      eledef.nnodes  = 18;
      eledef.exptype = 13;
      //   const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,
      //                            12,14,13,15,17,16};
      const unsigned int nodes[] = {5,13,4,14,12,3,11,16,10,17,15,9,2,7,1,8,6,0};
      std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );

      eletypes_exp[PRISM18] = eledef;
      eletypes_imp[13]      = eledef;
    }

    // PYRAMID5
    {
      eledef.type    = PYRAMID5;
      eledef.dim     = 3;
      eledef.nnodes  = 5;
      eledef.exptype = 7;
     // eledef.nodes.clear();
      
     const unsigned int nodes[] = {0,1,2,3,4};
     std::vector<unsigned int> ( nodes, nodes+eledef.nnodes ).swap ( eledef.nodes );
      
      eletypes_exp[PYRAMID5] = eledef;
      eletypes_imp[7]        = eledef;
    }

    //==============================
  }
}
} // end anonymous namespace

// ------------------------------------------------------------
// ------------------------------------------------------------


// GambitIO  members
void GambitIO::read ( const std::string& name )
{
  std::ifstream in ( name.c_str() );
  this->read_mesh ( in );
}



// -----------------------------------------
void GambitIO::read_mesh ( std::istream& in )
{

  libmesh_assert ( in.good() );
  // initialize the map with element types
  init_eletypes();

  // clear any data in the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
//  MeshBase& mesh = GambitIO::mesh();
  mesh.clear();

  // some variables
  const int  bufLen = 256;
  char  buf[bufLen+1];


  // reserve blocks    ************************
  while (strncmp ( buf,"NDFCD",5 ) != 0) in >> buf;

  in >>buf;

  unsigned int numElem = 0;
  unsigned int numNodes = 0;
  
  std::cout << " reading from standard mesh file (gambit format) "
            << "\n" << "  numNodes    " << "numElem      " << std::endl; //<< "ngroup   "
           // <<  "nbdary "  <<  "NDFCD  " <<  "NDFVL  " << std::endl;
  in >> numNodes >>  numElem;
  std::cout <<   "     " << numNodes << "           " <<  numElem << std::endl;

  in >>buf >>buf >>buf;
  
  std::cout << "The mesh dimension is " << buf << std::endl;
  
  int dimension = atoi(&buf[0]);

  std::cout << "********* Remember to check the consistency of your application dimension ********************** " << std::endl;
  std::cout << "********* MeshBase dimension is 3 by default ********************** " << mesh.spatial_dimension() << std::endl;

  mesh.reserve_nodes ( numNodes );
  mesh.reserve_elem ( numElem );
 

  // coordinates *****************************
  // map to hold the node numbers for translation
  // note that the nodes can be non-consecutive
  // Gambit counts the nodes starting from 1 instead of 0
  std::map<unsigned int, unsigned int> nodetrans;

  
  while ( strncmp ( buf,"COORDINATES",11 ) != 0 ) in >> buf;

  in >> buf;
  // read in the nodal coordinates and form points.
  Real x, y;
  unsigned int id;
  // add the nodal coordinates to the mesh
  for ( unsigned int i=0; i<numNodes; ++i ) {

    if ( dimension == 2 ) {
     in >> id >> x >> y;
     mesh.add_point(Point(x, y),i);
    }
    else if ( dimension == 3 ) {
     Real z;
     in >> id >> x >> y >> z;
     mesh.add_point(Point ( x, y, z ),i );
    }
     nodetrans[id] = i;
  }
  // read the $ENDNOD delimiter
  in >> buf;

  // elements *****************************
//In the .neu (.gam) file there are two parameters
  //NDFCD = number of coordinate directions (it's the dimension!)
  //NDFVL = number of velocity components
  while ( strncmp ( buf,"ELEMENTS/CELLS",14 ) != 0 ) in >> buf;
  in >> buf;
  
  // read the elements
      //the scheme is:
    //for every element, read
    // id face type CONNECTIVITY
  //face seems to be meaningless...
  double ElNumNodes;
  for ( unsigned int iel=0; iel<numElem; ++iel ) {
    unsigned int id, type, face; //elementary=0, partition = 1, numTags;
    in >> id >> face >>  type;

    // number of nodes  
    /*numNodes*/ElNumNodes=type;
    int itype=0;
    // consult the import element table which element to build

    if ( dimension == 2 ) {
    switch ( type ) {
      case 3: itype=2; break;  // tri 3
      case 6: itype=9; break;  // tri 6
      case 4: itype=3; break;  // quad 4
      case 8: itype=19; break; // quad 8
      case 9: itype=10; break; // quad 9
      default:
      std::cout << "error element" << std::endl;
      abort();
      break;
     } 
    }
    else {

      switch ( type ) {
    case 4:
      itype=4;
      break; // tet4
    case 10:
      itype=11;
      break; // tet10
    case 6:
      itype=6;
      break; // wedge6
    case 15:
      itype=16;
      break; // wedge15
    case 18:
      itype=13;
      break; // wedge18
    case 8:
      itype=5;
      break; // hex 8
    case 20:
      itype=20;
      break; // hex 20
    case 27:
      itype=12;
      break; // hex 27
    case 5:
      itype=13;
      break; // pyramid 5
    default:
      std::cout << "error element" << std::endl;
      abort();
    break;
    }
   }
   
    const elementDefinition& eletype = eletypes_imp[itype];
    // add the elements to the mesh
    Elem* elem = mesh.add_elem ( Elem::build ( eletype.type ).release() );

    // add node pointers to the elements
    //read connectivity

if ( dimension == 2 ) {
// ***********************************************************************
    int nod;
     for (unsigned int i=0; i < ElNumNodes; i++) {
       in>>nod;
       elem->set_node(eletype.nodes[i]) = mesh.node_ptr(nodetrans[nod]);
        std::cout << mesh.n_nodes() <<  "  " << eletype.nodes[i] << "  " << nodetrans[nod] << "  " << nod << "  " <<  "  " << std::endl;
     }
}
// ***********************************************************
else {
  // ***********************************************************
   
    int vecnode[27];
    for ( unsigned int i=0; i < ElNumNodes; i++ )   in >> vecnode[i];
    for ( unsigned int i=0; i < ElNumNodes; i++ ) {
      elem->set_node ( eletype.nodes[i] ) = mesh.node_ptr (vecnode[i]-1);
    }
    // ************************************************************************
}
  }
  // read the $ENDELM delimiter
  in >> buf;

}

