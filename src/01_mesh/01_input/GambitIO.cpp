/*=========================================================================

 Program: FEMUS
 Module: GambitIO
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//local include
#include "GambitIO.hpp"
#include "Mesh.hpp"

//C++ include
#include "cstdio"
#include "fstream"


namespace femus {

  const unsigned GambitIO::GambitToFemusVertexIndex[N_GEOM_ELS][MAX_EL_N_NODES] = {
    {
      4, 16, 0, 15, 23, 11, 7, 19, 3,
      12, 20, 8, 25, 26, 24, 14, 22, 10,
      5, 17, 1, 13, 21, 9, 6, 18, 2
    },
    {
      0, 4, 1, 6, 5,
      2, 7, 8, 9, 3
    },
    {
      3, 11, 5, 9, 10, 4,
      12, 17, 14, 15, 16, 13,
      0, 8, 2, 6, 7, 1
    },
    {0, 4, 1, 5, 2, 6, 3, 7, 8},
    {0, 3, 1, 4, 2, 5, 6},
    {0, 2, 1}
  };


  const unsigned GambitIO::GambitToFemusFaceIndex[N_GEOM_ELS][MAX_EL_N_FACES] = {
    {0, 4, 2, 5, 3, 1},
    {0, 1, 2, 3},
    {2, 1, 0, 4, 3},
    {0, 1, 2, 3},
    {0, 1, 2},
    {0, 1}
  };

 
//   const unsigned GambitIO::_numberOfMissedBiquadraticNodes[N_GEOM_ELS] = {0,5,3,0,1,0};
//   
//   const double GambitIO::_baricentricWeight[N_GEOM_ELS][5][18] = {
//     {},
//     { 
//       { -1./9., -1./ 9., -1./9.,  0    , 4./9., 4./9., 4./9., 0.   , 0.   , 0.   },
//       { -1./9., -1./ 9.,  0.   , -1./9., 4./9., 0.   , 0.   , 4./9., 4./9., 0.   },
//       {  0.   , -1./ 9., -1./9., -1./9., 0.   , 4./9., 0.   , 0.   , 4./9., 4./9.},
//       { -1./9.,  0.    , -1./9., -1./9., 0.   , 0.   , 4./9., 4./9., 0.   , 4./9.},
//       { -1./8., -1./ 8., -1./8., -1./8., 1./4., 1./4., 1./4., 1./4., 1./4., 1./4.}
//     },
//     {
//       { -1./9., -1./9., -1./9., 0.   ,  0.   ,  0.   , 4./9., 4./9., 4./9.},
//       {  0.   ,  0.   ,  0.   ,-1./9., -1./9., -1./9., 0.   , 0.   , 0.   , 4./9., 4./9., 4./9.},
//       {  0.   ,  0.   ,  0.   , 0.   ,  0.   ,  0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
//         -1./9., -1./9., -1./9., 4./9.,  4./9.,  4./9.},
//     },
//     {},
//     {{ -1. / 9., -1. / 9., -1. / 9., 4. / 9., 4. / 9., 4. / 9.}}
//   };

  void GambitIO::read(const std::string& name, vector < vector < double> > &coords, const double Lref, std::vector<bool> &type_elem_flag, const bool read_groups, const bool read_boundary_groups) {

    Mesh& mesh = GetMesh();

    std::ifstream inf;
    std::string str2;
    unsigned ngroup;
    unsigned nbcd;
    unsigned dim, dimNodes;
    double x, y, z;
    unsigned nvt;
    unsigned nvt0;
    unsigned nel;

    mesh.SetLevel(0);
    // read control data ******************** A
    inf.open(name.c_str());
    if(!inf) {
      std::cout << "Generic-mesh file " << name << " can not read parameters\n";
      exit(0);
    }
    str2 = "0";
    while(str2.compare("NDFVL") != 0) inf >> str2;
    inf >> nvt >> nel >>  ngroup >> nbcd >> dim >> dimNodes ;
    nvt0 = nvt;
    mesh.SetDimension(dim);
    mesh.SetNumberOfElements(nel);
    mesh.SetNumberOfNodes(nvt);
    inf >> str2;
    if(str2.compare("ENDOFSECTION") != 0) {
      std::cout << "error control data mesh" << std::endl;
      exit(0);
    }
//   std::cout << "***************" << _dimension << std::endl;
    inf.close();
    // end read control data **************** A
    // read ELEMENT/cell ******************** B
    inf.open(name.c_str());
    if(!inf) {
      std::cout << "Generic-mesh file " << name << " cannot read elements\n";
      exit(0);
    }
    mesh.el = new elem(nel);
    while(str2.compare("ELEMENTS/CELLS") != 0) inf >> str2;
    inf >> str2;
    for(unsigned iel = 0; iel < nel; iel++) {
      mesh.el->SetElementGroup(iel, 1);
      unsigned nve;
      unsigned elementType;
      inf >> str2 >> str2 >> nve;
      if(nve == 27) {
        type_elem_flag[0] = type_elem_flag[3] = true;
        mesh.el->AddToElementNumber(1, "Hex");
        mesh.el->SetElementType(iel, 0);
	elementType = 0;
      }
      else if(nve == 10) {
        type_elem_flag[1] = type_elem_flag[4] = true;
        mesh.el->AddToElementNumber(1, "Tet");
        mesh.el->SetElementType(iel, 1);
	elementType = 1;
      }
      else if(nve == 18) {
        type_elem_flag[2] = type_elem_flag[3] = type_elem_flag[4] = true;
        mesh.el->AddToElementNumber(1, "Wedge");
        mesh.el->SetElementType(iel, 2);
	elementType = 2;
      }
      else if(nve == 9) {
        type_elem_flag[3] = true;
        mesh.el->AddToElementNumber(1, "Quad");
        mesh.el->SetElementType(iel, 3);
	elementType = 3;
      }
      else if(nve == 6 && mesh.GetDimension() == 2) {
        type_elem_flag[4] = true;
        mesh.el->AddToElementNumber(1, "Triangle");
        mesh.el->SetElementType(iel, 4);
	elementType = 4;
      }
      else if(nve == 3 && mesh.GetDimension() == 1) {
        mesh.el->AddToElementNumber(1, "Line");
        mesh.el->SetElementType(iel, 5);
	elementType = 5;
      }
      else {
        std::cout << "Error! Invalid element type in reading Gambit File!" << std::endl;
        std::cout << "Error! Use a second order discretization" << std::endl;
        exit(0);
      }
      for(unsigned i = 0; i < nve; i++) {
        unsigned inode = GambitIO::GambitToFemusVertexIndex[mesh.el->GetElementType(iel)][i];
        unsigned value;
        inf >> value;
        mesh.el->SetElementDofIndex(iel, inode, value - 1u);
      }
    }
    
    //GambitIO::BiquadraticNodesNotInGambit(mesh);
    //nvt = mesh.GetNumberOfNodes();
       
    inf >> str2;
    if(str2.compare("ENDOFSECTION") != 0) {
      std::cout << "error element data mesh" << std::endl;
      exit(0);
    }
    inf.close();

    // end read  ELEMENT/CELL **************** B

    // read NODAL COORDINATES **************** C
    inf.open(name.c_str());
    if(!inf) {
      std::cout << "Generic-mesh file " << name << " cannot read nodes\n";
      exit(0);
    }
    while(str2.compare("COORDINATES") != 0) inf >> str2;
    inf >> str2;  // 2.0.4
    coords[0].resize(nvt);
    coords[1].resize(nvt);
    coords[2].resize(nvt);

    if(dimNodes == 3) {
      for(unsigned j = 0; j < nvt0; j++) {
        inf >> str2 >> x >> y >> z;
        coords[0][j] = x / Lref;
        coords[1][j] = y / Lref;
        coords[2][j] = z / Lref;
      }
    }
    else if(dimNodes == 2) {
      for(unsigned j = 0; j < nvt0; j++) {
        inf >> str2 >> x >> y;
        coords[0][j] = x / Lref;
        coords[1][j] = y / Lref;
        coords[2][j] = 0.;
      }
    }
    else if(dimNodes == 1) {
      for(unsigned j = 0; j < nvt0; j++) {
        inf >> str2 >> x;
        coords[0][j] = x / Lref;
        coords[1][j] = 0.;
        coords[2][j] = 0.;
      }
    }
    inf >> str2; // "ENDOFSECTION"
    if(str2.compare("ENDOFSECTION") != 0) {
      std::cout << "error node data mesh 1" << std::endl;
      exit(0);
    }
    
    // add the coordinates of the biquadratic nodes not included in gambit
/*    for(int iel = 0; iel < nel; iel++) {
      unsigned elementType = mesh.el->GetElementType(iel);
      
      unsigned ndof = mesh.el->GetElementDofNumber(iel,2);
      unsigned jstart = ndof - _numberOfMissedBiquadraticNodes[elementType];
            
      for(unsigned j = jstart; j < ndof; j++){
	
	unsigned jnode = mesh.el->GetElementDofIndex(iel, j);
	
	coords[0][jnode] = 0.;
        coords[1][jnode] = 0.;
        coords[2][jnode] = 0.;
        for(int i = 0; i < jstart; i++) {
          unsigned inode = mesh.el->GetElementDofIndex(iel, i);
          for(int k = 0; k < mesh.GetDimension(); k++) {
	    coords[k][jnode] += coords[k][inode] * _baricentricWeight[ elementType ][j - jstart][i];
          }          
        }  
      }     
    }*/   

    
    inf.close();
    // end read NODAL COORDINATES ************* C

    // read GROUP **************** E
    inf.open(name.c_str());
    if(!inf) {
      std::cout << "Generic-mesh file " << name << " cannot read group\n";
      exit(0);
    }
    std::vector < unsigned > materialElementCounter(3,0);
    mesh.el->SetElementGroupNumber(ngroup);
    for(unsigned k = 0; k < ngroup; k++) {
      int ngel;
      int gr_name;
      int gr_mat;
      while(str2.compare("GROUP:") != 0) inf >> str2;
      inf >> str2 >> str2 >> ngel >> str2 >> gr_mat >> str2 >> str2 >> gr_name >> str2;
      for(int i = 0; i < ngel; i++) {
        int iel;
        inf >> iel;
        mesh.el->SetElementGroup(iel - 1, gr_name);
        mesh.el->SetElementMaterial(iel - 1, gr_mat);
	if( gr_mat == 2) materialElementCounter[0] += 1;
	else if(gr_mat == 3 ) materialElementCounter[1] += 1;
	else materialElementCounter[2] += 1;
      }
      inf >> str2;
      if(str2.compare("ENDOFSECTION") != 0) {
        std::cout << "error group data mesh" << std::endl;
        exit(0);
      }
    }
    mesh.el->SetMaterialElementCounter(materialElementCounter);
    inf.close();
    // end read GROUP **************** E

    // read boundary **************** D
    inf.open(name.c_str());
    if(!inf) {
      std::cout << "Generic-mesh file " << name << " cannot read boudary\n";
      exit(0);
    }
    for(unsigned k = 0; k < nbcd; k++) {
      while(str2.compare("CONDITIONS") != 0) inf >> str2;
      inf >> str2;
      int value;
      unsigned nface;
      inf >> value >> str2 >> nface >> str2 >> str2;
      value = -value - 1;
      for(unsigned i = 0; i < nface; i++) {
        unsigned iel, iface;
        inf >> iel >> str2 >> iface;
        iel--;
        iface = GambitIO::GambitToFemusFaceIndex[mesh.el->GetElementType(iel)][iface - 1u];
        mesh.el->SetFaceElementIndex(iel, iface, value);
      }
      inf >> str2;
      if(str2.compare("ENDOFSECTION") != 0) {
        std::cout << "error boundary data mesh" << std::endl;
        exit(0);
      }
    }
    inf.close();
    // end read boundary **************** D

  };
  
  
//   void GambitIO::BiquadraticNodesNotInGambit(Mesh& mesh) {
// 
//     unsigned int nnodes = mesh.GetNumberOfNodes();
// //     std::cout << " ********************************** "<< std::endl;
// //     std::cout << "nnodes before = "  << nnodes << std::endl;
// 
//     //intialize to UINT_MAX
//     for( unsigned iel = 0; iel < mesh.el->GetElementNumber(); iel++ ) {
//       unsigned elementType = mesh.el->GetElementType(iel);
//       if( elementType == 1 || elementType == 2 ) {
//         for( unsigned inode = mesh.el->GetElementDofNumber( iel, 2 ) - _numberOfMissedBiquadraticNodes[elementType]; 
// 	    inode < mesh.el->GetElementDofNumber( iel, 2 ); inode++ ) {
//           mesh.el->SetElementDofIndex( iel, inode, UINT_MAX );
//         }
//       }
//     }
// 
//     // generate face dofs for tet and wedge elements
//     for( unsigned iel = 0; iel < mesh.el->GetElementNumber(); iel++ ) {
//       unsigned elementType = mesh.el->GetElementType(iel);
//       if( elementType == 1 || elementType == 2 ) {
//         for( unsigned iface = mesh.el->GetElementFaceNumber( iel, 0 ); iface < mesh.el->GetElementFaceNumber( iel, 1 ); iface++ ) { //on all the faces that are triangles
//           unsigned inode = mesh.el->GetElementDofNumber( iel, 1 ) + iface; 
// 	  if( UINT_MAX == mesh.el->GetElementDofIndex( iel, inode ) ) {
//             mesh.el->SetElementDofIndex( iel, inode, nnodes);
// 	    unsigned i1 = mesh.el->GetFaceVertexIndex( iel, iface, 0 );
//             unsigned i2 = mesh.el->GetFaceVertexIndex( iel, iface, 1 );
// 	    unsigned i3 = mesh.el->GetFaceVertexIndex( iel, iface, 2 );
// 	    bool faceHasBeenFound = false;
// 	    for( unsigned jel = iel + 1; jel < mesh.el->GetElementNumber(); jel++ ) {
// 	      for( unsigned jface = mesh.el->GetElementFaceNumber( jel, 0 ); jface < mesh.el->GetElementFaceNumber( jel, 1 ); jface++ ) {
//                 unsigned jnode = mesh.el->GetElementDofNumber( jel, 1 ) + jface;
// 		if( UINT_MAX == mesh.el->GetElementDofIndex( jel, jnode ) ) {
//                   unsigned j1 = mesh.el->GetFaceVertexIndex( jel, jface, 0 );
//                   unsigned j2 = mesh.el->GetFaceVertexIndex( jel, jface, 1 );
//                   unsigned j3 = mesh.el->GetFaceVertexIndex( jel, jface, 2 );
// 		  if( ( i1 == j1 || i1 == j2 || i1 == j3 ) &&
//                       ( i2 == j1 || i2 == j2 || i2 == j3 ) &&
//                       ( i3 == j1 || i3 == j2 || i3 == j3 ) ) {
//                     mesh.el->SetElementDofIndex( jel, jnode, nnodes);
// 		    faceHasBeenFound = true;
// 		    break;
//                   }
// 		}
// 	      }
// 	      if(faceHasBeenFound){
// 		break;
// 	      }
// 	    }
// 	    ++nnodes; 
// 	  }    
// 	}
//       }
//     }
// 
//     // generates element dofs for tet, wedge and triangle elements
//     for( unsigned iel = 0; iel < mesh.el->GetElementNumber(); iel++ ) {
//       if( 1 == mesh.el->GetElementType( iel ) ) { //tet
//           mesh.el->SetElementDofIndex( iel, 14, nnodes);
// 	  ++nnodes;
//         }
//         else if( 2 == mesh.el->GetElementType( iel ) ) { //wedge
//          mesh.el->SetElementDofIndex( iel, 20, nnodes);
// 	 ++nnodes;
//         }
//         else if( 4 == mesh.el->GetElementType( iel ) ) { //triangle
//           mesh.el->SetElementDofIndex( iel, 6, nnodes);
// 	  ++nnodes;
//         }
//     }
// 
//     mesh.el->SetNodeNumber( nnodes );
//     mesh.SetNumberOfNodes( nnodes );
// //     std::cout <<"nnodes after="<< nnodes << std::endl;
//   }
 
}
