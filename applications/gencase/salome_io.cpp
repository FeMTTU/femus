
// This file was massively overhauled

#include "salome_io.hpp"

//C++ includes
#include <fstream>
#include <set>
#include <cstring>

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"

#include "hdf5.h"

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
//       const unsigned int nodes[] =  {0,2,4,6,1,3,5,7,8};
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

    // HEX27 Salome
    {
      eledef.type    = HEX27;
      eledef.dim     = 3;
      eledef.nnodes  = 27;
      eledef.exptype = 12;
      //const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,
      //                            15,16,19,17,18,20,21,24,22,23,25,26
      //                         };
      const unsigned int nodes[] = {0,1,5,4,3,2,6,7,8,13,16,12,10,14,18,15,11,9,17,19,21,20,22,25,24,23,26};
      // const unsigned int nodes[] = {4,16,5,12,21,13,0,8,1,19,25,17,24,26,22,11,20,9,7,18,6,15,23,14,3,10,2};
      //const unsigned int nodes[] = {5,17,6,13,22,14,1,9,2,16,25,18,21,26,23,8,20,10,4,19,7,12,24,15,0,11,3};
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
      eledef.nodes.clear();

      eletypes_exp[PYRAMID5] = eledef;
      eletypes_imp[7]        = eledef;
    }

    //==============================
  }
}
} // end anonymous namespace


// ===================================
void SalomeIO::read_mesh ( std::istream& /*in*/ ){std::cout<< "Use read";abort();}



// ===========================================
// SalomeIO  members
void SalomeIO::read ( const std::string& namefile )
{
  // initialize the map with element types
  init_eletypes();

  // clear any data in the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();
  
  // ***********************************************************
  // ************  start reading  hdf5  ******************************
  // Open an existing file. ---------------
  
  hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dims[2];

  char el_name[4]; uint Node_el;
  H5Lget_name_by_idx(file_id, "/ENS_MAA/Mesh_1/MAI", H5_INDEX_NAME, H5_ITER_INC,0, el_name, 4, H5P_DEFAULT);
  if((std::string) el_name =="HE8") Node_el=8;
  else if((std::string) el_name == "H20") Node_el=20;
  else {std::cout << "Error gmsh_io.C wrong element type"; exit(0);}
 
  std::cout << " Reading data from= " << namefile  <<" with el type " << el_name <<  std::endl;  
 
   // Getting  HEX ++++++++++++++++++++++++++++++++++++++++
   std::string group_name="/ENS_MAA/Mesh_1/MAI/"+(std::string) el_name+"/NOD";
   hid_t  dtset = H5Dopen(file_id,group_name.c_str(), H5P_DEFAULT);
  hid_t  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  H5Sget_simple_extent_dims(filespace, dims, NULL);
  // reading nodes
  unsigned int n_elements =dims[0]/Node_el; int   *conn_map5=new  int[n_elements*Node_el];
  std::cout << " Number of elements " <<  n_elements << " " <<  std::endl;
  H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,conn_map5);
  H5Dclose(dtset); 
   // connectivity table conversion
  uint *conn_map; conn_map = new uint[27*n_elements];
  for (uint j=0;j<n_elements;j++)
    for (uint k=0;k<Node_el;k++) conn_map[27*j+k]= conn_map5[j+k*n_elements]-1;//conn_map5[20*j+k]-1;
  delete[] conn_map5; 
  
  // Getting  xyz ++++++++++++++++++++++++++++++++++++++++
  dtset = H5Dopen(file_id, "/ENS_MAA/Mesh_1/NOE/COO", H5P_DEFAULT);
  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  H5Sget_simple_extent_dims(filespace, dims, NULL);
  // reading xyz
  unsigned int n_nodes =dims[0]/3; double   *xyz5=new double[3*n_nodes];
  std::cout << " Number of points " <<  n_nodes << " " <<  std::endl;
  H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xyz5);
  H5Dclose(dtset);
  // coordinate table conversion
  double *xyz; xyz=new double[3*(n_nodes+(27-Node_el)*n_elements)];
  for (uint j=0;j<n_nodes;j++) {
    xyz[3*j]= xyz5[j]; xyz[3*j+1]= xyz5[j+n_nodes];  xyz[3*j+2]= xyz5[j+2*n_nodes];
  }
  delete[] xyz5;

  

  // ************************* ************************************
  // ********* end reading ************************************


  //temporary global face table
  uint *faces_g; faces_g=new uint [6*n_elements*5];
  uint facecount = 0;
//temporary global edge table

  uint *edges_g; edges_g=new uint [12*n_elements*3];
  uint edgecount = 0;
  
  for (uint ie=0; ie < n_elements; ie++ ) {
    
    if(Node_el ==8){ //*****begin edge **********************
      //edge initialization
    uint edges_l[36] = {conn_map[ie*27+0],conn_map[ie*27+1],0,
			conn_map[ie*27+1],conn_map[ie*27+2],0,
                        conn_map[ie*27+2],conn_map[ie*27+3],0,
                        conn_map[ie*27+3],conn_map[ie*27+0],0,
                        conn_map[ie*27+4],conn_map[ie*27+5],0,
                        conn_map[ie*27+5],conn_map[ie*27+6],0,
                        conn_map[ie*27+6],conn_map[ie*27+7],0,
                        conn_map[ie*27+4],conn_map[ie*27+7],0,
                         conn_map[ie*27+0],conn_map[ie*27+4],0,
                        conn_map[ie*27+1],conn_map[ie*27+5],0,
                        conn_map[ie*27+2],conn_map[ie*27+6],0,
                        conn_map[ie*27+3],conn_map[ie*27+7],0
                       };
    //edge ordering
    for (int i=0;i<12;i++) 
          if (edges_l[3*i] > edges_l[3*i+1]) {
	    uint tmp=edges_l[3*i];edges_l[3*i]=edges_l[3*i+1];edges_l[3*i+1]=tmp;
          }
          
    //edge comparison
    for (int i=0;i<12;i++) {
    uint  edge_exists = 0;uint e=0;
      while(e<edgecount && edge_exists !=2) {
	edge_exists = 0;
        for (uint n=0;n<2;n++) if (edges_l[3*i+n] == edges_g[3*e+n]) edge_exists++;
     e++;
    }
    if(edge_exists == 2 )  edges_l[3*i+2]=edges_g[3*(e-1)+2];
    else {
	double sumx=0;double sumy=0;double sumz=0;
        for (uint n=0;n<2;n++) {
	  edges_g[3*edgecount+n]=edges_l[3*i+n];
	  sumx +=xyz[edges_l[3*i+n]*3];sumy +=xyz[edges_l[3*i+n]*3+1];sumz +=xyz[edges_l[3*i+n]*3+2];
	}
        edges_g[3*edgecount+2] = edges_l[3*i+2] =n_nodes;
	xyz[n_nodes*3]=0.5*sumx;xyz[n_nodes*3+1]=0.5*sumy;xyz[n_nodes*3+2]=0.5*sumz;
        n_nodes++; edgecount++;
      } 
    }

    //fill with edge centers
    for (uint e=0;e<12;e++) conn_map[ie*27+8+e]=edges_l[3*e+2]; 
    } //*************** end edge *********
    
     //************** begin face *********
    //face initialization
    uint faces_l[30] = {conn_map[ie*27+0],conn_map[ie*27+1],conn_map[ie*27+2],conn_map[ie*27+3],0,
                        conn_map[ie*27+0],conn_map[ie*27+1],conn_map[ie*27+5],conn_map[ie*27+4],0,
                        conn_map[ie*27+1],conn_map[ie*27+5],conn_map[ie*27+6],conn_map[ie*27+2],0,
                        conn_map[ie*27+2],conn_map[ie*27+6],conn_map[ie*27+7],conn_map[ie*27+3],0,
                        conn_map[ie*27+0],conn_map[ie*27+4],conn_map[ie*27+7],conn_map[ie*27+3],0,
                        conn_map[ie*27+4],conn_map[ie*27+5],conn_map[ie*27+6],conn_map[ie*27+7],0
                       };
		       
    //face ordering, bubble sort
    for (int i=0;i<6;i++) 
      for (int j=2;j>=0;j--) 
        for (int k=j;k<3;k++) 
          if (faces_l[5*i+k] > faces_l[5*i+k+1]) {
            uint tmp=faces_l[5*i+k];faces_l[5*i+k]=faces_l[5*i+k+1];faces_l[5*i+k+1]=tmp;
          }
          
    //face comparison
    for (int i=0;i<6;i++) {
    uint  face_exists = 0;uint f=0;
      while(f<facecount && face_exists !=4) {
	face_exists = 0;
        for (uint n=0;n<4;n++) if (faces_l[5*i+n] == faces_g[5*f+n]) face_exists++;
     f++;
    }
    if(face_exists == 4 )  faces_l[5*i+4]=faces_g[5*(f-1)+4];
    else {
	double sumx=0;double sumy=0;double sumz=0;
        for (uint n=0;n<4;n++) {
	  faces_g[5*facecount+n]=faces_l[5*i+n];
	  sumx +=xyz[faces_l[5*i+n]*3];sumy +=xyz[faces_l[5*i+n]*3+1];sumz +=xyz[faces_l[5*i+n]*3+2];
	}
        faces_g[5*facecount+4] = faces_l[5*i+4] =n_nodes;
	xyz[n_nodes*3]=0.25*sumx;xyz[n_nodes*3+1]=0.25*sumy;xyz[n_nodes*3+2]=0.25*sumz;
        n_nodes++; facecount++;
      } 
    }

    //fill with face centers
    for (uint f=0;f<6;f++) conn_map[ie*27+20+f]=faces_l[5*f+4];
    
    // fill with volume center
    double sumx=0;double sumy=0;double sumz=0;
    for (uint n=0;n<8;n++) {
	  sumx +=xyz[3*conn_map[ie*27+n]];sumy +=xyz[3*conn_map[ie*27+n]+1];sumz +=xyz[3*conn_map[ie*27+n]+2];
	}
    xyz[n_nodes*3]=0.125*sumx;xyz[n_nodes*3+1]=0.125*sumy;xyz[n_nodes*3+2]=0.125*sumz;
    conn_map[ie*27+26]= n_nodes; n_nodes++;
  } // end element loop

  // clean
  delete [] faces_g; delete [] edges_g;

  
  // libesh storage
  mesh.reserve_nodes (n_nodes);  mesh.reserve_elem (n_elements);

  // coordinates *****************************
  for ( unsigned int i=0; i<n_nodes; ++i ) {
    mesh.add_point ( Point (xyz[i*3],xyz[i*3+1],xyz[i*3+2]),i ); //3D
    // mesh.add_point(Point(xyz[i*3],xyz[i*3+1]),i); //2D
  }
 
  // elements *****************************
  // read the elements
  for ( unsigned int iel=0; iel<n_elements; ++iel ) {

    const elementDefinition& eletype = eletypes_imp[12];// hex27
    // add the elements to the mesh
    Elem* elem = mesh.add_elem ( Elem::build ( eletype.type ).release() );

    // add node pointers to the elements
    for (unsigned int i=0; i<27; i++) {
      elem->set_node(eletype.nodes[i])= mesh.node_ptr(conn_map[27*iel+i]);
     // std::cout <<mesh.n_nodes() << "  " << eletype.nodes[i]  << "  " << conn_map[27*iel+i] << "  " << std::endl;
    }
  }
  delete []xyz; delete []conn_map;

  mesh.print_info();

}

