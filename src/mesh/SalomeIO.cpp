/*=========================================================================

 Program: FEMUS
 Module: SalomeIO
 Authors: Sureka Pathmanathan, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//local include
#include "SalomeIO.hpp"
#include "Mesh.hpp"

//C++ include
#include <cassert>
#include <cstdio>
#include <fstream>

namespace femus {
  
   const std::string SalomeIO::group_name_begin = "FAS";
   const std::string SalomeIO::group_name_end   = "ELEME";
   const std::string SalomeIO::mesh_ensemble = "ENS_MAA";
   const std::string SalomeIO::aux_zeroone = "-0000000000000000001-0000000000000000001";
   const std::string SalomeIO::elem_list = "MAI";
   const std::string SalomeIO::connectivity = "NOD";
   const std::string SalomeIO::node_list = "NOE";
   const std::string SalomeIO::coord_list = "COO";
   const std::string SalomeIO::dofobj_indices = "NUM";
   const uint SalomeIO::max_length = 100;  ///@todo this length of the menu string is conservative enough...

  
  //How to determine a general connectivity:
  //you have to align the element with respect to the x-y-z (or xi-eta-zeta) reference frame,
  //and then look at the order in the med file.
  //For every node there is a location, and you have to put that index in that x-y-z location.
  //Look NOT at the NUMBERING, but at the ORDER!
  
// SALOME HEX27  
//         1------17-------5
//        /|              /|
//       / |             / |
//      8  |   21      12  |
//     /   9      22   /   13
//    /    |          /    |
//   0------16-------4     |
//   | 20  |   26    |  25 |
//   |     2------18-|-----6       zeta  
//   |    /          |    /          ^    
//  11   /  24       15  /           |   eta
//   | 10      23    |  14           |  /
//   | /             | /             | /
//   |/              |/              |/ 
//   3-------19------7               -------> xi  

// TEMPLATE HEX27  (for future uses)
//         X------X--------X
//        /|              /|
//       / |             / |
//      X  |   X        X  |
//     /   X      X    /   X
//    /    |          /    |
//   X------X--------X     |
//   | X   |   X     |  X  |
//   |     X------X--|-----X      zeta  
//   |    /          |    /         ^    
//   X   /  X        X   /          |   eta
//   |  X      X     |  X           |  /
//   | /             | /            | /
//   |/              |/             |/ 
//   X-------X-------X              -------> xi   

  
 const unsigned SalomeIO::SalomeToFemusVertexIndex[N_GEOM_ELS][MAX_EL_N_NODES]= 
   {
     
    {4,7,3,0,5,6,2,1,15,19,11,16,13,18,9,17,12,14,10,8,23,25,22,24,20,21,26},  //HEX27
    {0,1,2,3,4,5,6,7,8,9},  //TET10
    
    {
      3, 11,5, 9, 10,4,
      12,17,14,15,16,13,
      0, 8, 2, 6, 7, 1
    },  //WEDGE18
    
    {0,1,2,3,4,5,6,7,8}, //QUAD9
    {0,1,2,3,4,5},       //TRI6
    {0,1,2}              //EDGE3
  };

  
const unsigned SalomeIO::SalomeToFemusFaceIndex[N_GEOM_ELS][MAX_EL_N_FACES]= 
  {
    {0,4,2,5,3,1},
    {0,1,2,3},
    {2,1,0,4,3},
    {0,1,2,3},
    {0,1,2},
    {0,1}
  };

  /// @todo extend to Wegdes (aka Prisms)
void SalomeIO::read(const std::string& name, vector < vector < double> > &coords, const double Lref, std::vector<bool> &type_elem_flag) {
   
    Mesh& mesh = GetMesh();
    mesh.SetLevel(0);

    hsize_t dims[2];
    
   // compute number of menus ===============
    hid_t  file_id = H5Fopen(name.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
    
    hid_t  gid = H5Gopen(file_id,mesh_ensemble.c_str(),H5P_DEFAULT);
      
    hsize_t     n_menus;
    hid_t status= H5Gget_num_objs(gid, &n_menus);  // number of menus
    if(status !=0) { std::cout << "Number of mesh menus not found"; abort(); }
    

    // compute number of groups and number of meshes ===============
    std::vector<std::string>  mesh_menus;
    std::vector<std::string>  group_menus;
    
    unsigned n_groups = 0;
    unsigned n_meshes = 0;
   
    for (unsigned j=0; j<n_menus; j++) {
      
     char *  menu_names_j = new char[max_length];
     H5Gget_objname_by_idx(gid,j,menu_names_j,max_length); ///@deprecated see the HDF doc to replace this
     std::string tempj(menu_names_j);
     
  
     if  (tempj.substr(0,5).compare("Group") == 0) {
       n_groups++;
       group_menus.push_back(tempj);
     }
     else if  (tempj.substr(0,4).compare("Mesh") == 0) {
       n_meshes++;
       mesh_menus.push_back(tempj);
     }
     
    }
   // compute number of groups and number of meshes ===============
    
      unsigned int n_elements_b_bb = 0;
      
   // meshes ======================== 
    for (unsigned j=0; j< n_meshes; j++) {
      
     std::string tempj = mesh_menus[j];

// dimension ===============
     /// @todo this determination of the dimension from the mesh file would not work with a 2D mesh embedded in 3D
  std::string my_mesh_name_dir = mesh_ensemble +  "/" + tempj + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop
  
  hsize_t     n_fem_type;
  hid_t       gid = H5Gopen(file_id,my_mesh_name_dir.c_str(),H5P_DEFAULT);
  hid_t status0 = H5Gget_num_objs(gid, &n_fem_type);
  if(status0 !=0) {std::cout << "SalomeIO::read_fem_type:   H5Gget_num_objs not found"; abort();}
  FindDimension(gid,tempj,n_fem_type);
  H5Gclose(gid);

  if (mesh.GetDimension() != n_fem_type) { std::cout << "Mismatch between dimension and number of element types" << std::endl; abort(); }
  
// fem type ===============
  std::vector<std::string> el_fe_type(mesh.GetDimension());

   ReadFE(file_id,el_fe_type, n_fem_type, my_mesh_name_dir);
     
   //   // read NODAL COORDINATES **************** C
   std::string coord_dataset = mesh_ensemble +  "/" + tempj + "/" +  aux_zeroone + "/" + node_list + "/" + coord_list + "/";  ///@todo here we have to loop

  hid_t dtset = H5Dopen(file_id,coord_dataset.c_str(),H5P_DEFAULT);
  
  // SET NUMBER OF NODES
  hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status ==0) std::cerr << "SalomeIO::read dims not found";
  // reading xyz_med
  unsigned int n_nodes = dims[0]/mesh.GetDimension();
  double   *xyz_med = new double[dims[0]];
  std::cout << " Number of nodes in med file " <<  n_nodes << " " <<  std::endl;
  
  mesh.SetNumberOfNodes(n_nodes);
  
  // SET NODE COORDINATES  
    coords[0].resize(n_nodes);
    coords[1].resize(n_nodes);
    coords[2].resize(n_nodes);

  status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xyz_med);
  H5Dclose(dtset);
  
   if (mesh.GetDimension()==3) {
    for (unsigned j=0; j<n_nodes; j++) {
      coords[0][j] = xyz_med[j]/Lref;
      coords[1][j] = xyz_med[j+n_nodes]/Lref;
      coords[2][j] = xyz_med[j+2*n_nodes]/Lref;
    }
  }

  else if (mesh.GetDimension()==2) {
    for (unsigned j=0; j<n_nodes; j++) {
      coords[0][j] = xyz_med[j]/Lref;
      coords[1][j] = xyz_med[j+n_nodes]/Lref;
      coords[2][j] = 0.;
    }
  }

  else if (mesh.GetDimension()==1) {
    for (unsigned j=0; j<n_nodes; j++) {
      coords[0][j] = xyz_med[j]/Lref;
      coords[1][j] = 0.;
      coords[2][j] = 0.;
    }
  }
  
  delete[] xyz_med;
    
    //   // end read NODAL COORDINATES ************* C

    
    //   // read ELEMENT/cell ******************** B
  std::string node_name_dir = my_mesh_name_dir +  el_fe_type[mesh.GetDimension()-1] + "/" + connectivity;
  hid_t dtset2 = H5Dopen(file_id,node_name_dir.c_str(),H5P_DEFAULT);
  filespace = H5Dget_space(dtset2);
  status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status ==0) {std::cerr << "SalomeIO::read dims not found"; abort();}

  // DETERMINE NUMBER OF NODES PER ELEMENT
  unsigned Node_el = FindElemNodes( el_fe_type[mesh.GetDimension()-1] );
  
  const int dim_conn = dims[0];
  unsigned int n_elements = dim_conn/Node_el;
  
  // SET NUMBER OF ELEMENTS
   mesh.SetNumberOfElements(n_elements);
 
   
  int * conn_map = new  int[dim_conn];
  std::cout << " Number of elements in med file " <<  n_elements <<  std::endl;
  status=H5Dread(dtset2,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,conn_map);
  if(status !=0) {std::cout << "SalomeIO::read: connectivity not found"; abort();}
  H5Dclose(dtset2);

  mesh.el = new elem(n_elements);    ///@todo check where this is going to be deleted
 
   // BOUNDARY (and BOUNDARY of the BOUNDARY in 3D) =========================
 for (unsigned i=0; i < mesh.GetDimension()-1; i++) {
  hsize_t dims_i[2]; 
  std::string node_name_dir_i = my_mesh_name_dir + el_fe_type[i] + "/" + dofobj_indices;
  hid_t dtset_i = H5Dopen(file_id,node_name_dir_i.c_str(),H5P_DEFAULT);
  filespace = H5Dget_space(dtset_i);
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims_i, NULL);
  if(status ==0) { std::cerr << "SalomeIO::read dims not found"; abort(); }
  n_elements_b_bb += dims_i[0];
  H5Dclose(dtset_i);
  }
  // BOUNDARY =========================

  
  for (unsigned iel=0; iel<n_elements; iel++) {
    mesh.el->SetElementGroup(iel,1);
    unsigned nve = Node_el;  /// @todo this is only one element type
    if (nve==27) {
      type_elem_flag[0]=type_elem_flag[3]=true;
      mesh.el->AddToElementNumber(1,"Hex");
      mesh.el->SetElementType(iel,HEX);
    } else if (nve==10) {
      type_elem_flag[1]=type_elem_flag[4]=true;
      mesh.el->AddToElementNumber(1,"Tet");
      mesh.el->SetElementType(iel,TET);
    } else if (nve==18) {
      type_elem_flag[2]=type_elem_flag[3]=type_elem_flag[4]=true;
      mesh.el->AddToElementNumber(1,"Wedge");
      mesh.el->SetElementType(iel,WEDGE);
    } else if (nve==9) {
      type_elem_flag[3]=true;
      mesh.el->AddToElementNumber(1,"Quad");
      mesh.el->SetElementType(iel,QUAD);
    }
    else if (nve==6 && mesh.GetDimension()==2) {
      type_elem_flag[4]=true;
      mesh.el->AddToElementNumber(1,"Triangle");
      mesh.el->SetElementType(iel,TRI);
    }
    else if (nve==3 && mesh.GetDimension()==1) {
      mesh.el->AddToElementNumber(1,"Line");
      mesh.el->SetElementType(iel,LINE);
    } else {
      std::cout<<"Error! Invalid element type in reading  File!"<<std::endl;
      std::cout<<"Error! Use a second order discretization"<<std::endl;
      exit(0);
    }
    for (unsigned i=0; i<nve; i++) {
      unsigned inode = SalomeToFemusVertexIndex[mesh.el->GetElementType(iel)][i];
      mesh.el->SetElementVertexIndex(iel,inode,conn_map[iel+i*n_elements]);
    }
  }

    //   // end read  ELEMENT/CELL **************** B
    
  // clean
  delete [] conn_map;

        }   //end meshes

        
   mesh.el->SetElementGroupNumber(n_groups);
   
       // read GROUP **************** E   
   //we assume that these are VOLUME groups
   //in general, I'd say that a group can only have ONE element type (should study the possibility of hybrid mesh)
   
     for (unsigned j=0; j<n_groups; j++) {
       
       const uint n_fe_types_for_groups = 1; // so far we have this assumption
       
     std::string tempj = group_menus[j];
     std::string my_mesh_name_dir = mesh_ensemble +  "/" + tempj + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop
   
       /// @todo check the underscores according to our naming standard
       
       // strip the first number to get the group number 
       // strip the second number to get the group material
       int gr_name = atoi(tempj.substr(6,1).c_str());  
       int gr_mat =  atoi(tempj.substr(8,1).c_str());

  std::vector<std::string> el_fe_type(mesh.GetDimension());
      
     ReadFE(file_id, el_fe_type, n_fe_types_for_groups, my_mesh_name_dir);
       
    std::string group_dataset = mesh_ensemble +  "/" + tempj + "/" +  aux_zeroone + "/" + elem_list + "/" + el_fe_type[mesh.GetDimension()-1] + "/" + dofobj_indices;  ///@todo here we have to loop
    
  hid_t dtset = H5Dopen(file_id,group_dataset.c_str(),H5P_DEFAULT);
  hid_t filespace = H5Dget_space(dtset);
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  int * elem_indices = new int[dims[0]];
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,elem_indices);

  for (unsigned i=0; i < dims[0]; i++) {
           mesh.el->SetElementGroup(elem_indices[i] -1 - n_elements_b_bb, gr_name);
           mesh.el->SetElementMaterial(elem_indices[i] -1 - n_elements_b_bb ,gr_mat);
    }
	
	
   H5Dclose(dtset);
   delete [] elem_indices;
   
    }
     //   // end read GROUP **************** E
 
 
    status = H5Fclose(file_id);
    
    //loop over volume elements
    //extract faces
    
    
    // 
//   unsigned nbcd;
 
//   // read boundary **************** D
//   inf.open(name.c_str());
//   if (!inf) {
//     std::cout<<"Generic-mesh file "<< name << " cannot read boudary\n";
//     exit(0);
//   }
//   for (unsigned k=0; k<nbcd; k++) {
//     while (str2.compare("CONDITIONS") != 0) inf >> str2;
//     inf >> str2;
//     int value;
//     unsigned nface;
//     inf >>value>>str2>>nface>>str2>>str2;
//     value=-value-1;
//     for (unsigned i=0; i<nface; i++) {
//       unsigned iel,iface;
//       inf>>iel>>str2>>iface;
//       iel--;
//       iface=SalomeIO::SalomeToFemusFaceIndex[mesh.el->GetElementType(iel)][iface-1u];
//       mesh.el->SetFaceElementIndex(iel,iface,value);
//     }
//     inf >> str2;
//     if (str2.compare("ENDOFSECTION") != 0) {
//       std::cout<<"error boundary data mesh"<<std::endl;
//       exit(0);
//     }
//   }
//   inf.close();
//   // end read boundary **************** D
 
}




void  SalomeIO::ReadFE(
  hid_t file_id,
  std::vector<std::string> & fe_type_vec,
  hsize_t n_fem_types,
  const std::string  my_mesh_name_dir
    ) 
{
  
    Mesh& mesh = GetMesh();
   
  // Get the element name
  char **el_fem_type = new char*[n_fem_types];
  
  std::vector<int> index(n_fem_types);

    const uint fe_name_nchars = 4;
  
  for(int i=0; i<(int)n_fem_types; i++) {
    
    el_fem_type[i] = new char[fe_name_nchars];
    H5Lget_name_by_idx(file_id,my_mesh_name_dir.c_str(), H5_INDEX_NAME, H5_ITER_INC,i, el_fem_type[i], fe_name_nchars, H5P_DEFAULT);
    std::string temp_i(el_fem_type[i]);
    
    if (mesh.GetDimension() == 3) {
      
          if ( temp_i.compare("HE8") == 0 || 
	       temp_i.compare("H20") == 0 ||  
	       temp_i.compare("H27") == 0 || 
	       temp_i.compare("TE4") == 0 || 
	       temp_i.compare("T10") == 0   ) { index[mesh.GetDimension() -1] = i; fe_type_vec[mesh.GetDimension() -1] = el_fem_type[i];}
          else if  ( temp_i.compare("QU4") == 0 || 
	             temp_i.compare("QU8") == 0 ||  
	             temp_i.compare("QU9") == 0 || 
	             temp_i.compare("TR3") == 0 || 
	             temp_i.compare("TR6") == 0 ) { index[mesh.GetDimension() -1 -1] = i; fe_type_vec[mesh.GetDimension() -1 -1] = el_fem_type[i]; }
          else if  ( temp_i.compare("SE2") == 0 || 
	             temp_i.compare("SE3") == 0 ) { index[mesh.GetDimension() -1 -1 -1] = i; fe_type_vec[mesh.GetDimension() -1 -1 -1] =  el_fem_type[i]; }
	             
    }
    
    else if (mesh.GetDimension() == 2) {
      
           if  (     temp_i.compare("QU4") == 0 || 
	             temp_i.compare("QU8") == 0 ||  
	             temp_i.compare("QU9") == 0 || 
	             temp_i.compare("TR3") == 0 || 
	             temp_i.compare("TR6") == 0    ) {  index[mesh.GetDimension() -1] = i; fe_type_vec[mesh.GetDimension() -1] = el_fem_type[i]; }
          else if  ( temp_i.compare("SE2") == 0 || 
	             temp_i.compare("SE3") == 0 ) {  index[mesh.GetDimension() -1 -1] = i; fe_type_vec[mesh.GetDimension() -1 -1] = el_fem_type[i]; }
     
    }
    
    else if (mesh.GetDimension() == 1) {
               if  ( temp_i.compare("SE2") == 0 || 
	             temp_i.compare("SE3") == 0 ) { index[mesh.GetDimension() -1] = i; fe_type_vec[mesh.GetDimension() -1] = el_fem_type[i];}
    }    
    
  }

  // clean
  for(int i=0; i<(int)n_fem_types; i++) delete[] el_fem_type[i];
  delete[] el_fem_type;
  
  return;
}


void  SalomeIO::FindDimension(hid_t gid, const  std::string menu_name, hsize_t n_fem_type) {
  
  
    Mesh& mesh = GetMesh();
    uint mydim = 1;  //this is the initial value, then it will be updated below
    mesh.SetDimension(mydim);
    

  std::vector<char*> elem_list;  elem_list.resize(n_fem_type);
    for (unsigned j=0; j < elem_list.size(); j++) {
      elem_list[j] = new char[max_length];
      H5Gget_objname_by_idx(gid,j,elem_list[j],max_length); ///@deprecated see the HDF doc to replace this
      std::string tempj(elem_list[j]);
      
      if ( tempj.compare("HE8") == 0 || 
	   tempj.compare("H20") == 0 ||  
	   tempj.compare("H27") == 0 || 
	   tempj.compare("TE4") == 0 || 
	   tempj.compare("T10") == 0    )      { mydim = 3; } 
      else if ( tempj.compare("QU4") == 0 || 
	        tempj.compare("QU8") == 0 ||  
	        tempj.compare("QU9") == 0 || 
	        tempj.compare("TR3") == 0 || 
	        tempj.compare("TR6") == 0    ) { mydim = 2; } 
      
      if ( mydim > mesh.GetDimension() ) mesh.SetDimension(mydim); 
	
    }  //end for
    
  
  return;
}

unsigned  SalomeIO::FindElemNodes(const  std::string el_type) const {

  unsigned Node_el;
  
        if      ( el_type.compare("HE8") == 0 ) Node_el = 8;
	else if ( el_type.compare("H20") == 0 ) Node_el = 20;  
	else if ( el_type.compare("H27") == 0 ) Node_el = 27; 
	
	else if ( el_type.compare("TE4") == 0 ) Node_el = 4;     
	else if ( el_type.compare("T10") == 0 ) Node_el = 10;     
	
	else if ( el_type.compare("QU4") == 0 ) Node_el = 4;     
	else if ( el_type.compare("QU8") == 0 ) Node_el = 8;     
	else if ( el_type.compare("QU9") == 0 ) Node_el = 9;     
	
	else if ( el_type.compare("TR3") == 0 ) Node_el = 3;     
	else if ( el_type.compare("TR6") == 0 ) Node_el = 6;   

	else if ( el_type.compare("SE3") == 0 ) Node_el = 3;     
	else if ( el_type.compare("SE2") == 0 ) Node_el = 2;    
        else { std::cout << "SalomeIO::read: element not supported"; abort(); }

return Node_el;
}

} //end namespace femus

