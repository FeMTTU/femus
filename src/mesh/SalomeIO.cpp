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
#include <tuple>

namespace femus
{

  const std::string SalomeIO::group_name_begin = "FAS";
  const std::string SalomeIO::group_name_end   = "ELEME";
  const std::string SalomeIO::mesh_ensemble = "ENS_MAA";
  const std::string SalomeIO::aux_zeroone = "-0000000000000000001-0000000000000000001";
  const std::string SalomeIO::elem_list = "MAI";
  const std::string SalomeIO::group_fam = "FAM";
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


  const unsigned SalomeIO::SalomeToFemusVertexIndex[N_GEOM_ELS][MAX_EL_N_NODES] = {

    {4, 7, 3, 0, 5, 6, 2, 1, 15, 19, 11, 16, 13, 18, 9, 17, 12, 14, 10, 8, 23, 25, 22, 24, 20, 21, 26}, //HEX27
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, //TET10

    {
      3, 11, 5, 9, 10, 4,
      12, 17, 14, 15, 16, 13,
      0, 8, 2, 6, 7, 1
    },  //WEDGE18

    {0, 1, 2, 3, 4, 5, 6, 7, 8}, //QUAD9
    {0, 1, 2, 3, 4, 5},  //TRI6
    {0, 1, 2}            //EDGE3
  };


  const unsigned SalomeIO::SalomeToFemusFaceIndex[N_GEOM_ELS][MAX_EL_N_FACES] = {
    {0, 4, 2, 5, 3, 1},
    {0, 1, 2, 3},
    {2, 1, 0, 4, 3},
    {0, 1, 2, 3},
    {0, 1, 2},
    {0, 1}
  };
 
  
  
  std::pair<int,int>  SalomeIO::isolate_number_in_string(const std::string  string_in, const int begin_pos_to_investigate) {
      
      int str_pos = begin_pos_to_investigate;
      std::cout << " " <<        string_in.at(str_pos) << " ";
      std::string temp_buffer( &(string_in.at(str_pos)) );
      while ( temp_buffer.compare( "_" ) != 0 ) {   str_pos ++;   try { temp_buffer = string_in.at(str_pos); }  catch(const std::out_of_range& e) { std::cerr <<  "agg sbagliat "; abort(); }   }
      
      int str_pos_begin = str_pos;
      
      str_pos = str_pos_begin + 1;
      temp_buffer = string_in.at(str_pos);
      while ( temp_buffer.compare( "_" ) != 0 ) {   str_pos ++;   temp_buffer = string_in.at(str_pos);    }

      int str_pos_end = str_pos;
      
      std::cout <<  str_pos_begin << " " << " " << str_pos_end;    

      std::pair<int,int> output_pair(str_pos_begin,str_pos_end);
      
      return output_pair;
      
  }
  
  

  /// @todo extend to Wegdes (aka Prisms)
  void SalomeIO::read(const std::string& name, vector < vector < double> >& coords, const double Lref, std::vector<bool>& type_elem_flag) {

    Mesh& mesh = GetMesh();
    mesh.SetLevel(0);

    hsize_t dims[2];

    // compute number of menus ===============
    hid_t  file_id = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    hid_t  gid = H5Gopen(file_id, mesh_ensemble.c_str(), H5P_DEFAULT);

    hsize_t     n_meshes_ens;
    hid_t status = H5Gget_num_objs(gid, &n_meshes_ens); // number of meshes
    if(status != 0) {
      std::cout << "Number of mesh menus not found";
      abort();
    }


    // compute number of meshes ===============
    std::vector<std::string>  mesh_menus;

    unsigned n_meshes = 0;

    for(unsigned j = 0; j < n_meshes_ens; j++) {

      char*   menu_names_j = new char[max_length];
      H5Gget_objname_by_idx(gid, j, menu_names_j, max_length); ///@deprecated see the HDF doc to replace this
      std::string tempj(menu_names_j);


      if(tempj.substr(0, 4).compare("Mesh") == 0) {
        n_meshes++;
        mesh_menus.push_back(tempj);
      }
      else { std::cout << "Mesh MED fields must start with the word Mesh" << std::endl; abort();
      }

    }
      if (n_meshes != n_meshes_ens) { std::cout << "Meshes are called with Mesh"; abort(); }
    // compute number of meshes ===============


    // meshes ========================
    for(unsigned j = 0; j < n_meshes_ens; j++) {
        
      std::string tempj = mesh_menus[j];
        
    // ************** Groups of each Mesh *********************************
      std::string group_list = group_name_begin +  "/" + mesh_menus[j] + "/" + group_name_end;
     hsize_t n_groups = 0;
     hid_t  gid_groups   = H5Gopen(file_id, group_list.c_str(), H5P_DEFAULT);
     hid_t status_groups = H5Gget_num_objs(gid_groups, &n_groups);
    if(status_groups != 0) {
      std::cout << "Number of groups not found: You need groups at least for the boundary faces";
      abort();
    }
    
    std::vector<std::string>                     group_names(n_groups);
    std::vector< std::tuple<int,int,int,int> >   group_flags(n_groups);  //salome family; our name; our property; group size
    
     for(unsigned j = 0; j < n_groups; j++) {
        char*   group_names_char = new char[max_length];
        H5Gget_objname_by_idx(gid_groups, j, group_names_char, max_length); ///@deprecated see the HDF doc to replace this
              group_names[j] = group_names_char;
      std::cout << group_names[j] << std::endl;
        delete[] group_names_char;
    

    // read GROUP **************** E
    //we assume that these are VOLUME groups
    //in general, I'd say that a group can only have ONE element type (should study the possibility of hybrid mesh)
   
    std::vector < unsigned > materialElementCounter(3,0);

      const uint n_fe_types_for_groups = 1; // so far we have this assumption

      std::string my_mesh_name_dir = mesh_ensemble +  "/" + group_names[j] + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop

      /// @todo check the underscores according to our naming standard

      int str_pos = 0;
      std::pair<int,int> mypair = isolate_number_in_string(group_names[j],str_pos);
      int gr_family_in_salome = atoi(group_names[j].substr(mypair.first, mypair.second - mypair.first).c_str());
      mypair = isolate_number_in_string(group_names[j],mypair.second+1);
      int gr_name             = atoi(group_names[j].substr(mypair.first, mypair.second - mypair.first).c_str());  //at most 10 groups with this
//       mypair = isolate_number_in_string(group_names[j],mypair.second+1);
      int gr_property         = atoi(group_names[j].substr(mypair.first, mypair.second - mypair.first).c_str());
      int gr_size             = 0;      
      group_flags[j] = std::make_tuple(gr_family_in_salome, gr_name, gr_property,gr_size);

// 
//       std::vector<std::string> el_fe_type(mesh.GetDimension());
// 
//       ReadFE(file_id, el_fe_type, n_fe_types_for_groups, my_mesh_name_dir);
// 
//       std::string group_dataset = mesh_ensemble +  "/" + tempj + "/" +  aux_zeroone + "/" + elem_list + "/" + el_fe_type[mesh.GetDimension() - 1] + "/" + dofobj_indices; ///@todo here we have to loop
// 
//       hid_t dtset = H5Dopen(file_id, group_dataset.c_str(), H5P_DEFAULT);
//       hid_t filespace = H5Dget_space(dtset);
//       hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
//       int* elem_indices = new int[dims[0]];
//       status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, elem_indices);
// 
//       for(unsigned i = 0; i < dims[0]; i++) {
//         mesh.el->SetElementGroup(elem_indices[i] - 1 - n_elements_b_bb, gr_name);
//         mesh.el->SetElementMaterial(elem_indices[i] - 1 - n_elements_b_bb , gr_property);
// 	
// 	    if( gr_property == 2) materialElementCounter[0] += 1;
// 	else if(gr_property == 3 ) materialElementCounter[1] += 1;
// 	else materialElementCounter[2] += 1;
// 	
//       }
// 
//       H5Dclose(dtset);
//       delete [] elem_indices;
// 
//     mesh.el->SetElementGroupNumber(n_gr);
//     mesh.el->SetMaterialElementCounter(materialElementCounter);
    //   // end read GROUP **************** E      
    }
      
      // ************** Groups of each Mesh *********************************

        

      // ************** Mesh *********************************

// dimension ===============
      /// @todo this determination of the dimension from the mesh file would not work with a 2D mesh embedded in 3D
      std::string my_mesh_name_dir = mesh_ensemble +  "/" + tempj + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop

      hsize_t     n_fem_type;
      hid_t       gid = H5Gopen(file_id, my_mesh_name_dir.c_str(), H5P_DEFAULT);
      hid_t status = H5Gget_num_objs(gid, &n_fem_type);
      if(status != 0) {   std::cout << "SalomeIO::read_fem_type:   H5Gget_num_objs not found";  abort();   }
      
      FindDimension(gid, tempj, n_fem_type);
      const unsigned volume_pos = mesh.GetDimension() - 1;
      H5Gclose(gid);

      if(mesh.GetDimension() != n_fem_type) { std::cout << "Mismatch between dimension and number of element types" << std::endl;   abort();  }

// fem type ===============
      std::vector<std::string> el_fe_type(mesh.GetDimension());

      ReadFE(file_id, el_fe_type, n_fem_type, my_mesh_name_dir);

//   // read NODAL COORDINATES **************** C
      std::string coord_dataset = mesh_ensemble +  "/" + tempj + "/" +  aux_zeroone + "/" + node_list + "/" + coord_list + "/";  ///@todo here we have to loop

      hid_t dtset = H5Dopen(file_id, coord_dataset.c_str(), H5P_DEFAULT);

      // SET NUMBER OF NODES
      hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
      hid_t status_dims  = H5Sget_simple_extent_dims(filespace, dims, NULL);
      if(status_dims == 0) std::cerr << "SalomeIO::read dims not found";
      // reading xyz_med
      unsigned int n_nodes = dims[0] / 3; //mesh.GetDimension();
      double*   xyz_med = new double[dims[0]];
      std::cout << " Number of nodes in med file " <<  n_nodes << " " <<  std::endl;

      mesh.SetNumberOfNodes(n_nodes);

      // SET NODE COORDINATES
      coords[0].resize(n_nodes);
      coords[1].resize(n_nodes);
      coords[2].resize(n_nodes);

      H5Dread(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz_med);
      H5Dclose(dtset);

      if(mesh.GetDimension() == 3) {
        for(unsigned j = 0; j < n_nodes; j++) {
          coords[0][j] = xyz_med[j] / Lref;
          coords[1][j] = xyz_med[j + n_nodes] / Lref;
          coords[2][j] = xyz_med[j + 2 * n_nodes] / Lref;
        }
      }

      else if(mesh.GetDimension() == 2) {
        for(unsigned j = 0; j < n_nodes; j++) {
          coords[0][j] = xyz_med[j] / Lref;
          coords[1][j] = xyz_med[j + n_nodes] / Lref;
          coords[2][j] = 0.;
        }
      }

      else if(mesh.GetDimension() == 1) {
        for(unsigned j = 0; j < n_nodes; j++) {
          coords[0][j] = xyz_med[j] / Lref;
          coords[1][j] = 0.;
          coords[2][j] = 0.;
        }
      }

      delete[] xyz_med;

      //   // end read NODAL COORDINATES ************* C


      //   // read ELEMENT/cell ******************** B
       std::vector< unsigned int > n_elems(mesh.GetDimension(),0);
       std::vector< unsigned int > Node_el(mesh.GetDimension(),0);
       
      for(unsigned i = 0; i < mesh.GetDimension(); i++) {
        hsize_t dims_i[2];
        // NOD ***************************
        std::string conn_name_dir_i = my_mesh_name_dir +  el_fe_type[i] + "/" + connectivity;
        hid_t dtset_conn = H5Dopen(file_id, conn_name_dir_i.c_str(), H5P_DEFAULT);
        filespace = H5Dget_space(dtset_conn);
        hid_t status_els_i  = H5Sget_simple_extent_dims(filespace, dims_i, NULL);
        if(status_els_i == 0) { std::cerr << "SalomeIO::read dims not found";   abort();   }
        
            
      // DETERMINE NUMBER OF NODES PER ELEMENT
      Node_el[i] = FindElemNodes(el_fe_type[i]);

      const int dim_conn = dims_i[0];
      n_elems[i] = dim_conn / Node_el[i];
      std::cout << " Number of elements of dimension " << (i+1) << " in med file: " <<  n_elems[i] <<  std::endl;

      // SET NUMBER OF VOLUME ELEMENTS
        if ( i == volume_pos ) { 
      mesh.SetNumberOfElements(n_elems[i]);
      mesh.el = new elem(n_elems[i]);    ///@todo check where this is going to be deleted

      // READ CONNECTIVITY MAP
      int* conn_map = new  int[dim_conn];
      hid_t status_conn = H5Dread(dtset_conn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, conn_map);
      if(status_conn != 0) {     std::cout << "SalomeIO::read: connectivity not found";   abort();   }
      
      
            for(unsigned iel = 0; iel < n_elems[volume_pos]; iel++) {
        mesh.el->SetElementGroup(iel, 1);
        unsigned nve = Node_el[volume_pos];  /// @todo this is only one element type
        if(nve == 27) {
          type_elem_flag[0] = type_elem_flag[3] = true;
          mesh.el->AddToElementNumber(1, "Hex");
          mesh.el->SetElementType(iel, HEX);
        }
        else if(nve == 10) {
          type_elem_flag[1] = type_elem_flag[4] = true;
          mesh.el->AddToElementNumber(1, "Tet");
          mesh.el->SetElementType(iel, TET);
        }
        else if(nve == 18) {
          type_elem_flag[2] = type_elem_flag[3] = type_elem_flag[4] = true;
          mesh.el->AddToElementNumber(1, "Wedge");
          mesh.el->SetElementType(iel, WEDGE);
        }
        else if(nve == 9) {
          type_elem_flag[3] = true;
          mesh.el->AddToElementNumber(1, "Quad");
          mesh.el->SetElementType(iel, QUAD);
        }
        else if(nve == 6 && mesh.GetDimension() == 2) {
          type_elem_flag[4] = true;
          mesh.el->AddToElementNumber(1, "Triangle");
          mesh.el->SetElementType(iel, TRI);
        }
        else if(nve == 3 && mesh.GetDimension() == 1) {
          mesh.el->AddToElementNumber(1, "Line");
          mesh.el->SetElementType(iel, LINE);
        }
        else {
          std::cout << "Error! Invalid element type in reading  File!" << std::endl;
          std::cout << "Error! Use a second order discretization" << std::endl;
          exit(0);
        }
        for(unsigned i = 0; i < nve; i++) {
          unsigned inode = SalomeToFemusVertexIndex[mesh.el->GetElementType(iel)][i];
          mesh.el->SetElementDofIndex(iel, inode, conn_map[iel + i * n_elems[volume_pos]] - 1u);
         }
       }
      
           // clean
      delete [] conn_map;

     }
        
      H5Dclose(dtset_conn);

        
        //NUM ***************************
        std::string node_name_dir_i = my_mesh_name_dir + el_fe_type[i] + "/" + dofobj_indices;
        hid_t dtset_num = H5Dopen(file_id, node_name_dir_i.c_str(), H5P_DEFAULT);
        filespace = H5Dget_space(dtset_num);
        hid_t status_bdry  = H5Sget_simple_extent_dims(filespace, dims_i, NULL);
        if(status_bdry == 0) {    std::cerr << "SalomeIO::read dims not found";  abort();  }
        H5Dclose(dtset_num);
        
        //FAM ***************************
        std::string fam_name_dir_i = my_mesh_name_dir + el_fe_type[i] + "/" + group_fam;
        hid_t dtset_fam = H5Dopen(file_id, fam_name_dir_i.c_str(), H5P_DEFAULT);
        filespace = H5Dget_space(dtset_fam);
        hid_t status_fam  = H5Sget_simple_extent_dims(filespace, dims_i, NULL);
        if(status_fam == 0) {     std::cerr << "SalomeIO::read dims not found";  abort();  }
        
        int* fam_map = new  int[dims_i[0]];
        hid_t status_conn = H5Dread(dtset_conn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fam_map);

        for(unsigned i = 0; i < dims_i[0]; i++) {
          for(unsigned j = 0; j < n_groups; j++) {
          if ( fam_map[i] == std::get<0>(group_flags[j]) ) {
//               std::cout << "Current flag " << fam_map[i] << " matches " << std::get<1>(group_flags[j]) << std::endl; 
              std::get<3>(group_flags[j]) ++;
              std::cout << "Group number for group " << j << " " << std::get<3>(group_flags[j]) << std::endl; 
              
        }
          }  
        }
        
        
        if ( i == (volume_pos - 1) ) {
                    
    //loop over volume elements
    //extract faces

//   // read boundary **************** D
  for (unsigned k=0; k<n_groups/*nbcd*/; k++) { //@todo these should be the groups that are "boundary groups"
           int value = std::get<1>(group_flags[k]);       //flag of the boundary portion
      unsigned nface = std::get<3>(group_flags[k]);  //number of elements in a portion of boundary
               value = - (value + 1);  //@todo so value in salome must be positive
    for (unsigned i = 0; i < nface; i++) {
//       unsigned iel =,   //volume element to which the face belongs
//       iel--;
//       unsigned iface; //index of the face in that volume element
//       iface = SalomeIO::SalomeToFemusFaceIndex[mesh.el->GetElementType(iel)][iface-1u];
//       mesh.el->SetFaceElementIndex(iel,iface,value);  //value is (-1) for element faces that are not boundary faces
    }
  }
//   // end read boundary **************** D
                    
        } //end (volume pos - 1)

        delete [] fam_map;
        H5Dclose(dtset_fam);
        
      }

       //   // end read  ELEMENT/CELL **************** B
      
    }   //end meshes

    H5Fclose(file_id);

  }




  void  SalomeIO::ReadFE(
    hid_t file_id,
    std::vector<std::string>& fe_type_vec,
    hsize_t n_fem_types,
    const std::string  my_mesh_name_dir
  )
  {

    Mesh& mesh = GetMesh();

    // Get the element name
    char** el_fem_type = new char*[n_fem_types];

    std::vector<int> index(n_fem_types);

    const uint fe_name_nchars = 4;

    for(int i = 0; i < (int)n_fem_types; i++) {

      el_fem_type[i] = new char[fe_name_nchars];
      H5Lget_name_by_idx(file_id, my_mesh_name_dir.c_str(), H5_INDEX_NAME, H5_ITER_INC, i, el_fem_type[i], fe_name_nchars, H5P_DEFAULT);
      std::string temp_i(el_fem_type[i]);

      if(mesh.GetDimension() == 3) {

        if(temp_i.compare("HE8") == 0 ||
            temp_i.compare("H20") == 0 ||
            temp_i.compare("H27") == 0 ||
            temp_i.compare("TE4") == 0 ||
            temp_i.compare("T10") == 0) {
          index[mesh.GetDimension() - 1] = i;
          fe_type_vec[mesh.GetDimension() - 1] = el_fem_type[i];
        }
        else if(temp_i.compare("QU4") == 0 ||
                temp_i.compare("QU8") == 0 ||
                temp_i.compare("QU9") == 0 ||
                temp_i.compare("TR3") == 0 ||
                temp_i.compare("TR6") == 0) {
          index[mesh.GetDimension() - 1 - 1] = i;
          fe_type_vec[mesh.GetDimension() - 1 - 1] = el_fem_type[i];
        }
        else if(temp_i.compare("SE2") == 0 ||
                temp_i.compare("SE3") == 0) {
          index[mesh.GetDimension() - 1 - 1 - 1] = i;
          fe_type_vec[mesh.GetDimension() - 1 - 1 - 1] =  el_fem_type[i];
        }

      }

      else if(mesh.GetDimension() == 2) {

        if(temp_i.compare("QU4") == 0 ||
            temp_i.compare("QU8") == 0 ||
            temp_i.compare("QU9") == 0 ||
            temp_i.compare("TR3") == 0 ||
            temp_i.compare("TR6") == 0) {
          index[mesh.GetDimension() - 1] = i;
          fe_type_vec[mesh.GetDimension() - 1] = el_fem_type[i];
        }
        else if(temp_i.compare("SE2") == 0 ||
                temp_i.compare("SE3") == 0) {
          index[mesh.GetDimension() - 1 - 1] = i;
          fe_type_vec[mesh.GetDimension() - 1 - 1] = el_fem_type[i];
        }

      }

      else if(mesh.GetDimension() == 1) {
        if(temp_i.compare("SE2") == 0 ||
            temp_i.compare("SE3") == 0) {
          index[mesh.GetDimension() - 1] = i;
          fe_type_vec[mesh.GetDimension() - 1] = el_fem_type[i];
        }
      }

    }

    // clean
    for(int i = 0; i < (int)n_fem_types; i++) delete[] el_fem_type[i];
    delete[] el_fem_type;

    return;
  }


  void  SalomeIO::FindDimension(hid_t gid, const  std::string menu_name, hsize_t n_fem_type)
  {


    Mesh& mesh = GetMesh();
    uint mydim = 1;  //this is the initial value, then it will be updated below
    mesh.SetDimension(mydim);


    std::vector<char*> elem_list;
    elem_list.resize(n_fem_type);
    for(unsigned j = 0; j < elem_list.size(); j++) {
      elem_list[j] = new char[max_length];
      H5Gget_objname_by_idx(gid, j, elem_list[j], max_length); ///@deprecated see the HDF doc to replace this
      std::string tempj(elem_list[j]);

      if(tempj.compare("HE8") == 0 ||
          tempj.compare("H20") == 0 ||
          tempj.compare("H27") == 0 ||
          tempj.compare("TE4") == 0 ||
          tempj.compare("T10") == 0)      {
        mydim = 3;
      }
      else if(tempj.compare("QU4") == 0 ||
              tempj.compare("QU8") == 0 ||
              tempj.compare("QU9") == 0 ||
              tempj.compare("TR3") == 0 ||
              tempj.compare("TR6") == 0) {
        mydim = 2;
      }

      if(mydim > mesh.GetDimension()) mesh.SetDimension(mydim);

    }  //end for


    return;
  }

  unsigned  SalomeIO::FindElemNodes(const  std::string el_type) const
  {

    unsigned Node_el;

         if(el_type.compare("HE8") == 0) Node_el = 8;
    else if(el_type.compare("H20") == 0) Node_el = 20;
    else if(el_type.compare("H27") == 0) Node_el = 27;

    else if(el_type.compare("TE4") == 0) Node_el = 4;
    else if(el_type.compare("T10") == 0) Node_el = 10;

    else if(el_type.compare("QU4") == 0) Node_el = 4;
    else if(el_type.compare("QU8") == 0) Node_el = 8;
    else if(el_type.compare("QU9") == 0) Node_el = 9;

    else if(el_type.compare("TR3") == 0) Node_el = 3;
    else if(el_type.compare("TR6") == 0) Node_el = 6;

    else if(el_type.compare("SE3") == 0) Node_el = 3;
    else if(el_type.compare("SE2") == 0) Node_el = 2;
    else {
      std::cout << "SalomeIO::read: element not supported";
      abort();
    }

    return Node_el;
  }

} //end namespace femus

