/*=========================================================================
 *
 * Program: FEMUS
 * Module: MED_IO
 * Authors: Giorgio Bornia, Sureka Pathmanathan
 *
 * Copyright (c) FEMTTU
 * All rights reserved.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the above copyright notice for more information.
 *
 * =========================================================================*/


//local include
#include "MED_IO.hpp"
#include "Mesh.hpp"
#include "GeomElemQuad4.hpp"
#include "GeomElemQuad9.hpp"
#include "GeomElemHex8.hpp"
#include "GeomElemHex27.hpp"
#include "GeomElemTri3.hpp"
#include "GeomElemTri6.hpp"
#include "GeomElemTet4.hpp"
#include "GeomElemTet10.hpp"
#include "GeomElemEdge2.hpp"
#include "GeomElemEdge3.hpp"

#include "H5Opublic.h"

//C++ include
#include <cassert>
#include <cstdio>
#include <fstream>
#include <tuple>
#include <algorithm>



namespace femus {


  const std::string MED_IO::mesh_ensemble                      = "ENS_MAA";
  const std::string MED_IO::aux_zeroone                        = "-0000000000000000001-0000000000000000001";
  const std::string MED_IO::elem_types_folder                  = "MAI";
  const std::string MED_IO::elems_connectivity                 = "NOD";
  const std::string MED_IO::nodes_folder                       = "NOE";
  const std::string MED_IO::nodes_coord_list                   = "COO";
  const std::string MED_IO::node_or_elem_salome_gui_global_num = "NUM";
  const std::string MED_IO::_node_or_elem_group_fam                          = "FAM";  //both for Elements, and for Nodes
  const std::string MED_IO::group_ensemble                     = "FAS";
  const std::string MED_IO::group_elements                     = "ELEME";
  const std::string MED_IO::group_nodes                        = "NOEUD";
  const uint MED_IO::_max_length_med_folder = 300;  ///@todo this length of the menu string is conservative enough...



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


  const unsigned MED_IO::MEDToFemusVertexIndex[N_GEOM_ELS][MAX_EL_N_NODES] = {

    {4, 7, 3, 0, 5, 6, 2, 1, 15, 19, 11, 16, 13, 18, 9, 17, 12, 14, 10, 8, 23, 25, 22, 24, 20, 21, 26}, //HEX27
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, //TET10

    {
      3, 11, 5, 9, 10, 4,
      12, 17, 14, 15, 16, 13,
      0, 8, 2, 6, 7, 1
    },                           //WEDGE18

    {0, 1, 2, 3, 4, 5, 6, 7, 8}, //QUAD9
    {0, 1, 2, 3, 4, 5, 6},       //TRI7
    {0, 1, 2}                    //EDGE3
  };


  const unsigned MED_IO::MEDToFemusFaceIndex[N_GEOM_ELS][MAX_EL_N_FACES] = {
    {0, 4, 2, 5, 3, 1}, //HEX
    {0, 1, 2, 3},       //TET
    {2, 1, 0, 4, 3},    //WEDGE
    {0, 1, 2, 3},       //QUAD
    {0, 1, 2},          //TRI
    {0, 1}              //EDGE
  }; 
  
  
  hid_t MED_IO::open_mesh_file(const std::string& name) {
      
    hid_t  file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    return file_id;
  }
 

 
  void MED_IO::close_mesh_file(hid_t file_id) {
      
        H5Fclose(file_id);
  }
 
 
 
  //template specialization - BEGIN
//template specialization in hpp file, explicit instantiation in cpp file
 template < >  
  void MED_IO::dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(hid_t file_id, std::vector< TYPE_FOR_INT_DATASET > & fam_map, const std::string fam_name_dir_i) const  {
      
       hid_t dtset_fam            = H5Dopen(file_id, fam_name_dir_i.c_str(), H5P_DEFAULT);
      hid_t filespace_fam        = H5Dget_space(dtset_fam);
      hsize_t dims_fam[2];
      hid_t status_fam           = H5Sget_simple_extent_dims(filespace_fam, dims_fam, NULL);
      if(status_fam == 0) {
        std::cerr << "MED_IO::read dims not found";
        abort();
      }

      const unsigned n_elements = dims_fam[0];
//       std::vector< TYPE_FOR_INT_DATASET > fam_map(n_elements);
      fam_map.resize(n_elements);
      hid_t status_conn = H5Dread(dtset_fam,/*ONLY DIFFERENCE*/H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fam_map.data());
      H5Dclose(dtset_fam);     
      
  }


  
 template < >  
  void MED_IO::dataset_open_and_close_store_in_vector<TYPE_FOR_REAL_DATASET>(hid_t file_id, std::vector< TYPE_FOR_REAL_DATASET > & fam_map, const std::string fam_name_dir_i) const  {
      
       hid_t dtset_fam            = H5Dopen(file_id, fam_name_dir_i.c_str(), H5P_DEFAULT);
      hid_t filespace_fam        = H5Dget_space(dtset_fam);
      hsize_t dims_fam[2];
      hid_t status_fam           = H5Sget_simple_extent_dims(filespace_fam, dims_fam, NULL);
      if(status_fam == 0) {
        std::cerr << "MED_IO::read dims not found";
        abort();
      }

      const unsigned n_elements = dims_fam[0];
//       std::vector< TYPE_FOR_REAL_DATASET > fam_map(n_elements);
      fam_map.resize(n_elements);
      hid_t status_conn = H5Dread(dtset_fam,/*ONLY DIFFERENCE*/ H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fam_map.data());
      H5Dclose(dtset_fam);     
      
  }  
   
  
  
// //explicit instantiation==============
//   
//  template < >  
//   void MED_IO::dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(hid_t file_id, std::vector< TYPE_FOR_INT_DATASET > & fam_map, const std::string fam_name_dir_i) const;
// 
//       template < >  
//   void MED_IO::dataset_open_and_close_store_in_vector<TYPE_FOR_REAL_DATASET>(hid_t file_id, std::vector< TYPE_FOR_REAL_DATASET > & fam_map, const std::string fam_name_dir_i) const;
  
    //template specialization - END


 
 
 

  
  /// @todo extend to Wegdes (aka Prisms)
  /// @todo why pass coords other than get it through the Mesh class pointer?
  void MED_IO::read(const std::string& name, 
                    std::vector < std::vector < double> >& coords, 
                    const double Lref, 
                    std::vector<bool>& type_elem_flag, 
                    const bool read_domain_groups_flag, 
                    const bool read_boundary_groups_flag) {


    Mesh& mesh = GetMesh();


// === Level - BEGIN =================
    mesh.SetLevel(0);
// === Level - END =================

    hid_t  file_id = open_mesh_file(name);

    // meshes, read - BEGIN ========================
    const std::vector< std::string > mesh_menus = get_mesh_names(file_id);
    // meshes, read - END ========================


    // dimension - BEGIN ===============
    const std::string elem_list_dir = get_element_info_all_dims_H5Group( mesh_menus[0] );  ///@todo here we have to loop
    
    const unsigned dim_from_attr = get_mesh_dimension_from_attributes(file_id,  mesh_menus[0] );
    
    const unsigned dim_from_elem_types = get_mesh_dimension_by_looping_over_geom_elem_types( file_id, elem_list_dir);

    if ( dim_from_attr != dim_from_elem_types) abort();
    
    mesh.SetDimension(dim_from_elem_types);
    
    const unsigned mesh_dim = mesh.GetDimension();
    // dimension - END ===============

    
    // geom_el types - BEGIN ===============
    
       //  std::vector std::vector< GeomElemBase* > geom_elem_per_dimension_vec(mydim);
    // geom_elem_per_dimension_vec[d]

    const std::vector< GeomElemBase* >  geom_elem_per_dimension = get_geom_elem_type_per_dimension(file_id, elem_list_dir);

    //        if(mesh.GetDimension() != n_fem_type) { std::cout << "Mismatch between dimension and number of element types" << std::endl;   abort();  }
    // ///@todo removed this check to allow 2d object in 3d

    // geom_el types - END ===============


    // meshes - BEGIN ========================
    for(unsigned j = 0; j < mesh_menus.size(); j++) {

      // nodes, coordinates - BEGIN
      set_node_coordinates(file_id, mesh_menus[j], coords, Lref);
      // nodes, coordinates - END

      // Elements, Volume connectivity - BEGIN
      //       for(unsigned i = 0; i < mesh_dim; i++) {
      unsigned i = mesh_dim - 1;
      set_elem_connectivity_and_initialize_elem_group(file_id, mesh_menus[j], i, geom_elem_per_dimension[i], type_elem_flag);  //type_elem_flag is to say "There exists at least one element of that type in the mesh"
      // Elements, Volume connectivity - END



       // Groups and Material of the mesh - BEGIN ===============
      
      if(read_domain_groups_flag == true || read_boundary_groups_flag == true)  {

        // Group info - BEGIN ===============
        std::vector< GroupInfo >     group_info = get_all_groups_per_mesh(file_id, mesh_menus[j]);
        // Group info - END ===============

        // Group Geom Elem and Size - BEGIN ===============
        for(unsigned i = 0; i < group_info.size(); i++) {
          node_or_elem_Compute_group_geometric_object_type_and_size(file_id, mesh_menus[j], group_info[i]);
        }
        // Group Geom Elem and Size - END ===============

        // Domain groups (dimension same as max mesh dimension) - BEGIN ===============
        if(read_domain_groups_flag == true) {
          if(i == (mesh_dim - 1)) {
//           for(unsigned i = 0; i < mesh_dim; i++) {
            set_elem_group_ownership_And_Material(file_id, mesh_menus[j], group_info, geom_elem_per_dimension[i], i);
          }
        }
        // Domain groups (dimension same as max mesh dimension) - END ===============

        // Boundary groups (dimension - 1) - BEGIN ===============
        if(read_boundary_groups_flag == true)  {
          set_elem_group_ownership_boundary(file_id, mesh_menus[j], group_info, geom_elem_per_dimension, mesh_dim, mesh.el->GetElementNearFaceArray());
        }
        // Boundary groups (dimension - 1) - END ===============

      }
       // Groups and Material of the mesh - END ===============


    }
    
    // meshes - END ========================


    close_mesh_file(file_id);

  }



  // This function reads all boundary groups for the boundary conditions
  // I need a function that reads only one given boundary group, characterized by a name
  // I should distinguish the boundary groups related to the boundary conditions with some name, such as "Boundary_1_0" instead of "Group_1_0"
  void MED_IO::set_elem_group_ownership_boundary(const hid_t  file_id,
                                                 const std::string mesh_menu,
                                                 const std::vector< GroupInfo > & group_info,
                                                 const std::vector< GeomElemBase* > geom_elem_per_dimension,
                                                 const unsigned mesh_dim,
                                                 MyMatrix <int> & element_faces_array
                                                ) {

    if(mesh_dim > 1)       find_boundary_faces_and_set_face_flags(file_id, mesh_menu, group_info, geom_elem_per_dimension[mesh_dim - 1 - 1], element_faces_array);
    else if(mesh_dim == 1)   find_boundary_nodes_and_set_node_flags(file_id, mesh_menu, group_info, element_faces_array);

  }





  unsigned int MED_IO::get_user_flag_from_med_flag(const std::vector< GroupInfo > & group_info, const TYPE_FOR_INT_DATASET med_flag_in) const {

    unsigned int user_flag;

    for(unsigned gv = 0; gv < group_info.size(); gv++) {
      if(group_info[gv]._med_flag == med_flag_in) user_flag = group_info[gv]._user_defined_flag;
    }

    return user_flag;

  }
  

  unsigned int MED_IO::get_med_flag_from_user_flag(const std::vector< GroupInfo > & group_info, const TYPE_FOR_INT_DATASET input_flag) const {

    unsigned int output_flag = 0;

    for(unsigned gv = 0; gv < group_info.size(); gv++) {
      if(group_info[gv]._user_defined_flag == input_flag) output_flag = group_info[gv]._med_flag;
    }

    return output_flag;

  }
  


  
  
  //here I need a routine to compute the group GeomElem and the group size

  //separate groups by dimension
  //as soon as an entry is equal to the group_med_flag, that means the dimension is that of the current element dataset
  void MED_IO::node_or_elem_Compute_group_geometric_object_type_and_size(const hid_t&  file_id, 
                                                const std::string mesh_menu,
                                                GroupInfo & group_info) const {
                                                    

    if (group_info._node_or_cell_group == group_nodes ) {
         
// ========= group geom elem - BEGIN ==========
        std::string  elem_types_str = nodes_folder;
        
       group_info._geom_el = get_geom_elem_from_med_name(elem_types_str);
// ========= group geom elem - END ==========
 
// ========= group size - BEGIN ==========
     std::string my_mesh_name_dir = get_node_info_H5Group(mesh_menu);
 
     std::string fam_name_dir = my_mesh_name_dir  + _node_or_elem_group_fam + "/";
      
     std::vector< TYPE_FOR_INT_DATASET > fam_map;
      
      dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(file_id, fam_map, fam_name_dir);
      
       int group_size = 0;
      for(unsigned k = 0; k < fam_map.size(); k++) {

        if(fam_map[k] == group_info._med_flag)   {
                group_size++;
          if(_print_info) {
            std::cout << "Current flag " << fam_map[k] << " matches " << group_info._med_flag << " node" << std::endl;
                }
           }
       }
       
       group_info._size = group_size;
// ========= group size - END ==========
         
     }
     
     
    else if (group_info._node_or_cell_group == group_elements)  {
      
      std::string my_mesh_name_dir = get_element_info_all_dims_H5Group(mesh_menu);  ///@todo here we have to loop
        
    hid_t       gid = H5Gopen(file_id, my_mesh_name_dir.c_str(), H5P_DEFAULT);
    
    hsize_t    n_geom_el_types = get_H5G_size(gid);

    std::vector< std::string > elem_types(n_geom_el_types);

    bool group_found = false;

    //loop over all FAM fields until the group is found
    unsigned j = 0;
    while(j < elem_types.size() && group_found == false) {

      std::string elem_types_str =  get_H5L_name_by_idx(gid, ".", j);
      elem_types[j] = elem_types_str;

      std::string fam_name_dir_i = my_mesh_name_dir + elem_types_str + "/" + _node_or_elem_group_fam;
      
      std::vector< TYPE_FOR_INT_DATASET > fam_map;
      
      dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(file_id, fam_map, fam_name_dir_i);

      

      int group_size = 0;
      for(unsigned k = 0; k < fam_map.size(); k++) {

        if(fam_map[k] == group_info._med_flag)   {
          group_found = true;
          group_size++;
          if(_print_info) {
            std::cout << "Current flag " << fam_map[k] << " matches " << group_info._med_flag << " " << elem_types[j] << std::endl;
          }
          group_info._geom_el = get_geom_elem_from_med_name(elem_types_str);
        }

      }
      group_info._size = group_size;


      j++;
    }

    H5Gclose(gid);

    
  }
  
  
    return;

  }


   std::string MED_IO::get_element_info_all_dims_H5Group(const std::string mesh_menu) const {
    
     return  mesh_ensemble +  "/" + mesh_menu + "/" +  aux_zeroone + "/" + elem_types_folder + "/";
       
   }
   

   std::string MED_IO::get_node_info_H5Group(const std::string mesh_menu) const {
    
     return mesh_ensemble +  "/" + mesh_menu + "/" +  aux_zeroone + "/" + nodes_folder + "/";
       
   }
   
   
   std::string MED_IO::get_group_info_H5Group(const std::string mesh_menu, const std::string geom_elem_type) const {
       
    std::string group_list =  group_ensemble +  "/" + mesh_menu + "/" + geom_elem_type;
    
    return  group_list;
    
   }

   
  void MED_IO::find_boundary_faces_and_set_face_flags(const hid_t&  file_id,
                                                      const std::string mesh_menu,
                                                      const std::vector< GroupInfo > & group_info,
                                                      const GeomElemBase* geom_elem_per_dimension,
                                                      MyMatrix <int> & element_faces_array)  {
    //basically you have certain elements on the boundary, and you have to find them in the list of all elements


    std::string my_mesh_name_dir = get_element_info_all_dims_H5Group(mesh_menu);  ///@todo here we have to loop

    //open the NOD field of the boundary element list (for the connectivities)
    std::string   conn_name_dir = my_mesh_name_dir + geom_elem_per_dimension->get_name_med() + "/" + elems_connectivity; ///@todo these boundary connectivities were not stored, so we need to read them now
    
      std::vector< TYPE_FOR_INT_DATASET > conn_map;
      
      dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(file_id, conn_map, conn_name_dir);
      

    //open the FAM field of the boundary element list (for the group flags)
    std::string fam_name_dir = my_mesh_name_dir + geom_elem_per_dimension->get_name_med() + "/" + _node_or_elem_group_fam;
    
      std::vector< TYPE_FOR_INT_DATASET > fam_map;
      
      dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(file_id, fam_map, fam_name_dir);
    

    //check that all boundary faces were set in the mesh file
    for(unsigned i = 0; i < fam_map.size(); i++) {
      if(fam_map[i] == 0) {
        std::cout << "Some boundary elements were not set in the mesh MED file. If you believe all boundary elements belong to a Group, the error could be due to the presence of inner boundary elements. Please remove them in your mesh. " << std::endl;
        abort();
      }
    }


    // loop over elements to find faces
    Mesh& mesh = GetMesh();

    unsigned count_found_face = 0;

    //loop over the volume connectivity and find the boundary faces

    for(unsigned iel = 0; iel < mesh.GetNumberOfElements(); iel++) {

      unsigned iel_geom_type = mesh.GetElementType(iel);

      for(unsigned f = 0; f < mesh.GetElementFaceNumber(iel); f++) {

        unsigned n_nodes_face = _geom_elems[iel_geom_type]->get_nodes_of_face(f).size();

        // on one hand I construct the boundary face connectivity from the volume connectivity, with the order that was given in our code
        std::vector<unsigned> face_nodes_from_vol_connectivity(n_nodes_face);

        for(unsigned nd = 0; nd < n_nodes_face; nd++) {
          unsigned nd_of_face = _geom_elems[iel_geom_type]->get_nodes_of_face(f)[nd];
          face_nodes_from_vol_connectivity[nd] = mesh.el->GetElementDofIndex(iel, nd_of_face);

        }

        // on the other I read the boundary face connectivity from the list of boundary faces, which is the one that has the FAM information
        //loop over all the bdry group elements

        std::vector<unsigned> face_nodes_from_bdry_group(n_nodes_face);

        for(unsigned k = 0; k < fam_map.size(); k++) {

          for(unsigned nd = 0; nd < n_nodes_face; nd++) {
            face_nodes_from_bdry_group[nd] = conn_map[ k + nd * fam_map.size() ] - 1;
          }

          bool are_faces_the_same = see_if_faces_from_different_lists_are_the_same(geom_elem_per_dimension, face_nodes_from_vol_connectivity, face_nodes_from_bdry_group);

                 const TYPE_FOR_INT_DATASET med_flag = fam_map[k];
                 
  if (med_flag > 0) abort(); //MED puts its own NEGATIVE flags for elements in dim >=1 (edges, faces), and POSITIVE for nodes (dim=0)

          if(are_faces_the_same)  {
            count_found_face++;


            int user_flag =  get_user_flag_from_med_flag(group_info, med_flag);   //flag of the boundary portion
            user_flag = - (user_flag + 1);  ///@todo these boundary indices need to be NEGATIVE,  so the user_flag in salome must be POSITIVE

            if(_print_info) {
              std::cout << "Found face " << k << " in element " << iel << " with MED flag " << med_flag << " and user flag " << user_flag << std::endl;
            }

            //       unsigned iface = MED_IO::MEDToFemusFaceIndex[mesh.el->GetElementType(iel)][iface-1u];//index of the face in that volume element
            element_faces_array[iel][f] = user_flag;   //user_flag is (-1) for element faces that are not boundary faces, SO WE MUST BE CAREFUL HERE!
            //                 mesh.el->SetFaceElementIndex(iel, f, user_flag);  //old version

          }

        }
        //loop over all the bdry group elements - end


      } //faces of volume elements

    } // end volume elements


    check_all_boundaries_have_some_condition(count_found_face, fam_map.size());

  }


  void MED_IO::check_all_boundaries_have_some_condition(const unsigned int count_found_face, const unsigned int fam_size) const {
    //Checking that all faces have some boundary condition specified

    if(_print_info) {
      std::cout << "Count found faces " << count_found_face << std::endl;
    }
    if(count_found_face < fam_size) {
      std::cout << "Found " << count_found_face << " faces out of " << fam_size << ", not enough: missing certain boundary conditions." << std::endl;
      abort();
    }

  }

  //this is for 1D domains
  void MED_IO::find_boundary_nodes_and_set_node_flags(const hid_t&  file_id,
                                                      const std::string mesh_menu,
                                                      const std::vector<GroupInfo> & group_info,
                                                      MyMatrix <int> & element_faces_array)  {

    std::string my_mesh_name_dir = get_node_info_H5Group(mesh_menu);  ///@todo here we have to loop

    //open the FAM field of NOE
    std::string fam_name_dir = my_mesh_name_dir /*+ geom_elem_per_dimension->get_name_med()*/ + "/" + _node_or_elem_group_fam;
    
    std::vector< TYPE_FOR_INT_DATASET > fam_map;
      
    dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(file_id, fam_map, fam_name_dir);


    // I have to loop over elements to then find faces of the elements, get the dof of those faces, and then with that dof go get the position in NOE/FAM 
       
       
            // loop over elements to find faces
       Mesh& mesh = GetMesh();
       
       //loop over the volume connectivity and find the boundary faces 
       //the boundary faces of each volume element have already been constructed after the mesh reading
       
        for(unsigned iel = 0; iel < mesh.GetNumberOfElements(); iel++) {
            
//             unsigned iel_geom_type = mesh.GetElementType(iel);
                               
            for(unsigned f = 0; f < mesh.GetElementFaceNumber(iel); f++) {
                
                unsigned n_nodes_in_face = 1 /*_geom_elems[iel_geom_type]->get_nodes_of_face(f).size()*/;  //here the faces are made of 1 node
                
                std::vector<unsigned> face_nodes(n_nodes_in_face);
                
               for(unsigned nd = 0; nd < face_nodes.size(); nd++) {
                   unsigned nd_of_face = f /*_geom_elems[iel_geom_type]->get_nodes_of_face(f)[nd]*/;
                   face_nodes[nd] = mesh.el->GetElementDofIndex(iel, nd_of_face);
                
               }
            
            // now I take this dof and read from NOE/FAM  
            // If the group is different from zero 
               for(unsigned k = 0; k < fam_map.size(); k++) {

                      const TYPE_FOR_INT_DATASET med_flag = fam_map[k];
                      
          if (med_flag < 0) abort(); //MED puts its own NEGATIVE flags for elements in dim >=1 (edges, faces), and POSITIVE for nodes (dim=0)
  
              if ( face_nodes[0] == k && med_flag != 0 )  {
                      
           int user_flag =  get_user_flag_from_med_flag(group_info, med_flag);   //flag of the boundary portion
               user_flag = - (user_flag + 1);  ///@todo these boundary indices need to be NEGATIVE,  so the user_flag in salome must be POSITIVE
               
                     if (_print_info) { std::cout << "Found face " << k << " in element " << iel << " with MED flag " << med_flag << " and user flag " << user_flag << std::endl; }

//       unsigned iface = MED_IO::MEDToFemusFaceIndex[mesh.el->GetElementType(iel)][iface-1u];//index of the face in that volume element
                  element_faces_array[iel][f] = user_flag;  //user_flag is (-1) for element faces that are not boundary faces, SO WE MUST BE CAREFUL HERE!
              //  mesh.el->SetFaceElementIndex(iel, f, user_flag);  //old version

               }
      
      
                } //end k  
                
           } //faces loop
           
      } 
       
       
   }

   
   bool MED_IO::boundary_of_boundary_3d_check_face_of_face_via_nodes(const std::vector < int > nodes_face_face_flags, const unsigned group_salome)  {

                   bool is_face_bdry_bdry = false;
              unsigned int i_bdry_bdry = 0;
              while ( i_bdry_bdry < nodes_face_face_flags.size() ) {
                  if (nodes_face_face_flags[i_bdry_bdry] == group_salome) { i_bdry_bdry++;}
                  else break;
               }
                if  (  i_bdry_bdry ==   nodes_face_face_flags.size() ) { is_face_bdry_bdry = true; }

           return  is_face_bdry_bdry;    
                
  }
  
  
  //the node global ordering is given by the mesh file, as well as the element global ordering
  void MED_IO::all_nodes_read_group_flag(const hid_t&  file_id, const std::string mesh_menu,  std::vector < TYPE_FOR_REAL_DATASET >  & node_group_map) {

    
    std::string node_group_dataset = get_node_info_H5Group(mesh_menu) + _node_or_elem_group_fam + "/";
              
    dataset_open_and_close_store_in_vector<TYPE_FOR_REAL_DATASET>(file_id, node_group_map, node_group_dataset);
 
  }

  
  /// @todo Pay attention that here it is all in double
   std::vector< TYPE_FOR_REAL_DATASET > MED_IO::node_based_flag_read_from_file(const std::string& name, const std::vector< unsigned > & mapping) {
       
       if (GetMesh().GetDimension() != 3 ) abort();
       
       
        hid_t  file_id = open_mesh_file(name);
        
        const std::vector< std::string > mesh_menus = get_mesh_names(file_id);

        std::vector< TYPE_FOR_REAL_DATASET > node_group_map_with_med_ordering;

        all_nodes_read_group_flag(file_id, mesh_menus[0],  node_group_map_with_med_ordering);
               
        close_mesh_file(file_id);
        
        //reorder with femus mapping - BEGIN
       std::vector < TYPE_FOR_REAL_DATASET > node_group_map_with_femus_ordering(GetMesh().GetNumberOfNodes());


        for(unsigned j = 0; j < node_group_map_with_med_ordering.size(); j++) {
          node_group_map_with_femus_ordering[ mapping[j] ] = node_group_map_with_med_ordering[j];
        }
        //reorder with femus mapping - END
       
        
        return node_group_map_with_femus_ordering;
    } 
 
  
  //this is for 3D domains
   void MED_IO::boundary_of_boundary_3d_via_nodes(const std::string& name, const unsigned group_user) {
       
       
      // ======= FILE READ - BEGIN  ==================
       
        hid_t  file_id = open_mesh_file(name);
        
        const std::vector< std::string > mesh_menus = get_mesh_names(file_id);

     std::vector< TYPE_FOR_INT_DATASET > node_group_map;

    std::string node_group_dataset = get_node_info_H5Group(mesh_menus[0]) + _node_or_elem_group_fam + "/";
              
    dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(file_id, node_group_map, node_group_dataset);
        
        
      // ======= group that one wants to find (have to put the user flag and convert it) ==================
       // Groups of the mesh ===============
        std::vector< GroupInfo >     group_info = get_all_groups_per_mesh(file_id, mesh_menus[0]);
        
        const unsigned group_salome = get_med_flag_from_user_flag(group_info, group_user);
     // =========================

               
        close_mesh_file(file_id);
      // ======= FILE READ - END ==================

        
        
         Mesh& mesh = GetMesh();
       if (mesh.GetDimension() != 3 ) abort();
       
   
 unsigned solType_coords = 2;
       
        for (unsigned iel = 0; iel < mesh.GetNumberOfElements(); iel++) {
            
      unsigned iel_geom_type = mesh.GetElementType(iel);
      unsigned iel_n_faces = mesh.GetElementFaceNumber(iel);
      
      for (unsigned f = 0; f < iel_n_faces; f++) {

          
      unsigned iel_geom_type_face = mesh.GetElementFaceType(iel, f);
 
      unsigned f_n_faces_faces  =  mesh.el->GetNFC(iel_geom_type, iel_geom_type_face); /* ElementFaceFaceNumber */
//           unsigned n_nodes_face = _geom_elems[iel_geom_type]->get_nodes_of_face(f).size();
 
      
      for (unsigned f_f = 0; f_f < f_n_faces_faces; f_f++) {
          
          
                  unsigned n_nodes_face_face = _geom_elems[iel_geom_type_face]->get_nodes_of_face(f_f).size();

                  std::vector < int > nodes_face_face_flags(n_nodes_face_face, 0); 
          
                                               
          
		  for (unsigned i_bdry_bdry = 0; i_bdry_bdry < n_nodes_face_face; i_bdry_bdry++) {
               
    unsigned LocalFaceFaceVertexIndex = mesh.el->GetIG(iel_geom_type_face, f_f, i_bdry_bdry);    //from n-2 to n-1
    unsigned LocalFaceVertexIndex    = mesh.el->GetIG(iel_geom_type, f, LocalFaceFaceVertexIndex); //from n-1 to n
 
                unsigned int i_vol_iel = LocalFaceVertexIndex; //mesh.GetLocalFaceVertexIndex(iel, f, i_bdry);
                
                //here is where I want to go from LOCAL to GLOBAL mesh node
  unsigned node_global = mesh.el->GetElementDofIndex(iel, i_vol_iel);
  
   if (node_group_map[node_global] == group_salome) std::cout << node_global << std::endl;
    
       nodes_face_face_flags[i_bdry_bdry] = node_group_map[node_global];
     
              }
              
                  bool is_face_bdry_bdry  =  boundary_of_boundary_3d_check_face_of_face_via_nodes( nodes_face_face_flags, group_salome);
                  
                 if  ( is_face_bdry_bdry ) { std::cout << " Face in iel ======= " << iel << " face " << f << " face_face " <<  f_f << std::endl;  }

                
              
            }
            
            
      }
      
    }
    
    
       
   }
  
   
   
  bool MED_IO::see_if_faces_from_different_lists_are_the_same( const GeomElemBase* geom_elem_per_dimension, 
                                                   const std::vector< unsigned > & face_nodes_from_vol_connectivity, 
                                                   const std::vector< unsigned > & face_nodes_from_bdry_group) {

                     // check any possible order of faces 
                  
                  //just look for the initial linear element (maybe even the 1st three only); if this is aligned, all the Quad9 will be aligned
                  //The problem, is that we don't know if the order corresponds to the OUTWARD NORMAL or not.
                  //How many ways are there? If the nodes are 0 1 2 3, it could only be 0123, or 1230, or 2301, or 3012, or the REVERSE of each of them
                  std::vector<unsigned>  face_nodes_from_vol_connectivity_linear(face_nodes_from_vol_connectivity.begin(),face_nodes_from_vol_connectivity.begin() + geom_elem_per_dimension->n_nodes_linear() );
                  std::vector<unsigned>  face_nodes_from_bdry_group_linear(face_nodes_from_bdry_group.begin(),face_nodes_from_bdry_group.begin() + geom_elem_per_dimension->n_nodes_linear() );
                  
                  unsigned n_alternatives = 2 * geom_elem_per_dimension->n_nodes_linear();
                  std::vector< std::vector<unsigned> > face_alternatives(n_alternatives);
                  
               for(unsigned alt = 0; alt < n_alternatives/2; alt++) {
                   unsigned face_length = geom_elem_per_dimension->n_nodes_linear();
                   face_alternatives[alt].resize( face_length );
                 for(unsigned i = 0; i < face_length; i++) {
                     unsigned mod_index = (i+alt)%face_length;
                          face_alternatives[alt][i] = face_nodes_from_bdry_group_linear[mod_index];
                   }
               }
              
               for(unsigned alt = n_alternatives/2; alt < n_alternatives; alt++) {
                    face_alternatives[alt] = face_alternatives[alt-(n_alternatives/2)];
                    std::reverse(face_alternatives[alt].begin(),face_alternatives[alt].end());
                }

          
             std::vector<bool>  is_same_face(n_alternatives);
             bool bool_union = false;
                      for(unsigned alt = 0; alt < n_alternatives; alt++) {
                            is_same_face[alt] = (face_nodes_from_vol_connectivity_linear == face_alternatives[alt]);
                            if (is_same_face[alt] == true) bool_union = true;
                       }
                      
                      return bool_union;
                                        
  }




  //  we loop over all elements and see which ones are of that group
  //
  void MED_IO::set_elem_group_ownership_And_Material(
    const hid_t&  file_id,
    const std::string mesh_menu,
    const std::vector<GroupInfo> & group_info,
    const GeomElemBase* geom_elem_per_dimension,
    const int i
  )  {

    Mesh& mesh = GetMesh();

    //FAM ***************************
    
    std::string my_mesh_name_dir = get_element_info_all_dims_H5Group(mesh_menu);  ///@todo here we have to loop
    
    std::string fam_name_dir_i = my_mesh_name_dir + geom_elem_per_dimension->get_name_med() + "/" + _node_or_elem_group_fam;
    
    std::vector< TYPE_FOR_INT_DATASET > fam_map;
      
    dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(file_id, fam_map, fam_name_dir_i);
    



    // ****************** Volume *******************************************
    if(i == (mesh.GetDimension() - 1)) {     //volume
      
      //I have to split the groups with the dimensions
      //I have to compute the number of elements of each group

      for(unsigned gv = 0; gv < group_info.size(); gv++) {

        if(i > 0 && group_info[gv]._node_or_cell_group == group_elements)  {
          if(  i == group_info[gv]._geom_el->get_dimension() - 1) {
            for(/*unsigned*/ int g = 0; g < fam_map.size()/*group_info[gv]._size*//*number_of_group_elements*/; g++) {
              if(fam_map[g] == group_info[gv]._med_flag)   {
                mesh.el->SetElementGroup(g,  group_info[gv]._user_defined_flag /*fam_map[g]*/ /*gr_integer_name*/);  //I think that 1 is set later as the default  group number
                mesh.el->SetElementMaterial(g, group_info[gv]._user_defined_property);
              }

            }
          }  //groups of the current dimension
        }


      }  //loop over all groups

      mesh.el->SetElementGroupNumber(1/*n_groups_of_that_space_dimension*/);

    }
    // ****************** Volume, end *******************************************


// // //     // ****************** Boundary *******************************************
// // //     else   if(i == (mesh.GetDimension() - 1 - 1)) {    //boundary
// // // 
// // //       //loop over volume elements
// // //       //extract faces
// // // 
// // //       //   // read boundary **************** D
// // //       for(unsigned k = 0; k < group_info.size()/*nbcd*//*@todo these should be the groups that are "boundary groups"*/; k++) {  //
// // //         int value = group_info[k]._user_defined_flag;       //flag of the boundary portion
// // //         unsigned nface = group_info[k]._size;  //number of elements in a portion of boundary
// // //         value = - (value + 1);  ///@todo these boundary indices need to be NEGATIVE,  so the value in salome must be POSITIVE
// // //         for(unsigned f = 0; f < nface; f++) {
// // //           //       unsigned iel =,   //volume element to which the face belongs
// // //           //       iel--;
// // //           //       unsigned iface; //index of the face in that volume element
// // //           //       iface = MED_IO::MEDToFemusFaceIndex[mesh.el->GetElementType(iel)][iface-1u];
// // //           //       mesh.el->SetFaceElementIndex(iel,iface,value);  //value is (-1) for element faces that are not boundary faces
// // //         }
// // //       }
// // //       //   // end read boundary **************** D
// // // 
// // //     } //end (volume pos - 1)
// // //     // ****************** Boundary end *******************************************


  }
  
  

  // Connectivities in MED files are stored on a per-node basis: first all 1st nodes, then all 2nd nodes, and so on.
  // Instead, in Gambit they are stored on a per-element basis
  void MED_IO::set_elem_connectivity_and_initialize_elem_group(const hid_t&  file_id, const std::string mesh_menu, const unsigned i, const GeomElemBase*  geom_elem_per_dimension, std::vector<bool>& type_elem_flag) {

    Mesh& mesh = GetMesh();

 
    if(i == mesh.GetDimension() - 1) { 
        
        
        
    std::string my_mesh_name_dir = get_element_info_all_dims_H5Group(mesh_menu);  ///@todo here we have to loop

    // NOD ***************************
    std::string conn_name_dir_i = my_mesh_name_dir +  geom_elem_per_dimension->get_name_med() + "/" + elems_connectivity;
    

       // READ CONNECTIVITY MAP
    std::vector< TYPE_FOR_INT_DATASET > conn_map;
      
    dataset_open_and_close_store_in_vector<TYPE_FOR_INT_DATASET>(file_id, conn_map, conn_name_dir_i);
      
      
      // SET NUMBER OF VOLUME ELEMENTS
         const unsigned el_nodes_per_dimension = geom_elem_per_dimension->n_nodes();

    const unsigned n_elems_per_dimension = conn_map.size() / el_nodes_per_dimension;
    std::cout << " Number of elements of dimension " << (i + 1) << " in med file: " <<  n_elems_per_dimension <<  std::endl;
       
      mesh.SetNumberOfElements(n_elems_per_dimension);
      mesh.el = new elem(n_elems_per_dimension, mesh.GetDimension());



      

      for(unsigned iel = 0; iel < n_elems_per_dimension; iel++) {
                  
        unsigned nve = el_nodes_per_dimension;  /// @todo this is only one element type
        
        if(nve == 27) {
          type_elem_flag[0] = type_elem_flag[3] = true;
          mesh.el->AddToElementNumber(1, geom_elems[HEX]);
          mesh.el->SetElementType(iel, HEX);
        }
        else if(nve == 10) {
          type_elem_flag[1] = type_elem_flag[4] = true;
          mesh.el->AddToElementNumber(1, geom_elems[TET]);
          mesh.el->SetElementType(iel, TET);
        }
        else if(nve == 18) {
          type_elem_flag[2] = type_elem_flag[3] = type_elem_flag[4] = true;
          mesh.el->AddToElementNumber(1, geom_elems[WEDGE]);
          mesh.el->SetElementType(iel, WEDGE);
        }
        else if(nve == 9) {
          type_elem_flag[3] = true;
          mesh.el->AddToElementNumber(1, geom_elems[QUAD]);
          mesh.el->SetElementType(iel, QUAD);
        }
        else if(nve == 6 && mesh.GetDimension() == 2) {
          type_elem_flag[4] = true;
          mesh.el->AddToElementNumber(1, geom_elems[TRI]);
          mesh.el->SetElementType(iel, TRI);
        }
        else if(nve == 3 && mesh.GetDimension() == 1) {
          mesh.el->AddToElementNumber(1, geom_elems[LINE]);
          mesh.el->SetElementType(iel, LINE);
        }
        else {
          std::cout << "Error! Invalid element type in reading  File!" << std::endl;
          std::cout << "Error! Use a second order discretization" << std::endl;
          abort();
        }
        
        for(unsigned i = 0; i < nve; i++) {
          unsigned inode = MEDToFemusVertexIndex[mesh.el->GetElementType(iel)][i];
          const unsigned dof_value = conn_map[iel + i * n_elems_per_dimension ] - 1u;
          mesh.el->SetElementDofIndex(iel, inode, dof_value);  //MED connectivity is stored on a per-node basis, not a per-element basis
        }
        
      }

      
      
      // Initialize Element Group - BEGIN
      for(unsigned iel = 0; iel < n_elems_per_dimension; iel++) {
        mesh.el->SetElementGroup(iel, 1);
      }
      // Initialize Element Group - END
      
      
    }


  }



  //the node global ordering is given by the mesh file, as well as the element global ordering
  void MED_IO::set_node_coordinates(const hid_t&  file_id, const std::string mesh_menu, std::vector < std::vector < double> > & coords, const double Lref) {

    Mesh& mesh = GetMesh();

    std::string coord_dataset = get_node_info_H5Group(mesh_menu) + nodes_coord_list + "/";  ///@todo here we have to loop

    
    std::vector< TYPE_FOR_REAL_DATASET > xyz_med;
      
    dataset_open_and_close_store_in_vector< TYPE_FOR_REAL_DATASET >(file_id, xyz_med, coord_dataset);

    
    
    // SET NUMBER OF NODES
    unsigned int n_nodes = xyz_med.size() / 3; //mesh.GetDimension();
     std::cout << " Number of nodes in med file " <<  n_nodes << " " <<  std::endl;

     mesh.SetNumberOfNodes(n_nodes);

    

    // SET COORDINATES
     for(unsigned d = 0; d < 3; d++)      coords[d].resize(n_nodes);

   for(unsigned d = 0; d < 3; d++)  {
      for(unsigned j = 0; j < n_nodes; j++) {
        coords[d][j] = xyz_med[j + d * n_nodes] / Lref;
      }
    }


  }



  //salome family; our name; our property; group size
  /// @todo check the underscores according to our naming standard
  const GroupInfo MED_IO::get_group_flags_per_mesh(const std::string & group_name, const std::string geom_elem_type) const {

    GroupInfo  group_info;

    const int initial_str_pos = 0;
    std::pair<int, std::vector< int > > gr_family_in_salome_pair = isolate_number_in_string_between_underscores(group_name, initial_str_pos);
    std::pair<int, std::vector< int > > gr_name_pair             = isolate_number_in_string_between_underscores(group_name, gr_family_in_salome_pair.second[1] + 1);
    std::pair<int, std::vector< int > > gr_property_pair         = isolate_number_in_string_between_underscores(group_name, gr_name_pair.second[1] + 1);

    group_info._med_flag              = gr_family_in_salome_pair.first;
    group_info._user_defined_flag     = gr_name_pair.first;
    group_info._user_defined_property = gr_property_pair.first;

    const unsigned pos_group_string_init = gr_family_in_salome_pair.second[1] + 1;
    const unsigned pos_group_string_length = gr_name_pair.second[0] - pos_group_string_init;

    group_info._user_defined_group_string = group_name.substr(pos_group_string_init, pos_group_string_length);

    group_info._node_or_cell_group = geom_elem_type;
    
    if (group_info._node_or_cell_group == group_nodes && group_info._med_flag <= 0) { std::cout << "Node flags internal to Salome are positive" << std::endl; abort(); }

    return  group_info;

  }


  // ************** Groups of each Mesh *********************************
  const std::vector< GroupInfo > MED_IO::get_all_groups_per_mesh(const hid_t&  file_id, const std::string & mesh_menu) const {

    const Mesh & mesh = GetMesh();

    std::vector< std::string >             group_names;
    std::vector< GroupInfo >                group_info;
    
    std::vector< std::string > geom_elem_types = {group_elements, group_nodes};

    for(unsigned geom_elem_type = 0; geom_elem_type < geom_elem_types.size(); geom_elem_type++) {
    
    std::string group_list = get_group_info_H5Group(mesh_menu, geom_elem_types[geom_elem_type]);

    htri_t does_group_exist = H5Lexists(file_id, group_list.c_str(), H5P_DEFAULT);

    if (does_group_exist) {
    
    hid_t  gid      = H5Gopen(file_id, group_list.c_str(), H5P_DEFAULT);
    
    hsize_t    n_groups = get_H5G_size(gid);


    for(unsigned j = 0; j < n_groups; j++) {
        

      std::string group_names_string =  get_H5L_name_by_idx(gid, ".", j);

      group_names.push_back(group_names_string);

      group_info.push_back( get_group_flags_per_mesh(group_names.back(), geom_elem_types[geom_elem_type]) );

    }

    H5Gclose(gid);
    
    }

    }
    
    //======================= After filling all info, Geom elem and size of each group ========================
        for(unsigned i = 0; i < group_info.size(); i++) {
          node_or_elem_Compute_group_geometric_object_type_and_size(file_id, mesh_menu, group_info[i]);
        }
    
    
    return group_info;

  }


   hsize_t   MED_IO::get_H5G_size(const hid_t&  gid) const {

    H5G_info_t  g_info; 
    herr_t status_g = H5Gget_info( gid, & g_info);
    
    if(status_g < 0) {
      std::cout << "Error";
      abort();
    }
    
   return g_info.nlinks;

   }
   

   std::string   MED_IO::get_H5L_name_by_idx(const hid_t&  loc_id, const char *group_name, const unsigned j) const {
       
      char*   group_names_char = new char[_max_length_med_folder];
      ssize_t str_size = H5Lget_name_by_idx(loc_id, group_name/*"."*/, H5_INDEX_NAME, H5_ITER_INC, j, group_names_char, _max_length_med_folder, H5P_DEFAULT);
   
      const std::string link_name(group_names_char);
            
      delete[] group_names_char;

      assert( str_size == link_name.size() );
      
      return link_name;
      
   }
   
  
  // compute number of Mesh fields in Salome file ==============
  const std::vector<std::string> MED_IO::get_mesh_names(const hid_t&  file_id) const {

    hid_t  gid = H5Gopen(file_id, mesh_ensemble.c_str(), H5P_DEFAULT);
    
    hsize_t    n_meshes_ens = get_H5G_size(gid);

    std::vector<std::string>  mesh_menus;

    unsigned n_meshes = 0;

    for(unsigned j = 0; j < n_meshes_ens; j++) {

      std::string tempj =  get_H5L_name_by_idx(gid, ".", j);


      if(tempj.substr(0, 4).compare("Mesh") == 0) {
        n_meshes++;
        mesh_menus.push_back(tempj);
      }
      else {
        std::cout << "Mesh MED fields must start with the word Mesh" << std::endl;
        abort();
      }

    }
    if(n_meshes != n_meshes_ens) {
      std::cout << "Meshes are called with Mesh";
      abort();
    }

    H5Gclose(gid);
    
    
    if(mesh_menus.size() > 1) {
      std::cout << "Review the code because there is only one Mesh object and most likely things don't work" << std::endl;
      abort();
    }
    

    return mesh_menus;
  }
  

  ///@deprecated
  std::string  MED_IO::isolate_first_field_before_underscore(const std::string &  string_in, const int begin_pos_to_investigate) const {

    int str_pos = begin_pos_to_investigate;

    std::string temp_buffer;

    assert(str_pos < string_in.size());

    temp_buffer = string_in.at(str_pos);

    if(temp_buffer.compare("_") == 0)  {
      std::cout <<  "I don't want to start with an underscore" << std::endl;
      abort();
    }

    std::string first_field = "";
    //begin search for the 1st underscore -------------------------------
    while(temp_buffer.compare("_") != 0  && (str_pos < (string_in.size() - 1))) {
      str_pos++;
      first_field += temp_buffer;
      temp_buffer = string_in.at(str_pos);
    }

    return first_field;

  }


  // This function starts from a given point in a string,
  // finds the first two occurrences of underscores,
  // and gets the string in between them
  // If it finds only one underscore and it gets to end-of-file, I want to get that string there
  std::pair<int, std::vector<int> >  MED_IO::isolate_number_in_string_between_underscores(const std::string &  string_in, const int begin_pos_to_investigate) const {

    try {

      int str_pos = begin_pos_to_investigate;
      if(_print_info) {
        std::cout << "Start searching in string " << string_in << " from the character " << string_in.at(str_pos) << " in position " << str_pos <<  std::endl;
      }

      std::vector<int> two_adj_underscores_pos(2, 0);

      std::string temp_buffer;

      assert(str_pos < string_in.size());

      temp_buffer = string_in.at(str_pos);

      if(temp_buffer.compare("_") == 0)  {
        std::cout <<  "I don't want to start with an underscore" << std::endl;
        abort();
      }

      //begin search for the 1st underscore -------------------------------
      while(temp_buffer.compare("_") != 0  && (str_pos < (string_in.size() - 1))) {
        str_pos++;
        temp_buffer = string_in.at(str_pos);
      }

      two_adj_underscores_pos[0] = str_pos;

      //begin search for the 2nd underscore -------------------------------
      if(str_pos < (string_in.size() - 1)) {
        str_pos++;
        temp_buffer = string_in.at(str_pos);
        while(temp_buffer.compare("_") != 0  && (str_pos < (string_in.size() - 1))) {
          str_pos++;
          temp_buffer = string_in.at(str_pos);
        }
        two_adj_underscores_pos[1] = str_pos;

      }
      else if(str_pos == (string_in.size() - 1)) {     //if it reaches the end, it does so after the 1st iteration, because there is an EVEN number of underscores
        two_adj_underscores_pos[0] = begin_pos_to_investigate - 1; //this allows for numbers with more than 1 digit
        two_adj_underscores_pos[1] = str_pos + 1;
      }
      //end search for the 2 underscores -------------------------------

      std::pair<int, int>  delimiting_positions(two_adj_underscores_pos[0], two_adj_underscores_pos[1]);

      int string_to_extract_pos    = delimiting_positions.first + 1;
      int string_to_extract_length = delimiting_positions.second - delimiting_positions.first - 1;

      if(_print_info) {
        std::cout <<  string_to_extract_pos << " " << string_to_extract_length << " " << string_in.substr(string_to_extract_pos, string_to_extract_length).c_str() << " " << std::endl;
      }

      const int flag = atoi(string_in.substr(string_to_extract_pos, string_to_extract_length).c_str());

      std::pair<int, std::vector<int> /*int*/>  flag_and_flag_pos_in_string(flag, two_adj_underscores_pos/*[1]*/);

      return flag_and_flag_pos_in_string;

    }

    catch(const std::out_of_range& e) {
      std::cerr <<  "Reading out of range" << std::endl;
      abort();
    }


  }




  const std::vector< GeomElemBase* > MED_IO::get_geom_elem_type_per_dimension(
    const hid_t & file_id,
    const std::string  my_mesh_name_dir
  ) const {

    const Mesh & mesh = GetMesh();

    const unsigned int dim = mesh.GetDimension();
    std::cout << "No hybrid mesh for now: only 1 FE type per dimension" << std::endl;

    
    hid_t       gid = H5Gopen(file_id, my_mesh_name_dir.c_str(), H5P_DEFAULT);
    
    hsize_t    n_geom_el_types = get_H5G_size(gid);
    
    if (n_geom_el_types != dim) {     std::cout << "Hybrid mesh: implement it" << std::endl; abort();  }
    

    const uint fe_name_nchars = 4;

    std::vector< GeomElemBase* >  geom_elem_per_dimension( dim );

    for(int i = 0; i < (int) dim; i++) {

      std::string temp_i = get_H5L_name_by_idx(file_id, my_mesh_name_dir.c_str(), i);
      


      if( dim  == 3) {

        if( /*temp_i.compare("HE8") == 0 ||*/
          /*temp_i.compare("H20") == 0 ||*/
          temp_i.compare("H27") == 0)  geom_elem_per_dimension[ dim  - 1] = new GeomElemHex27();

        if( /*temp_i.compare("TE4") == 0 ||*/
          temp_i.compare("T10") == 0)  geom_elem_per_dimension[ dim  - 1] = new GeomElemTet10();

        if(/*temp_i.compare("QU4") == 0 ||*/
          /*temp_i.compare("QU8") == 0 ||*/
          temp_i.compare("QU9") == 0) geom_elem_per_dimension[ dim  - 1 - 1] = new GeomElemQuad9();
        if(/*temp_i.compare("TR3") == 0 ||*/
          temp_i.compare("TR6") == 0)  geom_elem_per_dimension[ dim  - 1 - 1] = new GeomElemTri6();

        if(/*temp_i.compare("SE2") == 0 ||*/
          temp_i.compare("SE3") == 0)  geom_elem_per_dimension[ dim  - 1 - 1 - 1] = new GeomElemEdge3();

      }

      else if( dim  == 2) {

        if( /*temp_i.compare("QU4") == 0 ||*/
          /*temp_i.compare("QU8") == 0 ||*/
          temp_i.compare("QU9") == 0)   geom_elem_per_dimension[ dim  - 1] = new GeomElemQuad9();
        if( /*temp_i.compare("TR3") == 0 ||*/
          temp_i.compare("TR6") == 0)   geom_elem_per_dimension[ dim  - 1] = new GeomElemTri6();

        if(/*temp_i.compare("SE2") == 0 ||*/
          temp_i.compare("SE3") == 0)   geom_elem_per_dimension[ dim  - 1 - 1] = new GeomElemEdge3();

      }

      else if( dim  == 1) {

        if( /*temp_i.compare("SE2") == 0 ||*/
          temp_i.compare("SE3") == 0) geom_elem_per_dimension[ dim  - 1] = new GeomElemEdge3();
      }

    }

    H5Gclose(gid);
    
    return geom_elem_per_dimension;
    
  }
  

  const unsigned  MED_IO::get_mesh_dimension_from_attributes(const hid_t &  file_id, const std::string  &  my_mesh_name_dir) const {
    
    const  std::string all_path = mesh_ensemble + "/" + my_mesh_name_dir;
    
    hid_t       gid = H5Gopen(file_id, all_path.c_str(), H5P_DEFAULT);

//     H5O_info_t * 	oinfo;
// H5Oget_info(gid,
//  oinfo
//  /*,H5O_INFO_ALL*/
// ); 


    
// ===============     attr - BEGIN
   int attr_value;
     
   const std::string attr_name = "DIM";
   // const std::string attr_name = "ESP";
 
    hid_t attr = H5Aopen(gid, attr_name.c_str(), H5P_DEFAULT); 

        if ( ( attr ) == H5I_INVALID_HID ) {   abort();    }
   

    const herr_t 	status = H5Aread(attr, H5T_NATIVE_INT, & attr_value);

        if ( status < 0 ) { abort(); }
        
        std::cout << attr_value << "--------------------";
        H5Aclose(attr);
// ===============     attr - END
        
    H5Gclose(gid);
    
    return attr_value;
        
  }
  
  
  // figures out the Mesh dimension by looping over element types
  /// @todo this determination of the dimension from the mesh file would not work with a 2D mesh embedded in 3D
  /// @todo I think I can fix that because I found out about H5 ATTRIBUTES!
  
  const unsigned  MED_IO::get_mesh_dimension_by_looping_over_geom_elem_types(const hid_t &  file_id, const std::string  &  elem_list_in) const  {
    
    
    hid_t       gid = H5Gopen(file_id, elem_list_in.c_str(), H5P_DEFAULT);
    hsize_t    n_geom_elem_types = get_H5G_size(gid);

    unsigned int mydim = 1;  //this is the initial value, then it will be updated below
    
    //this is basically the MANIFOLD DIMENSION of the domain    

    unsigned int dim_aux = mydim;

    std::vector< std::string > elem_types( n_geom_elem_types );


    for(unsigned j = 0; j < elem_types.size(); j++) {
        
      std::string elem_types_str =  get_H5L_name_by_idx(gid, ".", j);
      elem_types[j] = elem_types_str;

      if( /*elem_types_str.compare("HE8") == 0 ||*/
        /*elem_types_str.compare("H20") == 0 ||*/
        elem_types_str.compare("H27") == 0 ||
        /*elem_types_str.compare("TE4") == 0 ||*/
        elem_types_str.compare("T10") == 0)       mydim = 3;
      else if(/*elem_types_str.compare("QU4") == 0 ||*/
        /*elem_types_str.compare("QU8") == 0 ||*/
        elem_types_str.compare("QU9") == 0 ||
        /*elem_types_str.compare("TR3") == 0 ||*/
        elem_types_str.compare("TR6") == 0)   mydim = 2;
      else if(/*elem_types_str.compare("SE2") == 0 ||*/
        elem_types_str.compare("SE3") == 0)  mydim = 1;

      if (mydim > dim_aux) { dim_aux = mydim; }
      
    }  //end for

    H5Gclose(gid);

   
    return dim_aux;
    
  }
  


  GeomElemBase * MED_IO::get_geom_elem_from_med_name(const  std::string el_type) const {

    if(el_type.compare("HE8") == 0) return new GeomElemHex8();
    else if(el_type.compare("H20") == 0) abort(); ///@todo //return new FEHex20();
    else if(el_type.compare("H27") == 0) return new GeomElemHex27();

    else if(el_type.compare("TE4") == 0) return new GeomElemTet4();
    else if(el_type.compare("T10") == 0) return new GeomElemTet10();

    else if(el_type.compare("QU4") == 0) return new GeomElemQuad4();
    else if(el_type.compare("QU8") == 0) abort();
    else if(el_type.compare("QU9") == 0) return new GeomElemQuad9();

    else if(el_type.compare("TR3") == 0) return new GeomElemTri3();
    else if(el_type.compare("TR6") == 0) return new GeomElemTri6();
    else if(el_type.compare("TR7") == 0) abort();

    else if(el_type.compare("SE2") == 0) return new GeomElemEdge2();
    else if(el_type.compare("SE3") == 0) return new GeomElemEdge3();

    else if(el_type.compare(nodes_folder) == 0) return NULL;  //no basic Node class, yet
    else {
      std::cout << "MED_IO::read: element not supported";
      abort();
    }


  }

  

  

} //end namespace femus

