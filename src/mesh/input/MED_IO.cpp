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

//C++ include
#include <cassert>
#include <cstdio>
#include <fstream>
#include <tuple>


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



namespace femus {


  const std::string MED_IO::mesh_ensemble           = "ENS_MAA";
  const std::string MED_IO::aux_zeroone             = "-0000000000000000001-0000000000000000001";
  const std::string MED_IO::elem_list               = "MAI";
  const std::string MED_IO::group_fam               = "FAM";  //both for Elements, and for Nodes
  const std::string MED_IO::connectivity            = "NOD";
  const std::string MED_IO::node_or_elem_global_num = "NUM";  //we don't need to read these fields
  const std::string MED_IO::node_list               = "NOE";
  const std::string MED_IO::coord_list              = "COO";
  const std::string MED_IO::group_ensemble          = "FAS";
  const std::string MED_IO::group_elements          = "ELEME";
  const std::string MED_IO::group_nodes             = "NOEUD";
  const uint MED_IO::max_length = 100;  ///@todo this length of the menu string is conservative enough...



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
    {0, 1, 2, 3, 4, 5},          //TRI6
    {0, 1, 2}                    //EDGE3
  };


  const unsigned MED_IO::MEDToFemusFaceIndex[N_GEOM_ELS][MAX_EL_N_FACES] = {
    {0, 4, 2, 5, 3, 1}, //HEX27
    {0, 1, 2, 3},       //TET10
    {2, 1, 0, 4, 3},    //WEDGE18
    {0, 1, 2, 3},       //QUAD9
    {0, 1, 2},          //TRI6
    {0, 1}              //EDGE3
  };



  /// @todo extend to Wegdes (aka Prisms)
  /// @todo why pass coords other than get it through the Mesh class pointer?
  void MED_IO::read(const std::string& name, vector < vector < double> >& coords, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups_flag, const bool read_boundary_groups_flag) {

    _print_info = true;  
      

    Mesh& mesh = GetMesh();
    mesh.SetLevel(0);

    hid_t  file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    const std::vector< std::string > mesh_menus = get_mesh_names(file_id);

    if(mesh_menus.size() > 1) {
      std::cout << "Review the code because there is only one MultilevelMesh object and most likely things don't work" << std::endl;
      abort();
    }

    // dimension and geom_el types ===============

    const std::vector< GeomElemBase* > geom_elem_per_dimension = set_mesh_dimension_and_get_geom_elems_by_looping_over_element_types(file_id, mesh_menus[0]);

    const unsigned mesh_dim = mesh.GetDimension();

    // meshes ========================
    for(unsigned j = 0; j < mesh_menus.size(); j++) {

      // node coordinates
      set_node_coordinates(file_id, mesh_menus[j], coords, Lref);

      // Volume connectivity
      //       for(unsigned i = 0; i < mesh_dim; i++) {
      unsigned i = mesh_dim - 1;
      set_elem_connectivity(file_id, mesh_menus[j], i, geom_elem_per_dimension[i], type_elem_flag);  //type_elem_flag is to say "There exists at least one element of that type in the mesh"



      if(read_groups_flag == true || read_boundary_groups_flag == true)  {

        // Groups of the mesh ===============
        std::vector< GroupInfo >     group_info = get_group_flags_per_mesh_vector(file_id, mesh_menus[j]);

        for(unsigned i = 0; i < group_info.size(); i++) {
          compute_group_geom_elem_and_size(file_id, mesh_menus[j], group_info[i]);
        }

        // Groups ===============
        if(read_groups_flag == true) {
          for(unsigned i = 0; i < mesh_dim; i++) {
            set_elem_group_ownership(file_id, mesh_menus[j], group_info, geom_elem_per_dimension[i], i);
          }
        }

        // Boundary groups ===============
        if(read_boundary_groups_flag == true)  {
          set_elem_group_ownership_boundary(file_id, mesh_menus[j], group_info, geom_elem_per_dimension, mesh_dim, mesh.el->GetElementNearFaceArray());
        }

      }


    }


    H5Fclose(file_id);

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





  unsigned int MED_IO::get_user_flag_from_med_flag(const std::vector< GroupInfo > & group_info, const TYPE_FOR_FAM_FLAGS med_flag_in) const {

    unsigned int user_flag;

    for(unsigned gv = 0; gv < group_info.size(); gv++) {
      if(group_info[gv]._med_flag == med_flag_in) user_flag = group_info[gv]._user_defined_flag;
    }

    return user_flag;

  }



  //here I need a routine to compute the group GeomElem and the group size

  //separate groups by dimension
  //as soon as an entry is equal to the group_med_flag, that means the dimension is that of the current element dataset
  void MED_IO::compute_group_geom_elem_and_size(const hid_t&  file_id, const std::string mesh_menu, GroupInfo & group_info) const {

    std::string my_mesh_name_dir = mesh_ensemble +  "/" + mesh_menu + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop

    hsize_t     n_geom_el_types;
    hid_t       gid = H5Gopen(file_id, my_mesh_name_dir.c_str(), H5P_DEFAULT);
    hid_t status = H5Gget_num_objs(gid, &n_geom_el_types);

    std::vector<char*> elem_types(n_geom_el_types);

    bool group_found = false;

    //loop over all FAM fields until the group is found
    unsigned j = 0;
    while(j < elem_types.size() && group_found == false) {

      hsize_t dims_fam[2];
      elem_types[j] = new char[max_length];
      H5Gget_objname_by_idx(gid, j, elem_types[j], max_length); ///@deprecated see the HDF doc to replace this
      std::string elem_types_str(elem_types[j]);

      std::string fam_name_dir_i = my_mesh_name_dir + elem_types_str + "/" + group_fam;
      hid_t dtset_fam            = H5Dopen(file_id, fam_name_dir_i.c_str(), H5P_DEFAULT);
      hid_t filespace_fam        = H5Dget_space(dtset_fam);
      hid_t status_fam           = H5Sget_simple_extent_dims(filespace_fam, dims_fam, NULL);
      if(status_fam == 0) {
        std::cerr << "MED_IO::read dims not found";
        abort();
      }

      const unsigned n_elements = dims_fam[0];
      std::vector< TYPE_FOR_FAM_FLAGS > fam_map(n_elements);
      hid_t status_conn = H5Dread(dtset_fam, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fam_map.data());

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

      H5Dclose(dtset_fam);

      j++;
    }

    H5Gclose(gid);

    return;

  }



  void MED_IO::find_boundary_faces_and_set_face_flags(const hid_t&  file_id,
                                                      const std::string mesh_menu,
                                                      const std::vector< GroupInfo > & group_info,
                                                      const GeomElemBase* geom_elem_per_dimension,
                                                      MyMatrix <int> & element_faces_array)  {
    //basically you have certain elements on the boundary, and you have to find them in the list of all elements


    std::string my_mesh_name_dir = mesh_ensemble +  "/" + mesh_menu + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop

    //open the NOD field of the boundary element list (for the connectivities)
    hsize_t dims_conn[2];
    std::string   conn_name_dir = my_mesh_name_dir + geom_elem_per_dimension->get_name_med() + "/" + connectivity; ///@todo these boundary connectivities were not stored, so we need to read them now
    hid_t dtset_conn     = H5Dopen(file_id, conn_name_dir.c_str(), H5P_DEFAULT);
    hid_t filespace_conn = H5Dget_space(dtset_conn);
    hid_t status_conn    = H5Sget_simple_extent_dims(filespace_conn, dims_conn, NULL);

    std::vector<int> conn_map(dims_conn[0]);
    hid_t status2_conn = H5Dread(dtset_conn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, conn_map.data());
    H5Dclose(dtset_conn);

    //open the FAM field of the boundary element list (for the group flags)
    hsize_t dims_fam[2];
    std::string fam_name_dir = my_mesh_name_dir + geom_elem_per_dimension->get_name_med() + "/" + group_fam;
    hid_t dtset_fam     = H5Dopen(file_id, fam_name_dir.c_str(), H5P_DEFAULT);
    hid_t filespace_fam = H5Dget_space(dtset_fam);
    hid_t status_fam    = H5Sget_simple_extent_dims(filespace_fam, dims_fam, NULL);

    std::vector< TYPE_FOR_FAM_FLAGS > fam_map(dims_fam[0]);
    hid_t status2_fam = H5Dread(dtset_fam, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fam_map.data());
    H5Dclose(dtset_fam);

    //check that all boundary faces were set in the mesh file
    for(unsigned i = 0; i < fam_map.size(); i++) {
      if(fam_map[i] == 0) {
        std::cout << "Some boundary face was not set in the mesh MED file" << std::endl;
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

        unsigned n_nodes = _geom_elems[iel_geom_type]->get_face(f).size();

        // on one hand I construct the boundary face connectivity from the volume connectivity, with the order that was given in our code
        std::vector<unsigned> face_nodes_from_vol_connectivity(n_nodes);

        for(unsigned nd = 0; nd < n_nodes; nd++) {
          unsigned nd_of_face = _geom_elems[iel_geom_type]->get_face(f)[nd];
          face_nodes_from_vol_connectivity[nd] = mesh.el->GetElementDofIndex(iel, nd_of_face);

        }

        // on the other I read the boundary face connectivity from the list of boundary faces, which is the one that has the FAM information
        //loop over all the bdry group elements

        std::vector<unsigned> face_nodes_from_bdry_group(n_nodes);

        for(unsigned k = 0; k < fam_map.size(); k++) {

          for(unsigned nd = 0; nd < n_nodes; nd++) {
            face_nodes_from_bdry_group[nd] = conn_map[ k + nd * fam_map.size() ] - 1;
          }

          bool are_faces_the_same = see_if_faces_from_different_lists_are_the_same(geom_elem_per_dimension, face_nodes_from_vol_connectivity, face_nodes_from_bdry_group);

                 const TYPE_FOR_FAM_FLAGS med_flag = fam_map[k];
                 
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

    std::string my_mesh_name_dir = mesh_ensemble +  "/" + mesh_menu + "/" +  aux_zeroone + "/" + node_list/*elem_list*/ + "/";  ///@todo here we have to loop

    //open the FAM field of NOE
    hsize_t dims_fam[2];
    std::string fam_name_dir = my_mesh_name_dir /*+ geom_elem_per_dimension->get_name_med()*/ + "/" + group_fam;
    hid_t dtset_fam     = H5Dopen(file_id, fam_name_dir.c_str(), H5P_DEFAULT);
    hid_t filespace_fam = H5Dget_space(dtset_fam);
    hid_t status_fam    = H5Sget_simple_extent_dims(filespace_fam, dims_fam, NULL);

    std::vector< TYPE_FOR_FAM_FLAGS > fam_map(dims_fam[0]);
    hid_t status2_fam = H5Dread(dtset_fam, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fam_map.data());
    H5Dclose(dtset_fam);


    // I have to loop over elements to then find faces of the elements, get the dof of those faces, and then with that dof go get the position in NOE/FAM 
       
       
            // loop over elements to find faces
       Mesh& mesh = GetMesh();
       
       //loop over the volume connectivity and find the boundary faces 
       //the boundary faces of each volume element have already been constructed after the mesh reading
       
        for(unsigned iel = 0; iel < mesh.GetNumberOfElements(); iel++) {
            
//             unsigned iel_geom_type = mesh.GetElementType(iel);
                               
            for(unsigned f = 0; f < mesh.GetElementFaceNumber(iel); f++) {
                
                unsigned n_nodes_in_face = 1 /*_geom_elems[iel_geom_type]->get_face(f).size()*/;  //here the faces are made of 1 node
                
                std::vector<unsigned> face_nodes(n_nodes_in_face);
                
               for(unsigned nd = 0; nd < face_nodes.size(); nd++) {
                   unsigned nd_of_face = f /*_geom_elems[iel_geom_type]->get_face(f)[nd]*/;
                   face_nodes[nd] = mesh.el->GetElementDofIndex(iel, nd_of_face);
                
               }
            
            // now I take this dof and read from NOE/FAM  
            // If the group is different from zero 
               for(unsigned k = 0; k < fam_map.size(); k++) {

                      const TYPE_FOR_FAM_FLAGS med_flag = fam_map[k];
                      
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
                                        
                  
// =======
//     // I have to loop over elements to then find faces of the elements, get the dof of those faces, and then with that dof go get the position in NOE/FAM
// 
// 
//     // loop over elements to find faces
//     Mesh& mesh = GetMesh();
// 
//     //loop over the volume connectivity and find the boundary faces
//     //the boundary faces of each volume element have already been constructed after the mesh reading
// 
//     for(unsigned iel = 0; iel < mesh.GetNumberOfElements(); iel++) {
// 
//       //             unsigned iel_geom_type = mesh.GetElementType(iel);
// 
//       for(unsigned f = 0; f < mesh.GetElementFaceNumber(iel); f++) {
// 
//         unsigned n_nodes = 1 /*_geom_elems[iel_geom_type]->get_face(f).size()*/;  //here the faces are made of 1 node
// 
//         std::vector<unsigned> face_nodes(n_nodes);
// 
//         for(unsigned nd = 0; nd < n_nodes; nd++) {
//           unsigned nd_of_face = 0 /*_geom_elems[iel_geom_type]->get_face(f)[nd]*/;
//           face_nodes[nd] = mesh.el->GetElementDofIndex(iel, nd_of_face);
// 
//         }
// 
//         // now I take this dof and read from NOE/FAM
//         // If the group is different from zero
//         for(unsigned k = 0; k < fam_map.size(); k++) {
// 
//           const TYPE_FOR_FAM_FLAGS med_flag = fam_map[k];
// 
//           if(med_flag != 0)  {
// 
//             int user_flag =  get_user_flag_from_med_flag(group_info, med_flag);  //flag of the boundary portion
//             user_flag = - (user_flag + 1);  ///@todo these boundary indices need to be NEGATIVE,  so the user_flag in salome must be POSITIVE
// 
//             std::cout << "Found face " << k << " in element " << iel << " with MED flag " << med_flag << " and user flag " << user_flag << std::endl;
// 
//             //       unsigned iface = MED_IO::MEDToFemusFaceIndex[mesh.el->GetElementType(iel)][iface-1u];//index of the face in that volume element
//             element_faces_array[iel][f] = user_flag;  //user_flag is (-1) for element faces that are not boundary faces, SO WE MUST BE CAREFUL HERE!
//             //  mesh.el->SetFaceElementIndex(iel, f, user_flag);  //old version
// 
// >>>>>>> main/master_gcc_7
//           }
// 
// 
//         } //end k
// 
//       } //faces loop
// 
//     }
// 
// 
//   }
// 
// 
// 
//   bool MED_IO::see_if_faces_from_different_lists_are_the_same(const GeomElemBase* geom_elem_per_dimension,
//                                                               const std::vector< unsigned > & face_nodes_from_vol_connectivity,
//                                                               const std::vector< unsigned > & face_nodes_from_bdry_group) {
// 
//     // check any possible order of faces
// 
//     //just look for the initial linear element (maybe even the 1st three only); if this is aligned, all the Quad9 will be aligned
//     //The problem, is that we don't know if the order corresponds to the OUTWARD NORMAL or not.
//     //How many ways are there? If the nodes are 0 1 2 3, it could only be 0123, or 1230, or 2301, or 3012, or the REVERSE of each of them
//     std::vector<unsigned>  face_nodes_from_vol_connectivity_linear(face_nodes_from_vol_connectivity.begin(), face_nodes_from_vol_connectivity.begin() + geom_elem_per_dimension->n_nodes_linear());
//     std::vector<unsigned>  face_nodes_from_bdry_group_linear(face_nodes_from_bdry_group.begin(), face_nodes_from_bdry_group.begin() + geom_elem_per_dimension->n_nodes_linear());
// 
//     unsigned n_alternatives = 2 * geom_elem_per_dimension->n_nodes_linear();
//     std::vector< std::vector<unsigned> > face_alternatives(n_alternatives);
// 
//     for(unsigned alt = 0; alt < n_alternatives / 2; alt++) {
//       unsigned face_length = geom_elem_per_dimension->n_nodes_linear();
//       face_alternatives[alt].resize(face_length);
//       for(unsigned i = 0; i < face_length; i++) {
//         unsigned mod_index = (i + alt) % face_length;
//         face_alternatives[alt][i] = face_nodes_from_bdry_group_linear[mod_index];
//       }
//     }
// 
//     for(unsigned alt = n_alternatives / 2; alt < n_alternatives; alt++) {
//       face_alternatives[alt] = face_alternatives[alt - (n_alternatives / 2)];
//       std::reverse(face_alternatives[alt].begin(), face_alternatives[alt].end());
//     }
// 
// 
//     std::vector<bool>  is_same_face(n_alternatives);
//     bool bool_union = false;
//     for(unsigned alt = 0; alt < n_alternatives; alt++) {
//       is_same_face[alt] = (face_nodes_from_vol_connectivity_linear == face_alternatives[alt]);
//       if(is_same_face[alt] == true) bool_union = true;
//     }
// 
//     return bool_union;


  }




  //  we loop over all elements and see which ones are of that group
  //
  void MED_IO::set_elem_group_ownership(
    const hid_t&  file_id,
    const std::string mesh_menu,
    const std::vector<GroupInfo> & group_info,
    const GeomElemBase* geom_elem_per_dimension,
    const int i
  )  {

    Mesh& mesh = GetMesh();

    //FAM ***************************
    hsize_t dims_fam[2];
    std::string my_mesh_name_dir = mesh_ensemble +  "/" + mesh_menu + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop
    std::string fam_name_dir_i = my_mesh_name_dir + geom_elem_per_dimension->get_name_med() + "/" + group_fam;
    hid_t dtset_fam = H5Dopen(file_id, fam_name_dir_i.c_str(), H5P_DEFAULT);
    hid_t filespace_fam = H5Dget_space(dtset_fam);
    hid_t status_fam  = H5Sget_simple_extent_dims(filespace_fam, dims_fam, NULL);
    if(status_fam == 0) {
      std::cerr << "MED_IO::read dims not found";
      abort();
    }

    const unsigned n_elements = dims_fam[0];
    std::vector< TYPE_FOR_FAM_FLAGS > fam_map(n_elements);
    hid_t status_conn = H5Dread(dtset_fam, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fam_map.data());


    // ****************** Volume *******************************************
    if(i == (mesh.GetDimension() - 1)) {     //volume
      std::vector < unsigned > materialElementCounter(3, 0); ///@todo I think this counts who is fluid, who is solid, who whatever else, need to double check with the Gambit input files
      const unsigned group_property_fluid_probably          = 2;
      const unsigned group_property_something_else_probably = 3;
      const unsigned group_property_solid_probably          = 4;

      //I have to split the groups with the dimensions
      //I have to compute the number of elements of each group

      for(unsigned gv = 0; gv < group_info.size(); gv++) {

        if(i > 0)  {
          if(i == group_info[gv]._geom_el->get_dimension() - 1) {
            for(/*unsigned*/ int g = 0; g < fam_map.size()/*group_info[gv]._size*//*number_of_group_elements*/; g++) {
              if(fam_map[g] == group_info[gv]._med_flag)   {
                mesh.el->SetElementGroup(g,  group_info[gv]._user_defined_flag /*fam_map[g]*/ /*gr_integer_name*/);  //I think that 1 is set later as the default  group number
                mesh.el->SetElementMaterial(g, group_info[gv]._user_defined_property);
                //         mesh.el->SetElementMaterial(elem_indices[g] - 1 - n_elements_b_bb, group_info[gv]._user_defined_property /*gr_material*/);
              }

              if(group_info[gv]._user_defined_property/*gr_material*/ == group_property_fluid_probably) materialElementCounter[0] += 1;
              else if(group_info[gv]._user_defined_property/*gr_material*/ == group_property_something_else_probably) materialElementCounter[1] += 1;
              else                                                            materialElementCounter[2] += 1;

            }
          }  //groups of the current dimension
        }


      }  //loop over all groups

      mesh.el->SetElementGroupNumber(1/*n_groups_of_that_space_dimension*/);
      mesh.el->SetMaterialElementCounter(materialElementCounter);
    }
    // ****************** Volume, end *******************************************


    // ****************** Boundary *******************************************
    else   if(i == (mesh.GetDimension() - 1 - 1)) {    //boundary

      //loop over volume elements
      //extract faces

      //   // read boundary **************** D
      for(unsigned k = 0; k < group_info.size()/*nbcd*//*@todo these should be the groups that are "boundary groups"*/; k++) {  //
        int value = group_info[k]._user_defined_flag;       //flag of the boundary portion
        unsigned nface = group_info[k]._size;  //number of elements in a portion of boundary
        value = - (value + 1);  ///@todo these boundary indices need to be NEGATIVE,  so the value in salome must be POSITIVE
        for(unsigned f = 0; f < nface; f++) {
          //       unsigned iel =,   //volume element to which the face belongs
          //       iel--;
          //       unsigned iface; //index of the face in that volume element
          //       iface = MED_IO::MEDToFemusFaceIndex[mesh.el->GetElementType(iel)][iface-1u];
          //       mesh.el->SetFaceElementIndex(iel,iface,value);  //value is (-1) for element faces that are not boundary faces
        }
      }
      //   // end read boundary **************** D

    } //end (volume pos - 1)
    // ****************** Boundary end *******************************************

    H5Dclose(dtset_fam);

  }

  // Connectivities in MED files are stored on a per-node basis: first all 1st nodes, then all 2nd nodes, and so on.
  // Instead, in Gambit they are stored on a per-element basis
  void MED_IO::set_elem_connectivity(const hid_t&  file_id, const std::string mesh_menu, const unsigned i, const GeomElemBase*  geom_elem_per_dimension, std::vector<bool>& type_elem_flag) {

    Mesh& mesh = GetMesh();

    const unsigned el_nodes_per_dimension = geom_elem_per_dimension->n_nodes();

    std::string my_mesh_name_dir = mesh_ensemble +  "/" + mesh_menu + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop
    hsize_t dims_i[2];
    // NOD ***************************
    std::string conn_name_dir_i = my_mesh_name_dir +  geom_elem_per_dimension->get_name_med() + "/" + connectivity;
    hid_t dtset_conn = H5Dopen(file_id, conn_name_dir_i.c_str(), H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dtset_conn);
    hid_t status_els_i  = H5Sget_simple_extent_dims(filespace, dims_i, NULL);
    if(status_els_i == 0) {
      std::cerr << "MED_IO::read dims not found";
      abort();
    }


    const int dim_conn = dims_i[0];
    const unsigned n_elems_per_dimension = dim_conn / el_nodes_per_dimension;
    std::cout << " Number of elements of dimension " << (i + 1) << " in med file: " <<  n_elems_per_dimension <<  std::endl;

    // SET NUMBER OF VOLUME ELEMENTS
    if(i == mesh.GetDimension() - 1) {
      mesh.SetNumberOfElements(n_elems_per_dimension);
      mesh.el = new elem(n_elems_per_dimension);    ///@todo check where this is going to be deleted. It is in the Destructor of Mesh

      // READ CONNECTIVITY MAP
      int* conn_map = new  int[dim_conn];
      hid_t status_conn = H5Dread(dtset_conn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, conn_map);
      if(status_conn != 0) {
        std::cout << "MED_IO::read: connectivity not found";
        abort();
      }


      for(unsigned iel = 0; iel < n_elems_per_dimension; iel++) {
        mesh.el->SetElementGroup(iel, 1);
        unsigned nve = el_nodes_per_dimension;  /// @todo this is only one element type
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
          abort();
        }
        for(unsigned i = 0; i < nve; i++) {
          unsigned inode = MEDToFemusVertexIndex[mesh.el->GetElementType(iel)][i];
          const unsigned dof_value = conn_map[iel + i * n_elems_per_dimension ] - 1u;
          mesh.el->SetElementDofIndex(iel, inode, dof_value);  //MED connectivity is stored on a per-node basis, not a per-element basis
        }
      }

      // clean
      delete [] conn_map;

    }

    H5Dclose(dtset_conn);

  }




  void MED_IO::set_node_coordinates(const hid_t&  file_id, const std::string mesh_menu, vector < vector < double> > & coords, const double Lref) {

    Mesh& mesh = GetMesh();
    hsize_t dims[2];

    std::string coord_dataset = mesh_ensemble +  "/" + mesh_menu + "/" +  aux_zeroone + "/" + node_list + "/" + coord_list + "/";  ///@todo here we have to loop

    hid_t dtset = H5Dopen(file_id, coord_dataset.c_str(), H5P_DEFAULT);

    // SET NUMBER OF NODES
    hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
    hid_t status_dims  = H5Sget_simple_extent_dims(filespace, dims, NULL);
    if(status_dims == 0) std::cerr << "MED_IO::read dims not found";
    // reading xyz_med
    unsigned int n_nodes = dims[0] / 3; //mesh.GetDimension();
    double*   xyz_med = new double[dims[0]];
    std::cout << " Number of nodes in med file " <<  n_nodes << " " <<  std::endl;

    mesh.SetNumberOfNodes(n_nodes);

    for(unsigned d = 0; d < 3; d++)      coords[d].resize(n_nodes);

    H5Dread(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz_med);
    H5Dclose(dtset);

    for(unsigned d = 0; d < 3; d++)  {
      for(unsigned j = 0; j < n_nodes; j++) {
        coords[d][j] = xyz_med[j + d * n_nodes] / Lref;
      }
    }

    delete[] xyz_med;

  }



  //salome family; our name; our property; group size
  /// @todo check the underscores according to our naming standard
  const GroupInfo MED_IO::get_group_flags_per_mesh(const std::string & group_name) const {

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

    return  group_info;

  }


  // ************** Groups of each Mesh *********************************
  const std::vector< GroupInfo > MED_IO::get_group_flags_per_mesh_vector(const hid_t&  file_id, const std::string & mesh_menu) const {

    const Mesh & mesh = GetMesh();

    std::string group_list;
    if(mesh.GetDimension() == 2 || mesh.GetDimension() == 3) group_list = group_ensemble +  "/" + mesh_menu + "/" + group_elements;
    else if(mesh.GetDimension() == 1)                        group_list = group_ensemble +  "/" + mesh_menu + "/" + group_nodes;

    hid_t  gid_groups      = H5Gopen(file_id, group_list.c_str(), H5P_DEFAULT);

    hsize_t n_groups = 0;
    hid_t status_groups = H5Gget_num_objs(gid_groups, &n_groups);

    std::vector< std::string >             group_names(n_groups);
    std::vector< GroupInfo >                group_info(n_groups);

    for(unsigned j = 0; j < n_groups; j++) {

      char*   group_names_char = new char[max_length];
      H5Gget_objname_by_idx(gid_groups, j, group_names_char, max_length); ///@deprecated see the HDF doc to replace this
      group_names[j] = group_names_char;
      delete[] group_names_char;

      group_info[j] = get_group_flags_per_mesh(group_names[j]);

    }

    H5Gclose(gid_groups);

    return group_info;

  }


  // compute number of Mesh fields in Salome file ==============
  const std::vector<std::string> MED_IO::get_mesh_names(const hid_t&  file_id) const {

    hid_t  gid = H5Gopen(file_id, mesh_ensemble.c_str(), H5P_DEFAULT);

    hsize_t     n_meshes_ens;
    hid_t status = H5Gget_num_objs(gid, &n_meshes_ens); // number of meshes
    if(status != 0) {
      std::cout << "Number of mesh menus not found";
      abort();
    }

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
  ) {

    Mesh& mesh = GetMesh();

    const int n_fem_types = mesh.GetDimension();
    std::cout << "No hybrid mesh for now: only 1 FE type per dimension" << std::endl;


    // Get the element name
    char** el_fem_type = new char*[n_fem_types];

    const uint fe_name_nchars = 4;

    std::vector< GeomElemBase* >  geom_elem_per_dimension(mesh.GetDimension());

    for(int i = 0; i < (int) n_fem_types; i++) {

      el_fem_type[i] = new char[fe_name_nchars];
      H5Lget_name_by_idx(file_id, my_mesh_name_dir.c_str(), H5_INDEX_NAME, H5_ITER_INC, i, el_fem_type[i], fe_name_nchars, H5P_DEFAULT);
      std::string temp_i(el_fem_type[i]);

      if(mesh.GetDimension() == 3) {

        if( /*temp_i.compare("HE8") == 0 ||*/
          /*temp_i.compare("H20") == 0 ||*/
          temp_i.compare("H27") == 0)  geom_elem_per_dimension[mesh.GetDimension() - 1] = new GeomElemHex27();

        if( /*temp_i.compare("TE4") == 0 ||*/
          temp_i.compare("T10") == 0)  geom_elem_per_dimension[mesh.GetDimension() - 1] = new GeomElemTet10();

        if(/*temp_i.compare("QU4") == 0 ||*/
          /*temp_i.compare("QU8") == 0 ||*/
          temp_i.compare("QU9") == 0) geom_elem_per_dimension[mesh.GetDimension() - 1 - 1] = new GeomElemQuad9();
        if(/*temp_i.compare("TR3") == 0 ||*/
          temp_i.compare("TR6") == 0)  geom_elem_per_dimension[mesh.GetDimension() - 1 - 1] = new GeomElemTri6();

        if(/*temp_i.compare("SE2") == 0 ||*/
          temp_i.compare("SE3") == 0)  geom_elem_per_dimension[mesh.GetDimension() - 1 - 1 - 1] = new GeomElemEdge3();

      }

      else if(mesh.GetDimension() == 2) {

        if( /*temp_i.compare("QU4") == 0 ||*/
          /*temp_i.compare("QU8") == 0 ||*/
          temp_i.compare("QU9") == 0)   geom_elem_per_dimension[mesh.GetDimension() - 1] = new GeomElemQuad9();
        if( /*temp_i.compare("TR3") == 0 ||*/
          temp_i.compare("TR6") == 0)   geom_elem_per_dimension[mesh.GetDimension() - 1] = new GeomElemTri6();

        if(/*temp_i.compare("SE2") == 0 ||*/
          temp_i.compare("SE3") == 0)   geom_elem_per_dimension[mesh.GetDimension() - 1 - 1] = new GeomElemEdge3();

      }

      else if(mesh.GetDimension() == 1) {

        if( /*temp_i.compare("SE2") == 0 ||*/
          temp_i.compare("SE3") == 0) geom_elem_per_dimension[mesh.GetDimension() - 1] = new GeomElemEdge3();
      }

    }

    // clean
    for(int i = 0; i < (int)n_fem_types; i++) delete[] el_fem_type[i];
    delete[] el_fem_type;

    return geom_elem_per_dimension;
  }

  // figures out the Mesh dimension by looping over element types
  /// @todo this determination of the dimension from the mesh file would not work with a 2D mesh embedded in 3D

  const std::vector< GeomElemBase* >  MED_IO::set_mesh_dimension_and_get_geom_elems_by_looping_over_element_types(const hid_t &  file_id, const std::string & mesh_menus)  {


    std::string my_mesh_name_dir = mesh_ensemble +  "/" + mesh_menus + "/" +  aux_zeroone + "/" + elem_list + "/";  ///@todo here we have to loop

    hsize_t     n_fem_type;
    hid_t       gid = H5Gopen(file_id, my_mesh_name_dir.c_str(), H5P_DEFAULT);
    hid_t status = H5Gget_num_objs(gid, &n_fem_type);
    if(status != 0) {
      std::cout << "MED_IO::read_fem_type:   H5Gget_num_objs not found";
      abort();
    }

    Mesh& mesh = GetMesh();
    uint mydim = 1;  //this is the initial value, then it will be updated below
    mesh.SetDimension(mydim);  //this is basically the MANIFOLD DIMENSION of the domain


    std::vector<char*> elem_types(n_fem_type);


    for(unsigned j = 0; j < elem_types.size(); j++) {
      elem_types[j] = new char[max_length];
      H5Gget_objname_by_idx(gid, j, elem_types[j], max_length); ///@deprecated see the HDF doc to replace this
      std::string elem_types_str(elem_types[j]);

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

      if(mydim > mesh.GetDimension()) mesh.SetDimension(mydim);

    }  //end for

    H5Gclose(gid);


    const std::vector< GeomElemBase* >  geom_elem_per_dimension = get_geom_elem_type_per_dimension(file_id, my_mesh_name_dir);

    //        if(mesh.GetDimension() != n_fem_type) { std::cout << "Mismatch between dimension and number of element types" << std::endl;   abort();  }
    // ///@todo removed this check to allow 2d object in 3d

    return geom_elem_per_dimension;
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
    else {
      std::cout << "MED_IO::read: element not supported";
      abort();
    }


  }


} //end namespace femus

