/*=========================================================================

 Program: FEMUS
 Module: MED_IO
 Authors: Giorgio Bornia, Sureka Pathmanathan
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_mesh_MED_IO_hpp__
#define __femus_mesh_MED_IO_hpp__

// Local includes
#include "MeshInput.hpp"
#include "GeomElTypeEnum.hpp"
#include "GeomElemBase.hpp"
#include "GeomElemHex27.hpp"
#include "GeomElemTet10.hpp"
#include "GeomElemQuad9.hpp"
#include "GeomElemTri6.hpp"
#include "GeomElemEdge3.hpp"

#include "FemusConfig.hpp"


#include "MED_IO_Group_Info.hpp"


#ifdef HAVE_HDF5

#include "hdf5.h"



namespace femus
{

   // Forward declarations
class Mesh;


/**
 * This class implements reading meshes in the MED format.
 */

// ------------------------------------------------------------
class MED_IO : public MeshInput< Mesh >
{

// === Constructors / Destructor - BEGIN =================
 public:
    
  /**
   * Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements via the read() command.
   */
  explicit
  MED_IO (Mesh& mesh);
  
  ~MED_IO();
// === Constructors / Destructor - END =================
  
  
// === Debug - BEGIN =================
   
 private:
    
   bool _print_info;
   
// === Debug - END =================
   
   
// === Read, main function - BEGIN =================
 public:
    
  /**
   * Reads in a mesh in the  *.med format
   */
  virtual void read (const std::string& name, 
                     std::vector < std::vector < double> > &coords,
                     const double Lref, 
                     std::vector<bool> &type_elem_flag, 
                     const bool read_groups, 
                     const bool read_boundary_groups);
// === Read, main function - END =================
  
  
// === Mesh, file, H5File - BEGIN =================
 public:
    
  hid_t open_mesh_file(const std::string& name);
  
  void close_mesh_file(hid_t file_id);
// === Mesh, file, H5File - END =================


// === Mesh, file, Attributes of a Group, Dataset or Datatype - BEGIN =================
//           if ((attr = H5Aopen(file, attr_name, H5P_DEFAULT)) == H5I_INVALID_HID) {
//             ret_val = EXIT_FAILURE;
//             goto fail_attr;
//         }
//         // read the attribute value
//         if (H5Aread(attr, H5T_NATIVE_INT, &value) < 0)
//             ret_val = EXIT_FAILURE;
//  
//         // do something w/ the attribute value
//  
//         H5Aclose(attr);

// === Mesh, file, Attributes of a Group, Dataset or Datatype - END =================
  
  
// === Mesh, File, H5Groups (folders in the file) - BEGIN =================
 private:
  
  hsize_t  get_H5G_size(const hid_t&  gid) const;
   
// === Mesh, File, H5Groups (folders in the file) - END =================


// === Mesh, file, H5Groups, H5Links of a Group (subobjects of that Group) - BEGIN =================
 private:
    
  std::string  get_H5L_name_by_idx(const hid_t&  loc_id, const char *group_name, const unsigned j) const;
  
// === Mesh, file, H5Groups, H5Links of a Group (subobjects of that Group) - END =================
  

  
  
// === Mesh, file, H5Datasets (arrays of data in the file) - BEGIN =================
 private:
    
 template < class DATASET_TYPE >  
  void dataset_open_and_close_store_in_vector(hid_t file_id, std::vector< DATASET_TYPE > & fam_map, const std::string fam_name_dir_i) const;
  

// === Mesh, file, H5Datasets (arrays of data in the file) - END =================

  
   
  

// === Mesh, dimension and Geom Elems - BEGIN =================
 private:
    
   /** Determine mesh dimension from mesh file. It cannot be const because it sets the dimension in the mesh */
   const std::vector< GeomElemBase* >  set_mesh_dimension_and_get_geom_elems_by_looping_over_element_types(const hid_t &  file_id, const std::string & menu_name);
// === Mesh, dimension and Geom Elems - END =================


// === Geometric elements, Connectivities - BEGIN =================
 private:
     
   /** Map from Salome vertex index to Femus vertex index */
   static const unsigned MEDToFemusVertexIndex[N_GEOM_ELS][MAX_EL_N_NODES]; 
 
   /** Map from Salome face index to Femus face index */
   static const unsigned MEDToFemusFaceIndex[N_GEOM_ELS][MAX_EL_N_FACES];
// === Geometric elements, Connectivities - END =================


// === Geometric elements, types - BEGIN =================
 private:
    
   GeomElemBase * get_geom_elem_from_med_name(const  std::string el_type) const;

   /** Read FE type */
  const std::vector< GeomElemBase* > get_geom_elem_type_per_dimension(const hid_t & file_id, const std::string my_mesh_name_dir) const;
   
   std::vector< GeomElemBase* > _geom_elems;
// === Geometric elements, types - END =================


// === Mesh, Elements or Nodes or Groups - BEGIN =================
 public:
    
   /// menus (folder of the mesh)
   const std::vector<std::string>  get_mesh_names(const hid_t & file_id) const;
   
// === Mesh, Elements or Nodes or Groups - END =================
   
   

// === Mesh, Elements or Nodes - BEGIN =================
 private:
    
   //elements or nodes
   static const uint _max_length_med_folder;  //max length of a string of an H5Link
   static const std::string mesh_ensemble;             //ENS_MAA
   static const std::string aux_zeroone;               // -0000000000000000001-0000000000000000001

   //elements or nodes, gui global numbering
   static const std::string node_or_elem_salome_gui_global_num;   //NUM    //this is the MED global numbering (as you see in Salome) both for nodes (in NOE) and for elements of all dimensions (in MAI). 
     //numeration in the Salome GUI, we don't need to read these fields
     //      
     //Salome global Numbering of both Nodes and Elements starts at 1.
     //For Elements, lower dimensional elements are numbered first
     
// === Mesh, Elements or Nodes - END =================


// === Mesh, Elements - BEGIN =================
 private:
  
   void set_elem_connectivity(const hid_t&  file_id, const std::string mesh_menu, const unsigned i, const GeomElemBase* geom_elem_per_dimension, std::vector<bool>& type_elem_flag);

   //elements
   std::string get_element_info_all_dims_H5Group(const std::string mesh_menu) const;

   static const std::string elem_types_folder;         //MAI
   static const std::string elems_connectivity;        //NOD    //These are written based on the NOE/NUM numbering !
   
// === Mesh, Elements - END =================

   
// === Mesh, Nodes - BEGIN =================
 private:
  
  void set_node_coordinates(const hid_t&  file_id, const std::string mesh_menu, std::vector < std::vector < double> >& coords, const double Lref);
    
  //nodes
   std::string get_node_info_H5Group(const std::string mesh_menu) const;
   
   static const std::string nodes_folder;                 //NOE
   static const std::string nodes_coord_list;             //COO
// === Mesh, Nodes - END =================
  

  
// === Mesh, Groups - BEGIN =================


// === Mesh, Groups, Elements or Nodes - BEGIN =================
 public:
    
   const std::vector< GroupInfo > get_all_groups_per_mesh(const hid_t &  file_id, const std::string & mesh_menu) const;
    
 private:

   const GroupInfo                get_group_flags_per_mesh(const std::string & group_names, const std::string  geom_elem_type) const;
  
   void node_or_elem_Compute_group_geometric_object_type_and_size(const hid_t&  file_id, const std::string mesh_menu, GroupInfo & group_info)  const;

   
   // groups, elements or nodes
   std::string get_group_info_H5Group(const std::string mesh_menu,  const std::string geom_elem_type) const;
  
   static const std::string _node_or_elem_group_fam;   //FAM
   static const std::string group_ensemble;            //FAS
   // groups, elements
   static const std::string group_elements;            //ELEME
   // groups, nodes
   static const std::string group_nodes;               //NOEUD
   
// === Mesh, Groups, Elements or Nodes - END =================


// === Mesh, Groups, Elements or Nodes, string processing - BEGIN =================
 private:
     
   unsigned int get_user_flag_from_med_flag(const std::vector< GroupInfo > & group_info, const TYPE_FOR_INT_DATASET med_flag_in ) const;
   
   unsigned int get_med_flag_from_user_flag(const std::vector< GroupInfo > & group_info, const TYPE_FOR_INT_DATASET input_flag) const;

   std::pair<int, std::vector<int> >  isolate_number_in_string_between_underscores(const std::string & string_in, const int begin_pos_to_investigate) const;
      
   std::string  isolate_first_field_before_underscore(const std::string &  string_in, const int begin_pos_to_investigate) const;
// === Mesh, Groups, Elements or Nodes, string processing - END =================

   
// === Mesh, Groups, Volume - Elements - BEGIN =================
 private:
  
   void set_elem_group_ownership(const hid_t&  file_id,
                                 const std::string mesh_menu,
                                 const std::vector<GroupInfo> & group_info,
                                 const GeomElemBase* geom_elem_per_dimension,
                                 const int i );
   
// === Mesh, Groups, Volume - Elements - END =================
      

// === Mesh, Groups, Boundary - Elements or Nodes ( Elements(2D/3D) or Nodes(1D) ) - BEGIN =================
 private:
  
    void set_elem_group_ownership_boundary(const hid_t  file_id,
                                           const std::string mesh_menu,
                                           const std::vector< GroupInfo > & group_info,
                                           const std::vector< GeomElemBase* > geom_elem_per_dimension,
                                           const unsigned mesh_dim,
                                           MyMatrix <int> & element_faces_array
                              );
                                      
   void find_boundary_faces_and_set_face_flags(const hid_t&  file_id, 
                                               const std::string mesh_menu, 
                                               const std::vector<GroupInfo> & group_info, const GeomElemBase* geom_elem_per_dimension,
                                               MyMatrix <int> & element_faces_array);

   void find_boundary_nodes_and_set_node_flags(const hid_t&  file_id, 
                                               const std::string mesh_menu, 
                                               const std::vector<GroupInfo> & group_info,
                                               MyMatrix <int> & element_faces_array);
   
  void check_all_boundaries_have_some_condition(const unsigned int count_found_face, const unsigned int fam_size) const;

   
  bool  see_if_faces_from_different_lists_are_the_same( const GeomElemBase* geom_elem_per_dimension, 
                                            const std::vector< unsigned > & face_nodes_from_vol_connectivity, 
                                            const std::vector< unsigned > & face_nodes_from_bdry_group);
  
// === Mesh, Groups, Boundary - Elements or Nodes ( Elements(2D/3D) or Nodes(1D) ) - END =================


// === Mesh, Groups, Boundary of Boundary (some part) - Nodes (cannot select such elements) - BEGIN =================
 public:
  
  std::vector< TYPE_FOR_REAL_DATASET >  node_based_flag_read_from_file(const std::string& name, const std::vector< unsigned > & mapping);
  
  void boundary_of_boundary_3d_via_nodes(const std::string& name, const unsigned group_user);
  
  static bool boundary_of_boundary_3d_check_face_of_face_via_nodes(const std::vector < int > nodes_face_face_flags, const unsigned group_salome);

 private:
     
  void all_nodes_read_group_flag(const hid_t&  file_id, const std::string mesh_menu,  std::vector < TYPE_FOR_REAL_DATASET >  & node_group_map);
// === Mesh, Groups, Boundary of Boundary (some part) - Nodes (cannot select such elements) - END =================


  
// === Mesh, Groups - END =================

   
};




inline
MED_IO::MED_IO (Mesh& mesh) :
   MeshInput<Mesh>  (mesh)
{
    
       _geom_elems.resize(N_GEOM_ELS);
       
       _geom_elems[0] = new GeomElemHex27();
       _geom_elems[1] = new GeomElemTet10();
//        _geom_elems[2] = new GeomElemWedge18();
       _geom_elems[3] = new GeomElemQuad9();
       _geom_elems[4] = new GeomElemTri6();
       _geom_elems[5] = new GeomElemEdge3();
    
    _print_info = false;
      
}

inline
MED_IO::~MED_IO () {
  
       delete _geom_elems[0];
       delete _geom_elems[1];
//     delete _geom_elems[2];
       delete _geom_elems[3];
       delete _geom_elems[4];
       delete _geom_elems[5];

    std::vector< GeomElemBase* > ().swap(_geom_elems);
    
}

   
  

// MED_IO::~MED_IO() {
//     
//         for(unsigned g = 0; g < _geom_elems.size(); g++) delete _geom_elems[g];
// 
// }

// @todo hybrid meshes (MED can export them)
// @todo groups with more than one type of element (in the sense hybrid, but at the same dimension)
// @todo mesh with overlapping groups  (MED can export them, defining the intersection independently; now the name parsing wouldn't work with them; Salome splits into non-intersecting family, so one needs to read the GRO/NOM dataset )
// @todo how would I use a 0d element? or a ball?
// @todo triggering of the copy of new input files happens right after "Configure" of cmake
// @todo Redirect CERR together with COUT
// @todo Update Xdmf to Xdmf3
// @todo not running in parallel
// @todo I know I can make only groups of edges/faces/volumes, but can I do hybrid faces, say quad and tri together?
// @todo Generate MED format from other software such as Gmsh
// @todo FEHex20, FEQuad8, FETri7 missing

} // namespace femus

#endif 

#endif

