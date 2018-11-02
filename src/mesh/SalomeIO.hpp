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

#ifndef __femus_mesh_SalomeIO_hpp__
#define __femus_mesh_SalomeIO_hpp__


// Local includes
#include "MeshInput.hpp"
#include "GeomElTypeEnum.hpp"

#ifdef HAVE_HDF5
  #include "hdf5.h"
#endif

namespace femus
{

// Forward declarations
class Mesh;

/**
 * This class implements writing meshes in the Gmsh format.
 */

// ------------------------------------------------------------
// GMVIO class definition
class SalomeIO : public MeshInput<Mesh>
{
 public:

  /**
   * Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements via the read() command.
   */
  explicit
  SalomeIO (Mesh& mesh);

  /**
   * Reads in a mesh in the neutral gambit *.neu format
   * from the ASCII file given by name.
   *
   */
  virtual void read (const std::string& name, vector < vector < double> > &coords, const double Lref, std::vector<bool> &type_elem_flag);

 private:
   
   const std::vector< std::tuple<int,int,int,int> >  compute_group_flags_per_mesh(const std::vector<std::string> & group_names) const;
   
   const std::vector<std::string> compute_number_of_groups_per_mesh(const hid_t &  file_id, const std::string & mesh_menu) const;
   
   const std::vector<std::string>  compute_number_of_meshes(const hid_t & file_id) const;
     
   std::pair<int,int>  isolate_number_in_string(const std::string & string_in, const int begin_pos_to_investigate) const;
      
   /** Map from Salome vertex index to Femus vertex index */
   static const unsigned SalomeToFemusVertexIndex[N_GEOM_ELS][MAX_EL_N_NODES]; 
 
   /** Map from Salome face index to Femus face index */
   static const unsigned SalomeToFemusFaceIndex[N_GEOM_ELS][MAX_EL_N_FACES];

   /** Determine mesh dimension from mesh file */
   void  set_mesh_dimension_by_looping_over_element_types(const hid_t &  file_id, const std::vector<std::string> & menu_name, std::vector<std::string> & el_fe_type_per_dimension);    //this cannot be const because it sets the dimension in the mesh

   unsigned  FindNumberOfElemNodes(const  std::string el_type) const;

   /** Read FE type */
   void  ReadFE(const hid_t & file_id, std::vector<std::string> & fe_type_vec, const std::string my_mesh_name_dir) ;   //@todo this should be const
   
//    std::vector<char*> menu_names;
   static const std::string group_name_begin; //FAS
   static const std::string group_name_end;   //ELEME
   static const std::string mesh_ensemble;    // ENS_MAA
   static const std::string aux_zeroone;      // -0000000000000000001-0000000000000000001
   static const std::string elem_list;        //MAI
   static const std::string group_fam;        //FAM
   static const std::string connectivity;     //NOD
   static const std::string node_list;        //NOE
   static const std::string coord_list;       //COO
   static const std::string dofobj_indices;   //NUM
   static const uint max_length;

};


inline
SalomeIO::SalomeIO (Mesh& mesh) :
   MeshInput<Mesh>  (mesh)
{
}


} // namespace femus

#endif 
