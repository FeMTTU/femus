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
   
   /** Map from Gambit vertex index to Femus vertex index */
   static const unsigned GambitToFemusVertexIndex[6][27]; 
 
   /** Map from Gambit face index to Femus face index */
   static const unsigned GambitToFemusFaceIndex[6][6];

   /** Read FE type */
   int  ReadFE(hid_t file_id, std::string & el_fem_type_vol, std::string & el_fem_type_bd, const  std::string menu_name);   
   
   std::vector<char*> menu_names;
   static const std::string group_name_begin; //FAS
   static const std::string group_name_end;   //ELEME
   static const std::string mesh_ensemble;        // ENS_MAA
   static const std::string aux_zeroone;      // -0000000000000000001-0000000000000000001
   static const std::string connectivity;     // MAI
   static const std::string node_coord;       //NOE/COO
   static const uint max_length;

};


inline
SalomeIO::SalomeIO (Mesh& mesh) :
   MeshInput<Mesh>  (mesh)
{
}


} // namespace femus

#endif 