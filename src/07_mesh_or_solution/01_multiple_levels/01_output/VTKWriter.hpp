/*=========================================================================

 Program: FEMUS
 Module: VTKWriter
 Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_solution_VTKWriter_hpp__
#define __femus_solution_VTKWriter_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer.hpp"
#include "Mesh.hpp"

#include <b64/b64.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <iomanip>
#include <algorithm>


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class Mesh;


class VTKWriter : public Writer {


// === Constructors / Destructor  - BEGIN =================
public:

    /** Constructor. */
    VTKWriter(const MultiLevelSolution * ml_sol);

    /** Constructor. */
    VTKWriter(const MultiLevelMesh * ml_mesh);

// === Constructors / Destructor  - END =================


// === Write, at finest level - BEGIN =================
public:

    /** write output function */
    void Write(const std::string output_path, 
               const std::string order,
               const std::vector < std::string > & vars = std::vector < std::string > (), 
               const unsigned time_step = Writer_one_level::_time_step_index_default) ;
    
    /** write output function with fixed level and arbitrary initial string */
    void Write(const std::string filename_prefix, 
               const std::string output_path,
               const std::string order,
               const std::vector < std::string >& vars = std::vector < std::string > (), 
               const unsigned time_step = Writer_one_level::_time_step_index_default);
private:

  void Write(const std::string filename_prefix, 
                        const std::string output_path,
                        const std::string suffix_pre_extension,
                        const std::string order,
                        const std::vector < std::string >& vars, 
                        const unsigned time_step);
// === Write, at finest level - END =================


// === Write, at arbitrary level - BEGIN =================
private:

  
      void Write(const unsigned level_in,
               const std::string filename_prefix,
               const std::string output_path,
               const std::string suffix_pre_extension,
               const std::string order,
               const std::vector < std::string >& vars = std::vector < std::string > (),
               const unsigned time_step = Writer_one_level::_time_step_index_default);



      /** write output function with arbitrary level (starting at 1) and arbitrary initial string and arbitrary suffix before the extension */
    void Write(const unsigned level_in, 
               const std::string filename_prefix, 
               const std::string output_path,
               const std::string suffix_pre_extension, 
               const std::string order,
                        const Mesh * mesh_in,
                        const Solution * solution_in,
               const std::vector < std::string >& vars, 
               const unsigned time_step);
// === Write, at arbitrary level - END =================
    

  private:
      
    void vtk_unstructured_header_parallel_wrapper(std::ofstream & Pfout) const;
    
    void vtk_unstructured_footer_parallel_wrapper(std::ofstream & Pfout) const;
    
    void vtk_unstructured_header_iproc(std::ofstream & fout) const;
      
    void vtk_unstructured_footer_iproc(std::ofstream & fout) const;
    
    void piece_iproc_begin(std::ofstream & fout, const unsigned n_nodes, const unsigned n_elements) const;
  
    void piece_iproc_end(std::ofstream & fout) const;
    
    void pieces_list_sources(std::ofstream & Pfout,
                                        const std::string dirnamePVTK,
                                        const std::string filename_prefix,
                                        const std::string level_name,
                                        const unsigned time_step,
                                        const std::string order,
                                        const std::string suffix_pre_extension ) const;
                                        

    std::map < unsigned, unsigned > ghost_map_proc(const Mesh * mesh, const unsigned index) const;
    
    void fill_connectivity_proc(const Mesh * mesh, const unsigned index, const std::map<unsigned, unsigned> & ghostMap,  int * const var_conn) const;
  
    unsigned size_connectivity_proc(const Mesh * mesh, const unsigned index) const;
   
    bool print_all_sols(const std::vector < std::string >& vars) const;
  
    unsigned compute_print_sol_size(const bool print_all, const std::vector < std::string >& vars, const Solution * solution) const;
     
    void fill_sol_on_elements(const Mesh * mesh, 
                             const unsigned elementOffset, const unsigned elementOffsetp1, 
                             const Solution * solution, const unsigned name, const unsigned i,  float * const var_el) const;
                             
    template < class ARRAY_TYPE >     
    void print_Mesh_related_Element_based_fields(const std::string field_string,  const std::string field_datatype, std::ofstream & fout, std::ofstream & Pfout, void* buffer_void, const unsigned elemetOffset, const unsigned elemetOffsetp1, const unsigned * dim_array_elvar, const Mesh * mesh, const unsigned fe_index, std::vector <char> & enc ) const;
    
    template < class ARRAY_TYPE >     
    void print_data_array(const std::string field_string, 
                                   const std::string field_datatype,
                                   std::ofstream & fout, std::ofstream & Pfout,
                                   const unsigned * dim_array_elvar,
                                   const ARRAY_TYPE * var_el,
                                   std::vector <char> & enc) const;
                                       
    template < class ARRAY_TYPE >     
    void print_data_array_vector(const std::string field_string,
                               const std::string field_datatype,
                               const unsigned n_components,
                               std::ofstream & fout, std::ofstream & Pfout,
                               const unsigned * dim_array_elvar,
                               const ARRAY_TYPE * var_el,
                               std::vector <char> & enc) const;

    void print_coordinates(std::ofstream & fout, 
                                    std::ofstream & Pfout,
                                    NumericVector * num_vec_aux_for_node_fields,
                                    std::vector <char> & enc,
                                    void * buffer_void,
                                    const unsigned * dim_array_coord,
                                    const Mesh * mesh,
                                    const Solution * solution,
                                    const unsigned index,
                                    const unsigned nvtOwned,
                                    const std::map<unsigned, unsigned>  &  ghostMap);

    /** [Lagrange linear/quadratic/biquadratic][geom_elem_type]  */
    static short unsigned int femusToVtkCellType[3][6];

    
};





 template < class ARRAY_TYPE >     
  void VTKWriter::print_data_array(const std::string field_string, 
                                   const std::string field_datatype,
                                   std::ofstream & fout, std::ofstream & Pfout,
                                   const unsigned * dim_array_elvar,
                                   const ARRAY_TYPE * var_el,
                                   std::vector <char> & enc) const {  ///@todo do we really need to pass this guy?
                                                                   
                                                                   
    fout  << "       <DataArray type=\"" << field_datatype << "\" Name=\"" << field_string << "\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"" << field_datatype << "\" Name=\"" << field_string << "\" format=\"binary\"/>" << std::endl;
    
    //print solution on element dimension
    size_t cch = b64::b64_encode( &dim_array_elvar[0], sizeof( unsigned ) /*DO NOT USE THIS!!! sizeof( dim_array_elvar )*/, NULL, 0 );
    b64::b64_encode( &dim_array_elvar[0],  sizeof( unsigned )/*DO NOT USE THIS!!! sizeof( dim_array_elvar )*/, &enc[0], cch );
    char* pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) {  ///@todo do we have the guarantee that the std::vector is CONTIGUOUS in memory??
        fout << *pt_char;
    }
    
    //print solution on element array
    cch = b64::b64_encode( &var_el[0], dim_array_elvar[0] , NULL, 0 );
    b64::b64_encode( &var_el[0], dim_array_elvar[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) {
        fout << *pt_char;
    }
    
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    
   }
   
   
 template < class ARRAY_TYPE >     
  void VTKWriter::print_data_array_vector(const std::string field_string, 
                                   const std::string field_datatype,
                                   const unsigned n_components,
                                   std::ofstream & fout, std::ofstream & Pfout,
                                   const unsigned * dim_array_elvar,
                                   const ARRAY_TYPE * var_el,
                                   std::vector <char> & enc) const {

            
    fout  << "       <DataArray type=\"" << field_datatype << "\" Name=\"" << field_string  << "\" NumberOfComponents=\"" << n_components << "\" format=\"binary\">"  << std::endl;
    Pfout << "      <PDataArray type=\"" << field_datatype << "\" Name=\"" << field_string  << "\" NumberOfComponents=\"" << n_components << "\" format=\"binary\"/>" << std::endl;

    //print solution on element dimension
    size_t cch = b64::b64_encode( &dim_array_elvar[0], sizeof( unsigned ) /*DO NOT USE THIS!!! sizeof( dim_array_elvar )*/, NULL, 0 );
    b64::b64_encode( &dim_array_elvar[0],  sizeof( unsigned )/*DO NOT USE THIS!!! sizeof( dim_array_elvar )*/, &enc[0], cch );
    char* pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) {  ///@todo do we have the guarantee that the std::vector is CONTIGUOUS in memory??
        fout << *pt_char;
    }
    
    //print solution on element array
    cch = b64::b64_encode( &var_el[0], dim_array_elvar[0] , NULL, 0 );
    b64::b64_encode( &var_el[0], dim_array_elvar[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) {
        fout << *pt_char;
    }
    
    fout << std::endl;

    fout  << "        </DataArray>" << std::endl;
   
    
   }
   
   

 template < class ARRAY_TYPE >     
  void VTKWriter::print_Mesh_related_Element_based_fields(const std::string field_string, const std::string field_datatype,
                                                               std::ofstream & fout, std::ofstream & Pfout,
                                                               void* buffer_void,
                                                               const unsigned elemetOffset, const unsigned elemetOffsetp1,
                                                               const unsigned * dim_array_elvar,
                                                               const Mesh * mesh,
                                                               const unsigned fe_index,
                                                               std::vector <char> & enc) const {
                               

    // point pointer to common memory area buffer of void type;
    ARRAY_TYPE * var_el = static_cast< ARRAY_TYPE * >( buffer_void );
    
    int icount = 0;

    if      (field_string == "Material")           { for( unsigned iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {  var_el[icount] = mesh->GetElementMaterial(iel); icount++; }   }
    else if (field_string == "Group")              { for( unsigned iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {  var_el[icount] = mesh->GetElementGroup(iel); icount++; }      }
    else if (field_string == "TYPE")               { for( unsigned iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {  var_el[icount] = mesh->GetElementType(iel); icount++; }       }
    else if (field_string == "Level")              { for( unsigned iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {  var_el[icount] = mesh->GetMeshElements()->GetElementLevel(iel); icount++; }  }
    else if (field_string == "Metis partition")    { for( unsigned iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {  var_el[icount] = _writer_one_level.processor_id(); icount++; }  }
    
    else if (field_string == "types")              { for( unsigned iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {  short unsigned ielt = mesh->GetElementType( iel );
                                                                                                                        var_el[icount] = femusToVtkCellType[fe_index][ielt];
                                                                                                                        icount++; 
                                                                                                                     }  
                                                   }
    else if (field_string == "offsets")            { int offset_el = 0;
                                                     for( unsigned iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {   offset_el += mesh->GetElementDofNumber( iel, fe_index );  var_el[icount] = offset_el;  icount++; } }
    else    { abort(); }


      print_data_array< ARRAY_TYPE >(field_string, field_datatype, fout, Pfout, dim_array_elvar, var_el, enc);
  
  }


/**
 * VTK Cell type: from http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
 *
 * VTK_EMPTY_CELL                         = 0,
 * VTK_VERTEX                             = 1,
 * VTK_POLY_VERTEX                        = 2,
 * VTK_LINE                               = 3,
 * VTK_POLY_LINE                          = 4,
 * VTK_TRIANGLE                           = 5,
 * VTK_TRIANGLE_STRIP                     = 6,
 * VTK_POLYGON                            = 7,
 * VTK_PIXEL                              = 8,
 * VTK_QUAD                               = 9,
 * VTK_TETRA                              = 10,
 * VTK_VOXEL                              = 11,
 * VTK_HEXAHEDRON                         = 12,
 * VTK_WEDGE                              = 13,
 * VTK_PYRAMID                            = 14,
 * VTK_PENTAGONAL_PRISM                   = 15,
 * VTK_HEXAGONAL_PRISM                    = 16,
 * VTK_QUADRATIC_EDGE                     = 21,
 * VTK_QUADRATIC_TRIANGLE                 = 22,
 * VTK_QUADRATIC_QUAD                     = 23,  // 23:Serendipity(8-nodes)
 * VTK_QUADRATIC_POLYGON                  = 36,
 * VTK_QUADRATIC_TETRA                    = 24,
 * VTK_QUADRATIC_HEXAHEDRON               = 25,
 * VTK_QUADRATIC_WEDGE                    = 26,
 * VTK_QUADRATIC_PYRAMID                  = 27,
 * VTK_BIQUADRATIC_QUAD                   = 28,  // 28:Quad9-Biquadratic
 * VTK_TRIQUADRATIC_HEXAHEDRON            = 29,
 * VTK_QUADRATIC_LINEAR_QUAD              = 30,
 * VTK_QUADRATIC_LINEAR_WEDGE             = 31,
 * VTK_BIQUADRATIC_QUADRATIC_WEDGE        = 32,
 * VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON   = 33,
 * VTK_BIQUADRATIC_TRIANGLE               = 34,
 * VTK_CUBIC_LINE                         = 35,
 * VTK_CONVEX_POINT_SET                   = 41,
 * VTK_POLYHEDRON                         = 42,
 * VTK_PARAMETRIC_CURVE                   = 51,
 * VTK_PARAMETRIC_SURFACE                 = 52,
 * VTK_PARAMETRIC_TRI_SURFACE             = 53,
 * VTK_PARAMETRIC_QUAD_SURFACE            = 54,
 * VTK_PARAMETRIC_TETRA_REGION            = 55,
 * VTK_PARAMETRIC_HEX_REGION              = 56,
 * VTK_HIGHER_ORDER_EDGE                  = 60,
 * VTK_HIGHER_ORDER_TRIANGLE              = 61,
 * VTK_HIGHER_ORDER_QUAD                  = 62,
 * VTK_HIGHER_ORDER_POLYGON               = 63,
 * VTK_HIGHER_ORDER_TETRAHEDRON           = 64,
 * VTK_HIGHER_ORDER_WEDGE                 = 65,
 * VTK_HIGHER_ORDER_PYRAMID               = 66,
 * VTK_HIGHER_ORDER_HEXAHEDRON            = 67,
 */


} //end namespace femus



#endif
