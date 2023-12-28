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

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "VTKWriter.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "Files.hpp"



namespace femus {


  short unsigned int VTKWriter::femusToVtkCellType[3][6] = {{12, 10, 13, 9, 5, 3}, {25, 24, 26, 23, 22, 21}, {29, 24, 32, 28, 34, 21}};
  //http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html#ab1d6fd1f3177b8a2a32bb018807151f8aff535f3b1a33b5e51d1ef1e3aed69447

  VTKWriter::VTKWriter(const MultiLevelSolution* ml_sol ): Writer( ml_sol ) {  }

  VTKWriter::VTKWriter(const MultiLevelMesh* ml_mesh ): Writer( ml_mesh ) {  }

  
    void VTKWriter::vtk_unstructured_header_parallel_wrapper(std::ofstream & Pfout) const {
        
    Pfout << "<?xml version=\"1.0\"?>" << std::endl;
    Pfout << "<VTKFile type = \"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    Pfout << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
    
    }

    void VTKWriter::vtk_unstructured_footer_parallel_wrapper(std::ofstream & Pfout) const {
    
    Pfout << "  </PUnstructuredGrid>" << std::endl;
    Pfout << "</VTKFile>" << std::endl;
        
    }
    

    void VTKWriter::vtk_unstructured_header_iproc(std::ofstream & fout) const {
        
    // *********** write vtu header ************
    fout << "<?xml version=\"1.0\"?>" << std::endl;
    fout << "<VTKFile type = \"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    fout << "  <UnstructuredGrid>" << std::endl;
    
    }

    void VTKWriter::vtk_unstructured_footer_iproc(std::ofstream & fout) const {
    
    fout << "  </UnstructuredGrid>" << std::endl;
    fout << "</VTKFile>" << std::endl;
        
    }
    
    
    void VTKWriter::pieces_list_sources(std::ofstream & Pfout,
                                        const std::string dirnamePVTK,
                                        const std::string filename_prefix,
                                        const std::string level_name,
                                        const unsigned time_step,
                                        const std::string order,
                                        const std::string suffix_pre_extension
                                        ) const {
      
    for( int jproc = 0; jproc < _writer_one_level.n_processors(); jproc++ ) {
      Pfout << "    <Piece Source=\"" << dirnamePVTK
            << filename_prefix << level_name << "." << jproc << "." << time_step << "." << order <<  suffix_pre_extension << ".vtu"
            << "\"/>" << std::endl;
      }
    
    }
    
    
    void VTKWriter::piece_iproc_begin(std::ofstream & fout, const unsigned n_nodes, const unsigned n_elements) const {
         
        fout  << "    <Piece NumberOfPoints= \"" << n_nodes << "\" NumberOfCells= \"" << n_elements << "\" >" << std::endl;
  
     }
     
     
   void VTKWriter::piece_iproc_end(std::ofstream & fout) const {
         
    fout << "    </Piece>" << std::endl;

   }
     
 
    
    
    std::map < unsigned, unsigned > VTKWriter::ghost_map_proc(const Mesh * mesh, const unsigned index) const {

    const unsigned elementOffset   = mesh->GetElementOffset( _writer_one_level.processor_id() );
    const unsigned elementOffsetp1 = mesh->GetElementOffset( _writer_one_level.processor_id() + 1 );
    
    const unsigned dofOffset = mesh->dofmap_get_dof_offset(index, _writer_one_level.processor_id());
    
    unsigned ghostMapCounter = 0;
     
    std::map < unsigned, unsigned > ghostMap;
    
    for( int iel = elementOffset; iel < elementOffsetp1; iel++ ) {

        for( unsigned j = 0; j < mesh->GetElementDofNumber( iel, index ); j++ ) {

            unsigned jdof = mesh->GetSolutionDof( j, iel, index );
        if( jdof < dofOffset ) { // check if jdof is a ghost node
          if( ghostMap.find( jdof ) == ghostMap.end() ) {
            ghostMap[jdof] = ghostMapCounter;
            ghostMapCounter++;
          }
        }
      }
    }
 
      return ghostMap;
 
    }
 

  void VTKWriter::fill_connectivity_proc(const Mesh * mesh, const unsigned index, const std::map<unsigned, unsigned>  &  ghostMap,  int * const var_conn) const {
       
    const unsigned elementOffset = mesh->GetElementOffset( _writer_one_level.processor_id() );
    const unsigned elementOffsetp1 = mesh->GetElementOffset( _writer_one_level.processor_id() + 1 );
    
    const unsigned dofOffset = mesh->dofmap_get_dof_offset(index, _writer_one_level.processor_id());
    const unsigned nvtOwned = mesh->dofmap_get_own_size(index, _writer_one_level.processor_id());
    
    // point pointer to common memory area buffer of void type;
    unsigned icount = 0;
    
    for(unsigned iel = elementOffset; iel < elementOffsetp1; iel++ ) {
      for( unsigned j = 0; j < mesh->GetElementDofNumber( iel, index ); j++ ) {
        unsigned loc_vtk_conn = (mesh->GetElementType( iel ) == 0)? Writer_one_level::FemusToVTKorToXDMFConn[j] : j;
        unsigned jdof = mesh->GetSolutionDof( loc_vtk_conn, iel, index );
        var_conn[icount] = ( jdof >= dofOffset ) ? jdof - dofOffset : nvtOwned + ghostMap.find(jdof)->second;
        icount++;
        }
      }
    
  }
  
  
 
   void VTKWriter::fill_sol_on_elements(const Mesh * mesh, 
                                        const unsigned elementOffset, const unsigned elementOffsetp1, 
                                        const Solution * solution, const unsigned name, const unsigned i,  float * const var_el) const {
       
            unsigned icount = 0;
            
            for( int iel = elementOffset; iel < elementOffsetp1; iel++ ) {
                
              unsigned placeholder_index = 0;
              unsigned iel_Metis = mesh->GetSolutionDof(placeholder_index, iel, solution->GetSolutionType( i ) );
              
              if( name == Writer_one_level::_index_sol )                 var_el[icount] = ( *solution->_Sol[i] )( iel_Metis );
              else if( name == Writer_one_level::_index_bdc )            var_el[icount] = ( *solution->_Bdc[i] )( iel_Metis );
              else if( name == Writer_one_level::_index_res )            var_el[icount] = ( *solution->_Res[i] )( iel_Metis );
              else if( name == Writer_one_level::_index_eps )            var_el[icount] = ( *solution->_Eps[i] )( iel_Metis );
              
              icount++;
              
            }
            
            
   }
   
   
 
   unsigned VTKWriter::size_connectivity_proc(const Mesh * mesh, const unsigned index) const {
       
      unsigned counter = 0;
      
      for(unsigned iel = mesh->GetElementOffset(_writer_one_level.processor_id()); iel < mesh->GetElementOffset(_writer_one_level.processor_id() + 1); iel++ ) {

        for( unsigned j = 0; j < mesh->GetElementDofNumber( iel, index ); j++ ) {
            counter++;
        }
      }
      
        return counter;
        
   }
   
   
   

   
   

   
   
   bool VTKWriter::print_all_sols(const std::vector < std::string >& vars) const {
       
    bool print_all = 0;
    for( unsigned ivar = 0; ivar < vars.size(); ivar++ ) {
      print_all += !( vars[ivar].compare( "All" ) ) + !( vars[ivar].compare( "all" ) ) + !( vars[ivar].compare( "ALL" ) );
    }
    
    return print_all;
    
   }
   
   
   unsigned VTKWriter::compute_print_sol_size(const bool print_all, const std::vector < std::string >& vars, const Solution * solution) const {
       
      return ( !print_all ) * vars.size() + print_all * solution->GetSolutionSize();

   }

   

  void VTKWriter::print_coordinates(std::ofstream & fout, 
                                    std::ofstream & Pfout,
                                    NumericVector * num_vec_aux_for_node_fields,
                                    std::vector <char> & enc,
                                    void * buffer_void,
                                    const unsigned * dim_array_coord,
                                    const Mesh * mesh,
                                    const Solution * solution,
                                    const unsigned index,
                                    const unsigned nvtOwned,
                                    const std::map<unsigned, unsigned>  &  ghostMap)
                                    {
                                      

    fout  << "     <Points>" << std::endl;
    Pfout << "    <PPoints>" << std::endl;
    
    
    // point pointer to common memory area buffer of void type;
    float * var_coord = static_cast<float*>( buffer_void );
    
    //var_coord: add own nodes - BEGIN -------------------------

    const unsigned dofOffset = mesh->dofmap_get_dof_offset(index, _writer_one_level.processor_id());
    
    for( int i = 0; i < 3; i++ ) {
      
    // num_vec_aux_for_node_fields - BEGIN -------------------------
      if( ! _writer_one_level.is_surface() ) {
        num_vec_aux_for_node_fields->matrix_mult( *mesh->GetTopology()->_Sol[i],    * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, 2, * mesh )  );
        if( solution != NULL && _writer_one_level._graph && i == 2 ) {
          const unsigned indGraph = solution->GetIndex( _writer_one_level._graphVariable.c_str() );
          num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indGraph],  * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indGraph ), * mesh ) );
        }
      }
      
      else if (_writer_one_level.is_surface() && solution != NULL ) {
        const unsigned indSurfVar = solution->GetIndex( _writer_one_level._surfaceVariables[i].c_str() );
        num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indSurfVar], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indSurfVar ), * mesh ) );
      }
    // num_vec_aux_for_node_fields - END -------------------------
      
    // var_coord - BEGIN -------------------------
      for( unsigned ii = 0; ii < nvtOwned; ii++ ) {
        var_coord[ ii * 3 + i] = ( *num_vec_aux_for_node_fields )( ii +  dofOffset );
      }
    // var_coord - END -------------------------
      
      if( solution != NULL && _writer_one_level._moving_mesh  && _writer_one_level._moving_vars.size() > i) { // if moving mesh

    // num_vec_aux_for_node_fields - BEGIN -------------------------
        const unsigned indDXDYDZ = solution->GetIndex( _writer_one_level._moving_vars[i].c_str() );
        
        num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indDXDYDZ],  * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indDXDYDZ ), * mesh ) );
    // num_vec_aux_for_node_fields - END -------------------------

    // var_coord - BEGIN -------------------------
        for( unsigned ii = 0; ii < nvtOwned; ii++ ) {
          var_coord[ii * 3 + i] += ( *num_vec_aux_for_node_fields )( ii +  dofOffset );
        }
    // var_coord - END -------------------------
        
      }
      
      
    }
    //var_coord: add own nodes - END -------------------------
    
    
    //var_coord: add ghost nodes - BEGIN -------------------------
    const unsigned offset_ig = 3 * nvtOwned;

    for( int i = 0; i < 3; i++ ) {
      
    // num_vec_aux_for_node_fields - BEGIN -------------------------
      if( ! _writer_one_level.is_surface() ) {
          num_vec_aux_for_node_fields->matrix_mult( *mesh->GetTopology()->_Sol[i],  * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, 2, * mesh ) );
        if( solution != NULL &&  _writer_one_level._graph && i == 2 ) {
          const unsigned indGraphVar = solution->GetIndex( _writer_one_level._graphVariable.c_str() );
          num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indGraphVar], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indGraphVar ), * mesh ) );
        }
      }
      else if ( _writer_one_level.is_surface() && solution != NULL ) {
        const unsigned indSurfVar = solution->GetIndex( _writer_one_level._surfaceVariables[i].c_str() );
        num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indSurfVar],  * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indSurfVar ), * mesh ) );
      }
    // num_vec_aux_for_node_fields - END -------------------------
      
    // var_coord - BEGIN -------------------------
      for( std::map <unsigned, unsigned>::const_iterator it = ghostMap.begin(); it != ghostMap.end(); ++it ) {
        var_coord[ offset_ig + 3 * it->second + i ] = ( *num_vec_aux_for_node_fields )( it->first );
      }
    // var_coord - END -------------------------
      
    }

    for( int i = 0; i < 3; i++ ) { // if moving mesh

      if( solution != NULL && _writer_one_level._moving_mesh  && _writer_one_level._moving_vars.size() > i ) { //&& mesh->GetDimension() > i )  {
        
    // num_vec_aux_for_node_fields - BEGIN -------------------------
        const unsigned indDXDYDZ = solution->GetIndex( _writer_one_level._moving_vars[i].c_str() );

        num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indDXDYDZ], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indDXDYDZ ), * mesh ) );
    // num_vec_aux_for_node_fields - END -------------------------
        
    // var_coord - BEGIN -------------------------
        for( std::map <unsigned, unsigned>::const_iterator it = ghostMap.begin(); it != ghostMap.end(); ++it ) {
          var_coord[ offset_ig + 3 * it->second + i ] += ( *num_vec_aux_for_node_fields )( it->first );
        }
    // var_coord - END -------------------------
        
      }
      
    }
    //var_coord: add ghost nodes - END -------------------------
    

    print_data_array_vector< float >("coordinates", "Float32", 3, fout, Pfout, dim_array_coord, var_coord, enc);
    
    
    fout  << "      </Points>" << std::endl;
    Pfout << "    </PPoints>" << std::endl;
    // print coordinates - END ****************************************************************************************
                                      
                                      
  }





 void VTKWriter::Write(const std::string output_path, 
                         const std::string order,
                         const std::vector < std::string >& vars,
                         const unsigned time_step ) {
       
    const Solution * solution = get_solution(_gridn);

    const std::string filename_prefix = _writer_one_level.get_filename_prefix(solution);
    
    const std::string suffix_pre_extension = "";
    
    Write(_gridn, filename_prefix, output_path, suffix_pre_extension, order, vars, time_step );

   }
   
   
   
   
  void VTKWriter::Write(const std::string filename_prefix, 
                        const std::string output_path,
                        const std::string order,
                        const std::vector < std::string >& vars,
                        const unsigned time_step ) {
       
      const std::string suffix_pre_extension = "";
      
       Write(_gridn, filename_prefix, output_path, suffix_pre_extension, order, vars, time_step );
   }


  void VTKWriter::Write(const std::string filename_prefix, 
                        const std::string output_path,
                        const std::string suffix_pre_extension,
                        const std::string order,
                        const std::vector < std::string >& vars, 
                        const unsigned time_step) {
    
        Write(_gridn, filename_prefix, output_path, suffix_pre_extension, order, vars, time_step );
 }


  

  void VTKWriter::Write(const unsigned my_level, 
                        const std::string filename_prefix, 
                        const std::string output_path,
                        const std::string suffix_pre_extension,
                        const std::string order,
                        const std::vector < std::string >& vars, 
                        const unsigned time_step) {
    
    //------------- Mesh, def - BEGIN ----------------------------------------------------------------------------------
    const Mesh * mesh = get_mesh(my_level);
    //------------- Mesh, def - END ----------------------------------------------------------------------------------

    
    //------------- Solution, def - BEGIN ----------------------------------------------------------------------------------
    const Solution * solution = get_solution(my_level);
    //------------- Solution, def - END ----------------------------------------------------------------------------------

    Write(my_level, filename_prefix, output_path, suffix_pre_extension, order, mesh, solution, vars, time_step);
    
  }




  void VTKWriter::Write(const unsigned my_level, 
                        const std::string filename_prefix, 
                        const std::string output_path,
                        const std::string suffix_pre_extension,
                        const std::string order,
                        const Mesh * mesh_in,
                        const Solution * solution_in,
                        const std::vector < std::string >& vars, 
                        const unsigned time_step
) {

    //------------- Mesh, def - BEGIN ----------------------------------------------------------------------------------
    const Mesh * mesh = mesh_in;
    //------------- Mesh, def - END ----------------------------------------------------------------------------------

    
    //------------- Solution, def - BEGIN ----------------------------------------------------------------------------------
    const Solution * solution = solution_in;
    //------------- Solution, def - END ----------------------------------------------------------------------------------


    // *********** level - BEGIN ************      
    std::ostringstream level_name_stream;    
    level_name_stream << ".level" << my_level;
    std::string level_name(level_name_stream.str());   
    // *********** level - END ************
       

    
    // *********** FE index - BEGIN ************
    const unsigned index = Writer_one_level::fe_index(order);
    // *********** FE index - END ************


    // *********** open vtu streams - BEGIN *************
    const std::string dirnamePVTK = "VTKParallelFiles/";
    Files files;
    files.CheckDir( output_path, "" );
    files.CheckDir( output_path, dirnamePVTK );

    std::ofstream fout;

    std::ostringstream filename;
    filename << output_path << "./" << dirnamePVTK << filename_prefix << level_name << "." << _writer_one_level.processor_id() << "." << time_step << "." << order << suffix_pre_extension << ".vtu";

    fout.open( filename.str().c_str() );
    if( !fout.is_open() ) {
      std::cout << std::endl << " The output file " << filename.str() << " cannot be opened.\n";
      abort();
    }
    // *********** open vtu streams - END *************


    // *********** open pvtu stream - BEGIN *************
    std::ofstream Pfout;
    if( _writer_one_level.processor_id() != 0 ) {
      Pfout.rdbuf();   //redirect to dev_null
    }
    else {
      std::ostringstream Pfilename;
      Pfilename << output_path << "./" << filename_prefix << level_name << "." << time_step << "." << order <<  suffix_pre_extension << ".pvtu";
      Pfout.open( Pfilename.str().c_str() );
      if( Pfout.is_open() ) {
        std::cout << std::endl << " The output is printed to file " << Pfilename.str() << " in parallel VTK-XML (64-based) format" << std::endl;
      }
      else {
        std::cout << std::endl << " The output file " << Pfilename.str() << " cannot be opened.\n";
        abort();
      }
    }
    // *********** open pvtu stream - END *************
    

    // *********** write vtu header - BEGIN ************
    vtk_unstructured_header_iproc(fout);
    // *********** write vtu header - END ************
    
    // *********** write pvtu header - BEGIN ***********
    vtk_unstructured_header_parallel_wrapper(Pfout);
    // *********** write pvtu header - END ***********
    
    
    //----------- ALL IPROCS - BEGIN ------------------------------------------------------------------------------------
   pieces_list_sources(Pfout, dirnamePVTK, filename_prefix, level_name, time_step, order, suffix_pre_extension);
    //----------- ALL IPROCS - END ------------------------------------------------------------------------------------

    

    //------------- Mesh, NODE and ELEMENT INFO - BEGIN ----------------------------------------------------------------------------------
    // count the own element dofs on all levels -------------
    const unsigned elemetOffset = mesh->GetElementOffset(_writer_one_level.processor_id() );
    const unsigned elemetOffsetp1 = mesh->GetElementOffset(_writer_one_level.processor_id() + 1 );
    const unsigned nel = elemetOffsetp1 - elemetOffset;
    
    //count the own node dofs on all levels -------------
    const unsigned nvtOwned = mesh->dofmap_get_own_size(index, _writer_one_level.processor_id());
    
    // count the ghost node dofs on all levels -------------
    const std::map < unsigned, unsigned > ghostMap = ghost_map_proc(mesh, index);

    // count the total node dofs (own + ghost) on all levels -------------
    const unsigned nvt = nvtOwned + ghostMap.size();
    

    const unsigned size_conn = size_connectivity_proc(mesh, index);
    
    const unsigned dim_array_coord [] = { nvt * 3 * static_cast<unsigned>(sizeof( float )) };
    const unsigned dim_array_conn[]   = { size_conn * static_cast<unsigned>(sizeof( int )) };
    const unsigned dim_array_off []   = { nel * static_cast<unsigned>(sizeof( int )) };
    const unsigned dim_array_type []  = { nel * static_cast<unsigned>(sizeof( short unsigned )) };
    const unsigned dim_array_reg []   = { nel * static_cast<unsigned>(sizeof( short unsigned )) };
    const unsigned dim_array_elvar [] = { nel * static_cast<unsigned>(sizeof( float )) };
    const unsigned dim_array_ndvar [] = { nvt * static_cast<unsigned>(sizeof( float )) };

    // initialize common buffer_void memory
    const unsigned buffer_size = ( dim_array_coord[0] > dim_array_conn[0] ) ? dim_array_coord[0] : dim_array_conn[0];  // maximum size
    void* buffer_void = new char [buffer_size];   /// @todo is this deleted? I don't think so
    char* buffer_char = static_cast <char*>( buffer_void );

    // size of the base64-encoded data
    const size_t cch = b64::b64_encode( &buffer_char[0], buffer_size , NULL, 0 );
    std::vector <char> enc;
    enc.resize( cch );
    //------------- Mesh, NODE and ELEMENT INFO - END ----------------------------------------------------------------------------------
    
    
    
    //---- NumericVector used for node-based fields - BEGIN -------------------------------------------------------------------------------------------
    NumericVector* num_vec_aux_for_node_fields;
    num_vec_aux_for_node_fields = NumericVector::build().release();

    if( _writer_one_level.n_processors() == 1 ) { // IF SERIAL
      num_vec_aux_for_node_fields->init( mesh->dofmap_get_dof_offset(index, _writer_one_level.n_processors()),
                                         mesh->dofmap_get_dof_offset(index, _writer_one_level.n_processors()), false, SERIAL );
    }
    else { // IF PARALLEL
      num_vec_aux_for_node_fields->init( mesh->dofmap_get_dof_offset(index, _writer_one_level.n_processors()),
                                         mesh->dofmap_get_own_size(index, _writer_one_level.processor_id()),
                                         mesh->dofmap_get_ghost_dofs(index, _writer_one_level.processor_id()), false, GHOSTED );
    }
    //---- NumericVector used for node-based fields - END -------------------------------------------------------------------------------------------

    
    //----------- IPROC - BEGIN ------------------------------------------------------------------------------------
    piece_iproc_begin(fout, nvt, nel);

    
    // print coordinates - BEGIN ****************************************************************************************
    print_coordinates(fout, Pfout,
                      num_vec_aux_for_node_fields,
                      enc,
                      buffer_void,
                      dim_array_coord,
                      mesh,
                      solution,
                      index,
                      nvtOwned,
                      ghostMap);
    // print coordinates - END ****************************************************************************************
    
    
    //----- Printing of element connectivity - offset - format type  * - BEGIN ----------------------------------------

    fout  << "      <Cells>" << std::endl;
    Pfout << "    <PCells>" << std::endl;
    
    //-----------------------------------------------------------------------------------------------
    //print connectivity
    int * var_conn = static_cast <int*>( buffer_void );
    
     fill_connectivity_proc(mesh, index, ghostMap, var_conn);

     print_data_array< int >("connectivity", "Int32", fout, Pfout, dim_array_conn, var_conn, enc);


    //-------------------------------------------------------------------------------------------------
    //printing offset
    print_Mesh_related_Element_based_fields< int >("offsets", "Int32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_off, mesh, index, enc);

    //--------------------------------------------------------------------------------------------------
    //Element format type
    print_Mesh_related_Element_based_fields< unsigned short >("types", "UInt16", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_type, mesh, index, enc);
    

    fout  << "      </Cells>" << std::endl;
    Pfout << "    </PCells>" << std::endl;
    //----- Printing of element connectivity - offset - format type  * - END ------------------------------------------


    // /Print Cell Data - BEGIN ****************************************************************************
    fout  << "      <CellData Scalars=\"scalars\">" << std::endl;
    Pfout << "    <PCellData Scalars=\"scalars\">" << std::endl;
    

    //---------------------------- PARALLEL PARTITION, MATERIAL, GROUP, FE TYPE, LEVEL - BEGIN -----------------------
    
    print_Mesh_related_Element_based_fields< unsigned short >("Metis partition", "UInt16", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_reg, mesh, index, enc);

    print_Mesh_related_Element_based_fields< float >("Material", "Float32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, index, enc);

    print_Mesh_related_Element_based_fields< float >("Group", "Float32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, index, enc);

    print_Mesh_related_Element_based_fields< float >("TYPE", "Float32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, index, enc);

    print_Mesh_related_Element_based_fields< float >("Level", "Float32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, index, enc);

    //---------------------------- PARALLEL PARTITION, MATERIAL, GROUP, FE TYPE, LEVEL - END -----------------------
    


   const bool print_all = print_all_sols(vars);
    
    
    //------------------------------------------- SOLUTIONS on ELEMENTS - BEGIN ---------------------------------------------------------
    
    if( solution != NULL ) {
        
      //Print Solution (on elements) - BEGIN ***************************************************************
      for( unsigned i = 0; i < compute_print_sol_size(print_all, vars, solution); i++ ) {
          
        const unsigned solIndex = ( print_all == 0 ) ? solution->GetIndex( vars[i].c_str() ) : i;
        const unsigned sol_fe_type = solution->GetSolutionType( solIndex );
        
        if( 3 <= sol_fe_type ) {

          std::string solName =  solution->GetSolName_from_index( solIndex );

          for( int name = 0; name < _writer_one_level.compute_sol_bdc_res_eps_size(solution, i); name++ ) {
              
            std::string printName = Writer_one_level::print_sol_bdc_res_eps_name(solName, name);
            

            //--------- fill var_el ------------
            float* var_el = static_cast< float*>( buffer_void );
            
            fill_sol_on_elements(mesh, elemetOffset, elemetOffsetp1, solution, name, i,  var_el);
            //--------- fill var_el ------------
            

    print_data_array< float >(printName, "Float32", fout, Pfout, dim_array_elvar, var_el, enc);

              
          }
        }
      }
      //Print Solution (on elements) - END ***************************************************************
      
    } //end solution != NULL
    
    //------------------------------------------- SOLUTIONS on ELEMENTS - END ---------------------------------------------------------

    fout  << "      </CellData>" << std::endl;
    Pfout << "    </PCellData>" << std::endl;

    // /Print Cell Data - END ****************************************************************************


    //------------------------------------------- SOLUTIONS on NODES - BEGIN ---------------------------------------------------------

    if( solution != NULL ) {
        
      fout  << "      <PointData Scalars=\"scalars\"> " << std::endl;
      Pfout << "    <PPointData Scalars=\"scalars\"> " << std::endl;
      
      
      // / Print Solution (on nodes) - BEGIN ********************************************************************
      //Loop on variables

      // point pointer to common memory area buffer of void type;
      float* var_nd = static_cast<float*>( buffer_void );
      
      for( unsigned i = 0; i < compute_print_sol_size(print_all, vars, solution); i++ ) {
          
        const unsigned solIndex = ( print_all == 0 ) ? solution->GetIndex( vars[i].c_str() ) : i;
        const unsigned sol_fe_type = solution->GetSolutionType( solIndex );
        
        if( sol_fe_type < NFE_FAMS_C_ZERO_LAGRANGE ) {
            
          //BEGIN LAGRANGIAN Fem SOLUTION
          std::string solName =  solution->GetSolName_from_index( solIndex );

          for( int name = 0; name < _writer_one_level.compute_sol_bdc_res_eps_size(solution, i); name++ ) {
              
            const std::string printName = Writer_one_level::print_sol_bdc_res_eps_name(solName, name);
           

         //--------- fill var_nd ------------
            //print own dofs -------------------------
            unsigned offset_iprc = mesh->dofmap_get_dof_offset(index, _writer_one_level.processor_id());
            unsigned nvt_ig = mesh->dofmap_get_own_size(index, _writer_one_level.processor_id());

            if( name == Writer_one_level::_index_sol )
              num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[solIndex], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( solIndex ), * mesh ) );
            else if( name == Writer_one_level::_index_bdc )
              num_vec_aux_for_node_fields->matrix_mult( *solution->_Bdc[solIndex], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( solIndex ), * mesh ) );
            else if( name == Writer_one_level::_index_res )
              num_vec_aux_for_node_fields->matrix_mult( *solution->_Res[solIndex], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( solIndex ), * mesh ) );
            else if( name == Writer_one_level::_index_eps )
              num_vec_aux_for_node_fields->matrix_mult( *solution->_Eps[solIndex], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( solIndex ), * mesh ) );

            for( unsigned ii = 0; ii < nvt_ig; ii++ ) {
              var_nd[ ii ] = ( *num_vec_aux_for_node_fields )( ii + offset_iprc );
            }
            
            //print ghost dofs -------------------------
            unsigned offset_ig = nvtOwned;

            for( std::map <unsigned, unsigned>::const_iterator it = ghostMap.begin(); it != ghostMap.end(); ++it ) {
              var_nd[ offset_ig + it->second ] = ( *num_vec_aux_for_node_fields )( it->first );
            }
            //--------- fill var_nd ------------
            

            print_data_array< float >(printName, "Float32", fout, Pfout, dim_array_ndvar, var_nd, enc);
            
          }
        } //endif
      } // end for sol
      
      delete [] var_nd;
      // / Print Solution (on nodes) - END ********************************************************************

      
      fout  << "      </PointData>" << std::endl;
      Pfout << "    </PPointData>" << std::endl;
      
    }  //end solution != NULL

    //------------------------------------------- SOLUTIONS on NODES - END ---------------------------------------------------------



    piece_iproc_end(fout);
    //----------- IPROC - END ------------------------------------------------------------------------------------


    // *********** write vtu footer and close stream - BEGIN ************
    vtk_unstructured_footer_iproc(fout);
    fout.close();
    // *********** write vtu footer and close stream - END ************

    // *********** write pvtu footer and close stream - BEGIN ***********
    vtk_unstructured_footer_parallel_wrapper(Pfout);
    Pfout.close();
    // *********** write pvtu footer and close stream - END ***********


    //-----------------------------------------------------------------------------------------------------
    //free memory
    delete num_vec_aux_for_node_fields;

    //--------------------------------------------------------------------------------------------------------
    return;
  }
  
  
  

} //end namespace femus


