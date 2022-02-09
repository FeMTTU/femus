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
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Files.hpp"


namespace femus {


  short unsigned int VTKWriter::femusToVtkCellType[3][6] = {{12, 10, 13, 9, 5, 3}, {25, 24, 26, 23, 22, 21}, {29, 24, 32, 28, 34, 21}};
  //http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html#ab1d6fd1f3177b8a2a32bb018807151f8aff535f3b1a33b5e51d1ef1e3aed69447

  VTKWriter::VTKWriter( MultiLevelSolution* ml_sol ): Writer( ml_sol ) {
    _debugOutput = false;
  }

  VTKWriter::VTKWriter( MultiLevelMesh* ml_mesh ): Writer( ml_mesh ) {
    _debugOutput = false;
  }

  VTKWriter::~VTKWriter(){}

  
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
    
    
    void VTKWriter::piece_iproc_begin(std::ofstream & fout, const unsigned n_nodes, const unsigned n_elements) const {
         
        fout  << "    <Piece NumberOfPoints= \"" << n_nodes << "\" NumberOfCells= \"" << n_elements << "\" >" << std::endl;
  
     }
     
     
   void VTKWriter::piece_iproc_end(std::ofstream & fout) const {
         
    fout << "    </Piece>" << std::endl;

   }
     
 
    unsigned VTKWriter::fe_index(const std::string & order_str) const {
        
        unsigned index = 0;
        
        if( !strcmp( order_str.c_str(), "linear" ) ) 	 index = 0; //linear
        else if( !strcmp( order_str.c_str(), "quadratic" ) ) 	 index = 1; //quadratic
        else if( !strcmp( order_str.c_str(), "biquadratic" ) ) index = 2; //biquadratic
        
        return index;
    }
    
    
    std::map < unsigned, unsigned > VTKWriter::ghost_map_proc(const Mesh * mesh, const unsigned index) const {

    unsigned elementOffset = mesh->_elementOffset[_iproc];
    unsigned elementOffsetp1 = mesh->_elementOffset[_iproc + 1];
    
    unsigned dofOffset = mesh->dofmap_get_dof_offset(index, _iproc);
    
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
       
    const unsigned elementOffset = mesh->_elementOffset[_iproc];
    const unsigned elementOffsetp1 = mesh->_elementOffset[_iproc + 1];
    
    const unsigned dofOffset = mesh->dofmap_get_dof_offset(index, _iproc);
    const unsigned nvtOwned = mesh->dofmap_get_own_size(index, _iproc);
    
    // point pointer to common memory area buffer of void type;
    unsigned icount = 0;
    
    for(unsigned iel = elementOffset; iel < elementOffsetp1; iel++ ) {
      for( unsigned j = 0; j < mesh->GetElementDofNumber( iel, index ); j++ ) {
        unsigned loc_vtk_conn = (mesh->GetElementType( iel ) == 0)? FemusToVTKorToXDMFConn[j] : j;
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
              unsigned iel_Metis = mesh->GetSolutionDof(placeholder_index, iel, _ml_sol->GetSolutionType( i ) );
              
              if( name == 0 )                 var_el[icount] = ( *solution->_Sol[i] )( iel_Metis );
              else if( name == 1 )            var_el[icount] = ( *solution->_Bdc[i] )( iel_Metis );
              else if( name == 2 )            var_el[icount] = ( *solution->_Res[i] )( iel_Metis );
              else                            var_el[icount] = ( *solution->_Eps[i] )( iel_Metis );
              icount++;
              
            }
            
            
   }
   
   
 
   unsigned VTKWriter::size_connectivity_proc(const Mesh * mesh, const unsigned index) const {
       
      unsigned counter = 0;
      
      for(unsigned iel = mesh->_elementOffset[_iproc]; iel < mesh->_elementOffset[_iproc + 1]; iel++ ) {

        for( unsigned j = 0; j < mesh->GetElementDofNumber( iel, index ); j++ ) {
            counter++;
        }
      }
      
        return counter;
        
   }
   
   
   
   std::string VTKWriter::print_sol_bdc_res_eps_name(const std::string solName, const unsigned name) const {
       
            std::string printName;

            if( name == 0 ) printName = solName;
            else if( name == 1 ) printName = "Bdc" + solName;
            else if( name == 2 ) printName = "Res" + solName;
            else printName = "Eps" + solName;
            
       return printName;     
   }
   
   
   unsigned VTKWriter::compute_sol_bdc_res_eps_size(const Solution * solution, const unsigned i) const {
       
       const unsigned print_sol_size = 1 + 3 * _debugOutput * solution->_ResEpsBdcFlag[i];
       return  print_sol_size;
       
   }
   
   
   bool VTKWriter::print_all_sols(const std::vector < std::string >& vars) const {
       
    bool print_all = 0;
    for( unsigned ivar = 0; ivar < vars.size(); ivar++ ) {
      print_all += !( vars[ivar].compare( "All" ) ) + !( vars[ivar].compare( "all" ) ) + !( vars[ivar].compare( "ALL" ) );
    }
    
    return print_all;
    
   }
   
   
   unsigned VTKWriter::compute_print_sol_size(const bool print_all, const std::vector < std::string >& vars) const {
       
      return ( !print_all ) * vars.size() + print_all * _ml_sol->GetSolutionSize();

   }
   
   
   void VTKWriter::Write(const std::string output_path, const char order[], const std::vector < std::string >& vars, const unsigned time_step ) {
       
    std::string filename_prefix;
    if( _ml_sol != NULL ) filename_prefix = "sol";
    else filename_prefix = "mesh";
    
      const std::string suffix_pre_extension = "";
    
       Write(_gridn, filename_prefix, output_path, suffix_pre_extension, order, vars, time_step );

   }
   

  void VTKWriter::Write(const unsigned my_level, const std::string output_path, const char order[], const std::vector < std::string >& vars, const unsigned time_step ) {
       
    std::string filename_prefix;
    if( _ml_sol != NULL ) filename_prefix = "sol";
    else filename_prefix = "mesh";
    
      const std::string suffix_pre_extension = "";
    
       Write(my_level, filename_prefix, output_path, suffix_pre_extension, order, vars, time_step );
       
   }
   
   
   
  void VTKWriter::Write(const std::string filename_prefix, const std::string output_path, const char order[], const std::vector < std::string >& vars, const unsigned time_step ) {
       
      const std::string suffix_pre_extension = "";
      
       Write(_gridn, filename_prefix, output_path, suffix_pre_extension, order, vars, time_step );
   }
   
  
  void VTKWriter::Write(const unsigned my_level, const std::string filename_prefix, const std::string output_path, const std::string suffix_pre_extension, const char order[], const std::vector < std::string >& vars, const unsigned time_step ) {

      
    // *********** level ************
    if (my_level < 1 || my_level > _gridn) { std::cout << "Level index in this routine is from 1 to num_levels" << std::endl; abort(); }  
      
    std::ostringstream level_name_stream;    
    level_name_stream << ".level" << my_level;
    std::string level_name(level_name_stream.str());   
       
    // *********** FE index ************
    const std::string order_str(order);
    unsigned index = fe_index(order_str);


    // *********** open vtu streams *************
    std::string dirnamePVTK = "VTKParallelFiles/";
    Files files;
    files.CheckDir( output_path, "" );
    files.CheckDir( output_path, dirnamePVTK );

    std::ofstream fout;

    std::ostringstream filename;
    filename << output_path << "./" << dirnamePVTK << filename_prefix << level_name << "." << _iproc << "." << time_step << "." << order << suffix_pre_extension << ".vtu";

    fout.open( filename.str().c_str() );
    if( !fout.is_open() ) {
      std::cout << std::endl << " The output file " << filename.str() << " cannot be opened.\n";
      abort();
    }


    // *********** open pvtu stream *************
    std::ofstream Pfout;
    if( _iproc != 0 ) {
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
    

    // *********** write vtu header ************
    vtk_unstructured_header_iproc(fout);
    
    // *********** write pvtu header ***********
    vtk_unstructured_header_parallel_wrapper(Pfout);
    
    for( int jproc = 0; jproc < _nprocs; jproc++ ) {
      Pfout << "    <Piece Source=\"" << dirnamePVTK
            << filename_prefix << level_name << "." << jproc << "." << time_step << "." << order <<  suffix_pre_extension << ".vtu"
            << "\"/>" << std::endl;
    }
    // ****************************************

    
    //------------- NODE and ELEMENT INFO ----------------------------------------------------------------------------------
    Mesh* mesh = _ml_mesh->GetLevel( my_level - 1 );
    Solution* solution;     if( _ml_sol != NULL ) { solution = _ml_sol->GetSolutionLevel( my_level - 1 ); }
    

    // count the own element dofs on all levels -------------
    const unsigned elemetOffset = mesh->_elementOffset[_iproc];
    const unsigned elemetOffsetp1 = mesh->_elementOffset[_iproc + 1];
    const unsigned nel = elemetOffsetp1 - elemetOffset;
    
    //count the own node dofs on all levels -------------
    unsigned nvtOwned = mesh->dofmap_get_own_size(index, _iproc);
    
    // count the ghost node dofs on all levels -------------
    const std::map < unsigned, unsigned > ghostMap = ghost_map_proc(mesh, index);

    // count the total node dofs (own + ghost) on all levels -------------
    unsigned nvt = nvtOwned + ghostMap.size();
    

    unsigned size_conn = size_connectivity_proc(mesh, index);
    
    const unsigned dim_array_coord [] = { nvt * 3 * static_cast<unsigned>(sizeof( float )) };
    const unsigned dim_array_conn[]   = { size_conn * static_cast<unsigned>(sizeof( int )) };
    const unsigned dim_array_off []   = { nel * static_cast<unsigned>(sizeof( int )) };
    const unsigned dim_array_type []  = { nel * static_cast<unsigned>(sizeof( short unsigned )) };
    const unsigned dim_array_reg []   = { nel * static_cast<unsigned>(sizeof( short unsigned )) };
    const unsigned dim_array_elvar [] = { nel * static_cast<unsigned>(sizeof( float )) };
    const unsigned dim_array_ndvar [] = { nvt * static_cast<unsigned>(sizeof( float )) };

    // initialize common buffer_void memory
    unsigned buffer_size = ( dim_array_coord[0] > dim_array_conn[0] ) ? dim_array_coord[0] : dim_array_conn[0];  // maximum size
    void* buffer_void = new char [buffer_size];   /// @todo is this deleted? I don't think so
    char* buffer_char = static_cast <char*>( buffer_void );

    // size of the base64-encoded data
    size_t cch;
    cch = b64::b64_encode( &buffer_char[0], buffer_size , NULL, 0 );
    std::vector <char> enc;
    enc.resize( cch );
    //-----------------------------------------------------------------------------------------------
    
    
    
    //---- NumericVector used for node-based fields -------------------------------------------------------------------------------------------
    NumericVector* num_vec_aux_for_node_fields;
    num_vec_aux_for_node_fields = NumericVector::build().release();

    if( n_processors() == 1 ) { // IF SERIAL
      num_vec_aux_for_node_fields->init( mesh->dofmap_get_dof_offset(index, _nprocs),
                                         mesh->dofmap_get_dof_offset(index, _nprocs), false, SERIAL );
    }
    else { // IF PARALLEL
      num_vec_aux_for_node_fields->init( mesh->dofmap_get_dof_offset(index, _nprocs),
                                         mesh->dofmap_get_own_size(index, _iproc),
                                         mesh->dofmap_get_ghost_dofs(index, _iproc), false, GHOSTED );
    }
    //---- NumericVector used for node-based fields -------------------------------------------------------------------------------------------

    
    //-----------------------------------------------------------------------------------------------
    piece_iproc_begin(fout, nvt, nel);

    // print coordinates *********************************************Solu*******************************************
    fout  << "     <Points>" << std::endl;
    Pfout << "    <PPoints>" << std::endl;
    
            //--------- fill coord ------------
    // point pointer to common memory area buffer of void type;
    float* var_coord = static_cast<float*>( buffer_void );
    
    //print own nodes -------------------------

    unsigned dofOffset = mesh->dofmap_get_dof_offset(index, _iproc);
    
    for( int i = 0; i < 3; i++ ) {
      if( !_surface ) {
        num_vec_aux_for_node_fields->matrix_mult( *mesh->_topology->_Sol[i],
                            *mesh->GetQitoQjProjection( index, 2 ) );
        if( _graph && i == 2 ) {
          unsigned indGraph = _ml_sol->GetIndex( _graphVariable.c_str() );
          num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indGraph],
                              *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indGraph ) ) );
        }
      }
      else {
        unsigned indSurfVar = _ml_sol->GetIndex( _surfaceVariables[i].c_str() );
        num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indSurfVar],
                            *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indSurfVar ) ) );
      }
      for( unsigned ii = 0; ii < nvtOwned; ii++ ) {
        var_coord[ ii * 3 + i] = ( *num_vec_aux_for_node_fields )( ii +  dofOffset );
      }
      if( _ml_sol != NULL && _moving_mesh  && _moving_vars.size() > i){//_ml_mesh->GetLevel( 0 )->GetDimension() > i )  { // if moving mesh
        unsigned indDXDYDZ = _ml_sol->GetIndex( _moving_vars[i].c_str() );
        num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indDXDYDZ],
                            *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indDXDYDZ ) ) );
        for( unsigned ii = 0; ii < nvtOwned; ii++ ) {
          var_coord[ii * 3 + i] += ( *num_vec_aux_for_node_fields )( ii +  dofOffset );
        }
      }
    }
    
    
    //print ghost nodes -------------------------
    unsigned offset_ig = 3 * nvtOwned;

    for( int i = 0; i < 3; i++ ) {
      if( !_surface ) {
        num_vec_aux_for_node_fields->matrix_mult( *mesh-> _topology->_Sol[i],
                            *mesh-> GetQitoQjProjection( index, 2 ) );
        if( _graph && i == 2 ) {
          unsigned indGraphVar = _ml_sol->GetIndex( _graphVariable.c_str() );
          num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indGraphVar],
                              *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indGraphVar ) ) );
        }
      }
      else {
        unsigned indSurfVar = _ml_sol->GetIndex( _surfaceVariables[i].c_str() );
        num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indSurfVar],
                            *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indSurfVar ) ) );
      }
      for( std::map <unsigned, unsigned>::const_iterator it = ghostMap.begin(); it != ghostMap.end(); ++it ) {
        var_coord[ offset_ig + 3 * it->second + i ] = ( *num_vec_aux_for_node_fields )( it->first );
      }
    }

    for( int i = 0; i < 3; i++ ) { // if moving mesh
      if( _ml_sol != NULL && _moving_mesh  && _moving_vars.size() > i ){ //&& mesh->GetDimension() > i )  {
        unsigned indDXDYDZ = _ml_sol->GetIndex( _moving_vars[i].c_str() );
        num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[indDXDYDZ],
                            *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indDXDYDZ ) ) );
        for( std::map <unsigned, unsigned>::const_iterator it = ghostMap.begin(); it != ghostMap.end(); ++it ) {
          var_coord[ offset_ig + 3 * it->second + i ] += ( *num_vec_aux_for_node_fields )( it->first );
        }
      }
    }
            //--------- fill coord ------------


    print_data_array_vector< float >("coordinates", "Float32", 3, fout, Pfout, dim_array_coord, var_coord, enc);
    
    
    fout  << "      </Points>" << std::endl;
    Pfout << "    </PPoints>" << std::endl;
    //-----------------------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------------------
    // Printing of element connectivity - offset - format type  *
    fout  << "      <Cells>" << std::endl;
    Pfout << "    <PCells>" << std::endl;
    
    //-----------------------------------------------------------------------------------------------
    //print connectivity
    int * var_conn = static_cast <int*>( buffer_void );
    
     fill_connectivity_proc(mesh, index, ghostMap, var_conn);

     print_data_array< int >("connectivity", "Int32", fout, Pfout, dim_array_conn, var_conn, enc);


    //-------------------------------------------------------------------------------------------------
    //printing offset
    print_element_based_fields< int >("offsets", "Int32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_off, mesh, index, enc);

    //--------------------------------------------------------------------------------------------------
    //Element format type
    print_element_based_fields< unsigned short >("types", "UInt16", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_type, mesh, index, enc);
    

    fout  << "      </Cells>" << std::endl;
    Pfout << "    </PCells>" << std::endl;
    //--------------------------------------------------------------------------------------------------

    // /Print Cell Data ****************************************************************************
    fout  << "      <CellData Scalars=\"scalars\">" << std::endl;
    Pfout << "    <PCellData Scalars=\"scalars\">" << std::endl;
    

    //------------------------------------------- PARALLEL PARTITION, MATERIAL, GROUP, FE TYPE, LEVEL ---------------------------------------------------------
    
    print_element_based_fields< unsigned short >("Metis partition", "UInt16", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_reg, mesh, index, enc);

    print_element_based_fields< float >("Material", "Float32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, index, enc);

    print_element_based_fields< float >("Group", "Float32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, index, enc);

    print_element_based_fields< float >("TYPE", "Float32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, index, enc);

    print_element_based_fields< float >("Level", "Float32", fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, index, enc);

    


    //------------------------------------------- SOLUTIONS ---------------------------------------------------------
   const bool print_all = print_all_sols(vars);
    
    
    
    if( _ml_sol != NULL ) {
        
      //Print Solution (on elements) ***************************************************************
      for( unsigned i = 0; i < compute_print_sol_size(print_all, vars); i++ ) {
          
        const unsigned solIndex = ( print_all == 0 ) ? _ml_sol->GetIndex( vars[i].c_str() ) : i;
        const unsigned sol_fe_type = _ml_sol->GetSolutionType( solIndex );
        
        if( 3 <= sol_fe_type ) {

          std::string solName =  _ml_sol->GetSolutionName( solIndex );

          for( int name = 0; name < compute_sol_bdc_res_eps_size(solution, i); name++ ) {
              
            std::string printName = print_sol_bdc_res_eps_name(solName, name);
            

            //--------- fill var_el ------------
            float* var_el = static_cast< float*>( buffer_void );
            
            fill_sol_on_elements(mesh, elemetOffset, elemetOffsetp1, solution, name, i,  var_el);
            //--------- fill var_el ------------
            

    print_data_array< float >(printName, "Float32", fout, Pfout, dim_array_elvar, var_el, enc);

              
          }
        }
      }
      
    } //end _ml_sol != NULL

    fout  << "      </CellData>" << std::endl;
    Pfout << "    </PCellData>" << std::endl;
    //   //------------------------------------------------------------------------------------------------

    if( _ml_sol != NULL ) {
        
      fout  << "      <PointData Scalars=\"scalars\"> " << std::endl;
      Pfout << "    <PPointData Scalars=\"scalars\"> " << std::endl;
      
      
      // / Print Solution (on nodes) ********************************************************************
      //Loop on variables

      // point pointer to common memory area buffer of void type;
      float* var_nd = static_cast<float*>( buffer_void );
      
      for( unsigned i = 0; i < compute_print_sol_size(print_all, vars); i++ ) {
          
        const unsigned solIndex = ( print_all == 0 ) ? _ml_sol->GetIndex( vars[i].c_str() ) : i;
        const unsigned sol_fe_type = _ml_sol->GetSolutionType( solIndex );
        
        if( sol_fe_type < 3 ) {
            
          //BEGIN LAGRANGIAN Fem SOLUTION
          std::string solName =  _ml_sol->GetSolutionName( solIndex );

          for( int name = 0; name < compute_sol_bdc_res_eps_size(solution, i); name++ ) {
              
            std::string printName = print_sol_bdc_res_eps_name(solName, name);
           

         //--------- fill var_nd ------------
            //print own dofs -------------------------
            unsigned offset_iprc = mesh->dofmap_get_dof_offset(index, _iproc);
            unsigned nvt_ig = mesh->dofmap_get_own_size(index, _iproc);

            if( name == 0 )
              num_vec_aux_for_node_fields->matrix_mult( *solution->_Sol[solIndex],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( solIndex ) ) );
            else if( name == 1 )
              num_vec_aux_for_node_fields->matrix_mult( *solution->_Bdc[solIndex],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( solIndex ) ) );
            else if( name == 2 )
              num_vec_aux_for_node_fields->matrix_mult( *solution->_Res[solIndex],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( solIndex ) ) );
            else
              num_vec_aux_for_node_fields->matrix_mult( *solution->_Eps[solIndex],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( solIndex ) ) );

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
      
      fout  << "      </PointData>" << std::endl;
      Pfout << "    </PPointData>" << std::endl;
      
    }  //end _ml_sol != NULL

    //------------------------------------------------------------------------------------------------

    piece_iproc_end(fout);

    vtk_unstructured_footer_iproc(fout);
    fout.close();

    vtk_unstructured_footer_parallel_wrapper(Pfout);
    Pfout.close();


    //-----------------------------------------------------------------------------------------------------
    //free memory
    delete num_vec_aux_for_node_fields;

    //--------------------------------------------------------------------------------------------------------
    return;
  }
  
  
  

} //end namespace femus


