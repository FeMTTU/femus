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
#include <b64/b64.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <iomanip>
#include <algorithm>
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
      
    std::ostringstream level_name_stream;    
    level_name_stream << ".level" << my_level;
    std::string level_name(level_name_stream.str());   
       
    // *********** open vtu files *************
    std::ofstream fout;

    std::string dirnamePVTK = "VTKParallelFiles/";
    Files files;
    files.CheckDir( output_path, "" );
    files.CheckDir( output_path, dirnamePVTK );

    int icount;
    unsigned index = 0;
    if( !strcmp( order, "linear" ) ) 	 index = 0; //linear
    else if( !strcmp( order, "quadratic" ) ) 	 index = 1; //quadratic
    else if( !strcmp( order, "biquadratic" ) ) index = 2; //biquadratic


    std::ostringstream filename;
    filename << output_path << "/" << dirnamePVTK << filename_prefix << level_name << "." << _iproc << "." << time_step << "." << order << suffix_pre_extension << ".vtu";

    fout.open( filename.str().c_str() );
    if( !fout.is_open() ) {
      std::cout << std::endl << " The output file " << filename.str() << " cannot be opened.\n";
      abort();
    }

    // *********** write vtu header ************
    fout << "<?xml version=\"1.0\"?>" << std::endl;
    fout << "<VTKFile type = \"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    fout << "  <UnstructuredGrid>" << std::endl;

    // *********** open pvtu file *************
    std::ofstream Pfout;
    if( _iproc != 0 ) {
      Pfout.rdbuf();   //redirect to dev_null
    }
    else {
      std::ostringstream Pfilename;
      Pfilename << output_path << "/" << filename_prefix << level_name << "." << time_step << "." << order <<  suffix_pre_extension << ".pvtu";
      Pfout.open( Pfilename.str().c_str() );
      if( Pfout.is_open() ) {
        std::cout << std::endl << " The output is printed to file " << Pfilename.str() << " in parallel VTK-XML (64-based) format" << std::endl;
      }
      else {
        std::cout << std::endl << " The output file " << Pfilename.str() << " cannot be opened.\n";
        abort();
      }
    }

    // *********** write pvtu header ***********
    Pfout << "<?xml version=\"1.0\"?>" << std::endl;
    Pfout << "<VTKFile type = \"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    Pfout << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
    for( int jproc = 0; jproc < _nprocs; jproc++ ) {
      Pfout << "    <Piece Source=\"" << dirnamePVTK
            << filename_prefix << level_name << "." << jproc << "." << time_step << "." << order <<  suffix_pre_extension << ".vtu"
            << "\"/>" << std::endl;
    }
    // ****************************************

    Mesh* mesh = _ml_mesh->GetLevel( my_level - 1 );
    Solution* solution = _ml_sol->GetSolutionLevel( my_level - 1 );

    //count the own node dofs on all levels
    unsigned nvt = mesh->_ownSize[index][_iproc];

    // count the ghost node dofs and the own element dofs element on all levels
    unsigned ghostMapCounter = 0;
    map < unsigned, unsigned > ghostMap;
    unsigned counter = 0;
    unsigned elemetOffset = mesh->_elementOffset[_iproc];
    unsigned elemetOffsetp1 = mesh->_elementOffset[_iproc + 1];
    unsigned nel = elemetOffsetp1 - elemetOffset;
    unsigned dofOffset = mesh->_dofOffset[index][_iproc];
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      short unsigned ielt = mesh->GetElementType( iel );
      for( unsigned j = 0; j < mesh->GetElementDofNumber( iel, index ); j++ ) {
        counter++;
        unsigned jdof = mesh->GetSolutionDof( j, iel, index );
        if( jdof < dofOffset ) { // check if jdof is a ghost node
          if( ghostMap.find( jdof ) == ghostMap.end() ) {
            ghostMap[jdof] = ghostMapCounter;
            ghostMapCounter++;
          }
        }
      }
    }

    unsigned nvtOwned = nvt;
    nvt += ghostMap.size(); // total node dofs (own + ghost)

    const unsigned dim_array_coord [] = { nvt * 3 * static_cast<unsigned>(sizeof( float )) };
    const unsigned dim_array_conn[]   = { counter * static_cast<unsigned>(sizeof( int )) };
    const unsigned dim_array_off []   = { nel * static_cast<unsigned>(sizeof( int )) };
    const unsigned dim_array_type []  = { nel * static_cast<unsigned>(sizeof( short unsigned )) };
    const unsigned dim_array_reg []   = { nel * static_cast<unsigned>(sizeof( short unsigned )) };
    const unsigned dim_array_elvar [] = { nel * static_cast<unsigned>(sizeof( float )) };
    const unsigned dim_array_ndvar [] = { nvt * static_cast<unsigned>(sizeof( float )) };

    // initialize common buffer_void memory
    unsigned buffer_size = ( dim_array_coord[0] > dim_array_conn[0] ) ? dim_array_coord[0] : dim_array_conn[0];
    void* buffer_void = new char [buffer_size];
    char* buffer_char = static_cast <char*>( buffer_void );

    size_t cch;
    cch = b64::b64_encode( &buffer_char[0], buffer_size , NULL, 0 );
    vector <char> enc;
    enc.resize( cch );
    char* pt_char;

    fout  << "    <Piece NumberOfPoints= \"" << nvt << "\" NumberOfCells= \"" << nel << "\" >" << std::endl;

    //-----------------------------------------------------------------------------------------------
    // print coordinates *********************************************Solu*******************************************
    fout  << "      <Points>" << std::endl;
    fout  << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;

    Pfout << "    <PPoints>" << std::endl;
    Pfout << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\"/>" << std::endl;

    NumericVector* mysol;
    mysol = NumericVector::build().release();

    if( n_processors() == 1 ) { // IF SERIAL
      mysol->init( mesh->_dofOffset[index][_nprocs],
                   mesh->_dofOffset[index][_nprocs], false, SERIAL );
    }
    else { // IF PARALLEL
      mysol->init( mesh->_dofOffset[index][_nprocs],
                   mesh->_ownSize[index][_iproc],
                   mesh->_ghostDofs[index][_iproc], false, GHOSTED );
    }

    // point pointer to common mamory area buffer of void type;
    float* var_coord = static_cast<float*>( buffer_void );

    for( int i = 0; i < 3; i++ ) {
      if( !_surface ) {
        mysol->matrix_mult( *mesh->_topology->_Sol[i],
                            *mesh->GetQitoQjProjection( index, 2 ) );
        if( _graph && i == 2 ) {
          unsigned indGraph = _ml_sol->GetIndex( _graphVariable.c_str() );
          mysol->matrix_mult( *solution->_Sol[indGraph],
                              *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indGraph ) ) );
        }
      }
      else {
        unsigned indSurfVar = _ml_sol->GetIndex( _surfaceVariables[i].c_str() );
        mysol->matrix_mult( *solution->_Sol[indSurfVar],
                            *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indSurfVar ) ) );
      }
      for( unsigned ii = 0; ii < nvtOwned; ii++ ) {
        var_coord[ ii * 3 + i] = ( *mysol )( ii +  dofOffset );
      }
      if( _ml_sol != NULL && _moving_mesh  && _moving_vars.size() > i){//_ml_mesh->GetLevel( 0 )->GetDimension() > i )  { // if moving mesh
        unsigned indDXDYDZ = _ml_sol->GetIndex( _moving_vars[i].c_str() );
        mysol->matrix_mult( *solution->_Sol[indDXDYDZ],
                            *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indDXDYDZ ) ) );
        for( unsigned ii = 0; ii < nvtOwned; ii++ ) {
          var_coord[ii * 3 + i] += ( *mysol )( ii +  dofOffset );
        }
      }
    }
    unsigned offset_ig = 3 * nvtOwned;

    //print ghost nodes
    for( int i = 0; i < 3; i++ ) {
      if( !_surface ) {
        mysol->matrix_mult( *mesh-> _topology->_Sol[i],
                            *mesh-> GetQitoQjProjection( index, 2 ) );
        if( _graph && i == 2 ) {
          unsigned indGraphVar = _ml_sol->GetIndex( _graphVariable.c_str() );
          mysol->matrix_mult( *solution->_Sol[indGraphVar],
                              *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indGraphVar ) ) );
        }
      }
      else {
        unsigned indSurfVar = _ml_sol->GetIndex( _surfaceVariables[i].c_str() );
        mysol->matrix_mult( *solution->_Sol[indSurfVar],
                            *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indSurfVar ) ) );
      }
      for( std::map <unsigned, unsigned>::iterator it = ghostMap.begin(); it != ghostMap.end(); ++it ) {
        var_coord[ offset_ig + 3 * it->second + i ] = ( *mysol )( it->first );
      }
    }

    for( int i = 0; i < 3; i++ ) { // if moving mesh
      if( _ml_sol != NULL && _moving_mesh  && _moving_vars.size() > i ){ //&& mesh->GetDimension() > i )  {
        unsigned indDXDYDZ = _ml_sol->GetIndex( _moving_vars[i].c_str() );
        mysol->matrix_mult( *solution->_Sol[indDXDYDZ],
                            *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indDXDYDZ ) ) );
        for( std::map <unsigned, unsigned>::iterator it = ghostMap.begin(); it != ghostMap.end(); ++it ) {
          var_coord[ offset_ig + 3 * it->second + i ] += ( *mysol )( it->first );
        }
      }
    }

    cch = b64::b64_encode( &dim_array_coord[0], sizeof( dim_array_coord ), NULL, 0 );
    b64::b64_encode( &dim_array_coord[0], sizeof( dim_array_coord ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    //print coordinates array
    cch = b64::b64_encode( &var_coord[0], dim_array_coord[0] , NULL, 0 );
    b64::b64_encode( &var_coord[0], dim_array_coord[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    fout << std::endl;

    fout  << "        </DataArray>" << std::endl;
    fout  << "      </Points>" << std::endl;
    Pfout << "    </PPoints>" << std::endl;
    //-----------------------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------------------
    // Printing of element connectivity - offset - format type  *
    fout  << "      <Cells>" << std::endl;
    Pfout << "    <PCells>" << std::endl;
    //-----------------------------------------------------------------------------------------------
    //print connectivity
    fout  << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\"/>" << std::endl;

    // point pointer to common mamory area buffer of void type;
    int* var_conn = static_cast <int*>( buffer_void );
    icount = 0;
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      for( unsigned j = 0; j < mesh->GetElementDofNumber( iel, index ); j++ ) {
        unsigned loc_vtk_conn = (mesh->GetElementType( iel ) == 0)? FemusToVTKorToXDMFConn[j] : j;
        unsigned jdof = mesh->GetSolutionDof( loc_vtk_conn, iel, index );
        var_conn[icount] = ( jdof >= dofOffset ) ? jdof - dofOffset : nvtOwned + ghostMap[jdof];
        icount++;
      }
    }

    //print connectivity dimension
    cch = b64::b64_encode( &dim_array_conn[0], sizeof( dim_array_conn ), NULL, 0 );
    b64::b64_encode( &dim_array_conn[0], sizeof( dim_array_conn ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    //print connectivity array
    cch = b64::b64_encode( &var_conn[0], dim_array_conn[0] , NULL, 0 );
    b64::b64_encode( &var_conn[0], dim_array_conn[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    //------------------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------------------
    //printing offset
    fout  << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"Int32\" Name=\"offsets\" format=\"binary\"/>" << std::endl;

    // point pointer to common memory area buffer of void type;
    int* var_off = static_cast <int*>( buffer_void );
    icount = 0;
    int offset_el = 0;
    // print offset array
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      offset_el += mesh->GetElementDofNumber( iel, index );
      var_off[icount] = offset_el;
      icount++;
    }

    //print offset dimension
    cch = b64::b64_encode( &dim_array_off[0], sizeof( dim_array_off ), NULL, 0 );
    b64::b64_encode( &dim_array_off[0], sizeof( dim_array_off ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    //print offset array
    cch = b64::b64_encode( &var_off[0], dim_array_off[0] , NULL, 0 );
    b64::b64_encode( &var_off[0], dim_array_off[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    fout  << std::endl;

    fout  << "        </DataArray>" << std::endl;

    //--------------------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------------------

    //Element format type : 23:Serendipity(8-nodes)  28:Quad9-Biquadratic
    fout  << "        <DataArray type=\"UInt16\" Name=\"types\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"UInt16\" Name=\"types\" format=\"binary\"/>" << std::endl;

    // point pointer to common mamory area buffer of void type;
    unsigned short* var_type = static_cast <unsigned short*>( buffer_void );
    icount = 0;
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      short unsigned ielt = mesh->GetElementType( iel );
      var_type[icount] = femusToVtkCellType[index][ielt];
      icount++;
    }

    //print element format dimension
    cch = b64::b64_encode( &dim_array_type[0], sizeof( dim_array_type ), NULL, 0 );
    b64::b64_encode( &dim_array_type[0], sizeof( dim_array_type ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    //print element format array
    cch = b64::b64_encode( &var_type[0], dim_array_type[0] , NULL, 0 );
    b64::b64_encode( &var_type[0], dim_array_type[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    fout  << std::endl;
    fout  << "        </DataArray>" << std::endl;
    //----------------------------------------------------------------------------------------------------
//
    fout  << "      </Cells>" << std::endl;
    Pfout << "    </PCells>" << std::endl;
    //--------------------------------------------------------------------------------------------------

    // /Print Cell Data ****************************************************************************
    fout  << "      <CellData Scalars=\"scalars\">" << std::endl;
    Pfout << "    <PCellData Scalars=\"scalars\">" << std::endl;

    unsigned short* var_reg = static_cast <unsigned short*>( buffer_void );

    // Print Metis Partitioning
    fout  << "        <DataArray type=\"UInt16\" Name=\"Metis partition\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"UInt16\" Name=\"Metis partition\" format=\"binary\"/>" << std::endl;

    // point pointer to common mamory area buffer of void type;
    unsigned short* var_proc = static_cast <unsigned short*>( buffer_void );

    icount = 0;
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      var_proc[icount] = _iproc;
      icount++;
    }

    //print regions dimension
    cch = b64::b64_encode( &dim_array_reg[0], sizeof( dim_array_reg ), NULL, 0 );
    b64::b64_encode( &dim_array_reg[0], sizeof( dim_array_reg ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    //print regions array
    cch = b64::b64_encode( &var_proc[0], dim_array_reg[0] , NULL, 0 );
    b64::b64_encode( &var_proc[0], dim_array_reg[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    fout  << std::endl;
    fout  << "        </DataArray>" << std::endl;


    //BEGIN SARA&GIACOMO


    //-------------------------------------------MATERIAL---------------------------------------------------------

    //material(fout, Pfout, buffer_void, elemetOffset, elemetOffsetp1, dim_array_elvar, mesh, enc); //TODO

    //-------------------------------------------MATERIAL---------------------------------------------------------

    //NumericVector& material =  mesh->_topology->GetSolutionName( "Material" );

    fout  << "       <DataArray type=\"Float32\" Name=\"" << "Material" << "\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"Float32\" Name=\"" << "Material" << "\" format=\"binary\"/>" << std::endl;
    // point pointer to common memory area buffer of void type;
    float* var_el = static_cast< float*>( buffer_void );
    icount = 0;
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      var_el[icount] = mesh->GetElementMaterial(iel); 
      icount++;
    }
    
    //print solution on element dimension
    cch = b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), NULL, 0 );
    b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    
    //print solution on element array
    cch = b64::b64_encode( &var_el[0], dim_array_elvar[0] , NULL, 0 );
    b64::b64_encode( &var_el[0], dim_array_elvar[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    
    
    
    //------------------------------------------------------GROUP-----------------------------------------------------------

    //NumericVector& group =  mesh->_topology->GetSolutionName( "Group" );

    fout  << "        <DataArray type=\"Float32\" Name=\"" << "Group" << "\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"Float32\" Name=\"" << "Group" << "\" format=\"binary\"/>" << std::endl;
    // point pointer to common memory area buffer of void type;
    var_el = static_cast< float*>( buffer_void );
    icount = 0;
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      var_el[icount] = mesh->GetElementGroup(iel);
      icount++;
    }
    //print solution on element dimension
    cch = b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), NULL, 0 );
    b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    
    //print solution on element array
    cch = b64::b64_encode( &var_el[0], dim_array_elvar[0] , NULL, 0 );
    b64::b64_encode( &var_el[0], dim_array_elvar[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;

    //-------------------------------------------------------TYPE--------------------------------------------------
   // NumericVector& type =  mesh->_topology->GetSolutionName( "Type" );

    fout  << "        <DataArray type=\"Float32\" Name=\"" << "TYPE" << "\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"Float32\" Name=\"" << "TYPE" << "\" format=\"binary\"/>" << std::endl;
    // point pointer to common memory area buffer of void type;
    var_el = static_cast< float*>( buffer_void );
    icount = 0;
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      var_el[icount] = mesh->GetElementType(iel);
      icount++;
    }
    //print solution on element dimension
    cch = b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), NULL, 0 );
    b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    
    //print solution on element array
    cch = b64::b64_encode( &var_el[0], dim_array_elvar[0] , NULL, 0 );
    b64::b64_encode( &var_el[0], dim_array_elvar[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;

    
    //-------------------------------------------------------LEVEL--------------------------------------------------
    fout  << "      <DataArray type=\"Float32\" Name=\"" << "Level" << "\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"Float32\" Name=\"" << "Level" << "\" format=\"binary\"/>" << std::endl;
    // point pointer to common memory area buffer of void type;
    var_el = static_cast< float*>( buffer_void );
    icount = 0;
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      var_el[icount] = mesh->el->GetElementLevel(iel);
      icount++;
    }
    //print solution on element dimension
    cch = b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), NULL, 0 );
    b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    //print solution on element array
    cch = b64::b64_encode( &var_el[0], dim_array_elvar[0] , NULL, 0 );
    b64::b64_encode( &var_el[0], dim_array_elvar[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    
    
    //END SARA&GIACOMO

    if( _ml_sol == NULL ) {
      delete [] var_proc;
    }

    bool print_all = 0;
    for( unsigned ivar = 0; ivar < vars.size(); ivar++ ) {
      print_all += !( vars[ivar].compare( "All" ) ) + !( vars[ivar].compare( "all" ) ) + !( vars[ivar].compare( "ALL" ) );
    }
    if( _ml_sol != NULL ) {
      //Print Solution (on element) ***************************************************************
      for( unsigned i = 0; i < ( !print_all )*vars.size() + print_all * _ml_sol->GetSolutionSize(); i++ ) {
        unsigned solIndex = ( print_all == 0 ) ? _ml_sol->GetIndex( vars[i].c_str() ) : i;
        if( 3 <= _ml_sol->GetSolutionType( solIndex ) ) {

          std::string solName =  _ml_sol->GetSolutionName( solIndex );

          for( int name = 0; name < 1 + 3 * _debugOutput * solution->_ResEpsBdcFlag[i]; name++ ) {
            std::string printName;

            if( name == 0 ) printName = solName;
            else if( name == 1 ) printName = "Bdc" + solName;
            else if( name == 2 ) printName = "Res" + solName;
            else printName = "Eps" + solName;

            fout  << "        <DataArray type=\"Float32\" Name=\"" << printName << "\" format=\"binary\">" << std::endl;
            Pfout << "      <PDataArray type=\"Float32\" Name=\"" << printName << "\" format=\"binary\"/>" << std::endl;
            // point pointer to common memory area buffer of void type;
            float* var_el = static_cast< float*>( buffer_void );
            icount = 0;
            for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
              unsigned iel_Metis = mesh->GetSolutionDof( 0, iel, _ml_sol->GetSolutionType( i ) );
              if( name == 0 )
                var_el[icount] = ( *solution->_Sol[i] )( iel_Metis );
              else if( name == 1 )
                var_el[icount] = ( *solution->_Bdc[i] )( iel_Metis );
              else if( name == 2 )
                var_el[icount] = ( *solution->_Res[i] )( iel_Metis );
              else
                var_el[icount] = ( *solution->_Eps[i] )( iel_Metis );
              icount++;
            }

            //print solution on element dimension
            cch = b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), NULL, 0 );
            b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), &enc[0], cch );
            pt_char = &enc[0];
            for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

            //print solution on element array
            cch = b64::b64_encode( &var_el[0], dim_array_elvar[0] , NULL, 0 );
            b64::b64_encode( &var_el[0], dim_array_elvar[0], &enc[0], cch );
            pt_char = &enc[0];
            for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
            fout << std::endl;
            fout << "        </DataArray>" << std::endl;
          }
        }
      } //end _ml_sol != NULL
    }

    fout  << "      </CellData>" << std::endl;
    Pfout << "    </PCellData>" << std::endl;
    //   //------------------------------------------------------------------------------------------------

    if( _ml_sol != NULL ) {
      //   //------------------------------------------------------------------------------------------------
      // / Print Solution (on nodes) ********************************************************************
      fout  << "      <PointData Scalars=\"scalars\"> " << std::endl;
      Pfout << "    <PPointData Scalars=\"scalars\"> " << std::endl;
      //Loop on variables

      // point pointer to common memory area buffer of void type;
      float* var_nd = static_cast<float*>( buffer_void );
      for( unsigned i = 0; i < ( !print_all )*vars.size() + print_all * _ml_sol->GetSolutionSize(); i++ ) {
        unsigned solIndex = ( print_all == 0 ) ? _ml_sol->GetIndex( vars[i].c_str() ) : i;
        if( _ml_sol->GetSolutionType( solIndex ) < 3 ) {
          //BEGIN LAGRANGIAN Fem SOLUTION
          std::string solName =  _ml_sol->GetSolutionName( solIndex );

          for( int name = 0; name < 1 + 3 * _debugOutput * solution->_ResEpsBdcFlag[i]; name++ ) {
            std::string printName;

            if( name == 0 ) printName = solName;
            else if( name == 1 ) printName = "Bdc" + solName;
            else if( name == 2 ) printName = "Res" + solName;
            else printName = "Eps" + solName;

            fout  << "        <DataArray type=\"Float32\" Name=\"" << printName << "\" format=\"binary\">" << std::endl;
            Pfout << "      <PDataArray type=\"Float32\" Name=\"" << printName << "\" format=\"binary\"/>" << std::endl;

            //print solutions on nodes dimension
            cch = b64::b64_encode( &dim_array_ndvar[0], sizeof( dim_array_ndvar ), NULL, 0 );
            b64::b64_encode( &dim_array_ndvar[0], sizeof( dim_array_ndvar ), &enc[0], cch );
            pt_char = &enc[0];
            for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

            unsigned offset_iprc = mesh->_dofOffset[index][_iproc];
            unsigned nvt_ig = mesh->_ownSize[index][_iproc];

            if( name == 0 )
              mysol->matrix_mult( *solution->_Sol[solIndex],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( solIndex ) ) );
            else if( name == 1 )
              mysol->matrix_mult( *solution->_Bdc[solIndex],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( solIndex ) ) );
            else if( name == 2 )
              mysol->matrix_mult( *solution->_Res[solIndex],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( solIndex ) ) );
            else
              mysol->matrix_mult( *solution->_Eps[solIndex],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( solIndex ) ) );

            for( unsigned ii = 0; ii < nvt_ig; ii++ ) {
              var_nd[ ii ] = ( *mysol )( ii + offset_iprc );
            }
            unsigned offset_ig = nvtOwned;

            for( std::map <unsigned, unsigned>::iterator it = ghostMap.begin(); it != ghostMap.end(); ++it ) {
              var_nd[ offset_ig + it->second ] = ( *mysol )( it->first );
            }

            cch = b64::b64_encode( &var_nd[0], dim_array_ndvar [0], NULL, 0 );
            b64::b64_encode( &var_nd[0], dim_array_ndvar [0], &enc[0], cch );
            pt_char = &enc[0];
            for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
            fout << std::endl;

            fout  << "        </DataArray>" << std::endl;
          }
        } //endif
      } // end for sol
      fout  << "      </PointData>" << std::endl;
      Pfout << "    </PPointData>" << std::endl;
      delete [] var_nd;
    }  //end _ml_sol != NULL

    //------------------------------------------------------------------------------------------------

    fout << "    </Piece>" << std::endl;
    fout << "  </UnstructuredGrid>" << std::endl;
    fout << "</VTKFile>" << std::endl;
    fout.close();

    Pfout << "  </PUnstructuredGrid>" << std::endl;
    Pfout << "</VTKFile>" << std::endl;
    Pfout.close();


    //-----------------------------------------------------------------------------------------------------
    //free memory
    delete mysol;

    //--------------------------------------------------------------------------------------------------------
    return;
  }
  
  
  
  void VTKWriter::material(std::ofstream & fout, std::ofstream & Pfout,
                           void* buffer_void, 
                           const unsigned elemetOffset, const unsigned elemetOffsetp1, 
                           const unsigned * dim_array_elvar, 
                           const Mesh * mesh,
                           std::vector <char> & enc) const { //TODO this function does not work, it prints infty
                               
      //-------------------------------------------MATERIAL---------------------------------------------------------

    //NumericVector& material =  mesh->_topology->GetSolutionName( "Material" );

    fout  << "       <DataArray type=\"Float32\" Name=\"" << "Material" << "\" format=\"binary\">" << std::endl;
    Pfout << "      <PDataArray type=\"Float32\" Name=\"" << "Material" << "\" format=\"binary\"/>" << std::endl;
    // point pointer to common memory area buffer of void type;
    float* var_el = static_cast< float*>( buffer_void );
    int icount = 0;
    for( int iel = elemetOffset; iel < elemetOffsetp1; iel++ ) {
      var_el[icount] = mesh->GetElementMaterial(iel); 
      //std::cout<< var_el[icount] <<" ";
      icount++;
    }

    //print solution on element dimension
    size_t cch = b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), NULL, 0 );
    b64::b64_encode( &dim_array_elvar[0], sizeof( dim_array_elvar ), &enc[0], cch );
    char* pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    
    //print solution on element array
    cch = b64::b64_encode( &var_el[0], dim_array_elvar[0] , NULL, 0 );
    b64::b64_encode( &var_el[0], dim_array_elvar[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;

  }

} //end namespace femus


