/*=========================================================================

 Program: FEMUS
 Module: GMVWriter
 Authors: Eugenio Aulisa, Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "GMVWriter.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "FElemTypeEnum_list.hpp"
#include "Files.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstring>


namespace femus {

  GMVWriter::GMVWriter(const MultiLevelSolution* ml_sol ) : Writer( ml_sol ) {  }

  GMVWriter::GMVWriter(const MultiLevelMesh* ml_mesh ) : Writer( ml_mesh ) {  }

 

  void GMVWriter::Write( const std::string output_path,
                         const std::string order,
                         const std::vector<std::string>& vars,
                         const unsigned time_step ) {
  
    const Solution * solution = get_solution(_gridn);

    const std::string filename_prefix = _writer_one_level.get_filename_prefix(solution);
    
    const std::string suffix_pre_extension = "";
    
    Write(_gridn, filename_prefix, output_path, suffix_pre_extension, order, vars, time_step);

  }
  
  
  void GMVWriter::Write(const std::string filename_prefix, 
                        const std::string output_path,
                        const std::string order,
                        const std::vector<std::string>& vars, 
                        const unsigned time_step ) {

    const std::string suffix_pre_extension = "";

    Write(_gridn, filename_prefix, output_path, suffix_pre_extension, order, vars, time_step);

 }
  
  
  
  void GMVWriter::Write(const unsigned level_in, 
                        const std::string filename_prefix, 
                        const std::string output_path,
                        const std::string suffix_pre_extension, 
                         const std::string order,
                         const std::vector<std::string>& vars, 
                         const unsigned time_step ) {

    //------------- Mesh, def - BEGIN ----------------------------------------------------------------------------------
    const Mesh * mesh = get_mesh(level_in);
    //------------- Mesh, def - END ----------------------------------------------------------------------------------
    
    //------------- Solution, def - BEGIN ----------------------------------------------------------------------------------
    const Solution * solution = get_solution(level_in);
    //------------- Solution, def - END ----------------------------------------------------------------------------------


    //------------- File stream - BEGIN ----------------------------------------------------------------------------------
    std::ostringstream filename;
    filename << output_path << "/" << filename_prefix << ".level" << level_in << "." << time_step << "." << order << suffix_pre_extension << ".gmv";

    std::ofstream fout;

    if( _writer_one_level.processor_id() != 0 ) {
      fout.rdbuf();   //redirect to dev_null
    }
    else {
      fout.open( filename.str().c_str() );
      if( fout.is_open() ) {
        std::cout << std::endl << " The output is printed to file " << filename.str() << " in GMV format" << std::endl;
      }
      else {
        std::cout << std::endl << " The output file " << filename.str() << " cannot be opened.\n";
        abort();
      }
    }
    //------------- File stream - END ----------------------------------------------------------------------------------


    // ********** linear -> index==0 *** quadratic -> index==1 **********
    const unsigned index = ( strcmp( order.c_str(), fe_fams_for_files[ FILES_CONTINUOUS_LINEAR ].c_str() ) ) ? FILES_CONTINUOUS_QUADRATIC : FILES_CONTINUOUS_LINEAR;


    unsigned nvt = mesh->GetTotalNumberOfDofs( index );
    unsigned nel = mesh->GetNumberOfElements();
    unsigned dim = mesh->GetDimension();

    std::vector < double > vector1;
    vector1.reserve( nvt );

    std::vector < double > vector2;
    if( nvt > ( dim + 1 ) *nel )
      vector2.reserve( nvt );
    else
      vector2.reserve( ( dim + 1 ) *nel );

    NumericVector* numVector;
    numVector = NumericVector::build().release();
    numVector->init( mesh->dofmap_get_dof_offset(index, _writer_one_level.n_processors() ), mesh->dofmap_get_own_size(index, _writer_one_level.processor_id()), true, AUTOMATIC );


    //BEGIN GMV FILE PRINT
    char* buffer = new char[10];
    sprintf( buffer, "%s", "gmvinput" );
    fout.write( ( char* ) buffer, sizeof( char ) * 8 );
    sprintf( buffer, "%s", "ieeei4r8" );
    fout.write( ( char* ) buffer, sizeof( char ) * 8 );


    //BEGIN COORDINATES
    sprintf( buffer, "%s", "nodes" );
    fout.write( ( char* ) buffer, sizeof( char ) * 8 );
    fout.write( ( char* ) &nvt, sizeof( unsigned ) );

    for( int i = 0; i < 3; i++ ) {
      if( ! _writer_one_level.is_surface() ) {
        numVector->matrix_mult( *mesh->GetTopology()->_Sol[i],   * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, 2, * mesh) );
        if( _writer_one_level._graph && i == 2 ) {
          const unsigned indGraphVar = solution->GetIndex( _writer_one_level._graphVariable.c_str() );
          numVector->matrix_mult( *solution->_Sol[indGraphVar],  * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indGraphVar ), * mesh ) );
        }
      }
      else if ( _writer_one_level.is_surface() && solution != NULL ) {
        const unsigned indSurfVar = solution->GetIndex( _writer_one_level._surfaceVariables[i].c_str() );
        numVector->matrix_mult( *solution->_Sol[indSurfVar],  * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indSurfVar ), * mesh ) );
      }

      numVector->localize_to_one( vector1, 0 );
      if( solution != NULL && _writer_one_level._moving_mesh  && dim > i )  {
        const unsigned indDXDYDZ = solution->GetIndex( _writer_one_level._moving_vars[i].c_str() );
        numVector->matrix_mult( *solution->_Sol[indDXDYDZ], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( indDXDYDZ ), * mesh ) );
        numVector->localize_to_one( vector2, 0 );
        if( _writer_one_level.processor_id() == 0 ) {
          for( unsigned i = 0; i < nvt; i++ )
            vector1[i] += vector2[i];
        }
      }
      fout.write( ( char* ) &vector1[0], nvt * sizeof( double ) );

    }
    //END COORDINATES

    //BEGIN CONNETTIVITY
    const int eltp[2][6] = {{8, 4, 6, 4, 3, 2}, {20, 10, 15, 8, 6, 3}};
    sprintf( buffer, "%s", "cells" );
    fout.write( ( char* ) buffer, sizeof( char ) * 8 );

    fout.write( ( char* ) &nel, sizeof( unsigned ) );

    unsigned topology[27];

    //mesh->GetTopology()->_Sol[mesh->GetTypeIndex()]->localize_to_one( vector1, 0 );

    for( unsigned isdom = 0; isdom < _writer_one_level.n_processors(); isdom++ ) {
      mesh->GetMeshElements()->LocalizeElementDof( isdom );
      mesh->GetMeshElements()->LocalizeElement_Level_Type_Group_Material(isdom);
      if( _writer_one_level.processor_id() == 0 ) {
        for( unsigned ii = mesh->GetElementOffset(isdom); ii < mesh->GetElementOffset(isdom + 1); ii++ ) {
	  short unsigned ielt = mesh->GetElementType(ii);
          if( ielt == 0 )
            sprintf( buffer, "phex%d", eltp[index][0] );
          else if( ielt == 1 )
            sprintf( buffer, "ptet%d", eltp[index][1] );
          else if( ielt == 2 )
            sprintf( buffer, "pprism%d", eltp[index][2] );
          else if( ielt == 3 ) {
            if( eltp[index][3] == 8 )
              sprintf( buffer, "%dquad", eltp[index][3] );
            else
              sprintf( buffer, "quad" );
          }
          else if( ielt == 4 ) {
            if( eltp[index][4] == 6 )
              sprintf( buffer, "%dtri", eltp[index][4] );
            else
              sprintf( buffer, "tri" );
          }
          else if( ielt == 5 ) {
            if( eltp[index][5] == 3 )
              sprintf( buffer, "%dline", eltp[index][5] );
            else
              sprintf( buffer, "line" );
          }
          fout.write( ( char* ) buffer, sizeof( char ) * 8 );
          const unsigned nvertices = mesh->GetMeshElements()->GetNVE(ielt, index); 
          fout.write( ( char* ) &  nvertices, sizeof( unsigned ) );
          for( unsigned j = 0; j < nvertices; j++ ) {
            unsigned jnode_Metis = mesh->GetSolutionDof( j, ii, index );
            topology[j] = jnode_Metis + 1;
          }
          fout.write( ( char* ) topology, sizeof( unsigned ) *  nvertices  );
        }
      }
      mesh->GetMeshElements()->FreeLocalizedElementDof();
      mesh->GetMeshElements()->FreeLocalizedElement_Level_Type_Group_Material();
    }

    //END CONNETTIVITY

    //BEGIN VARIABLES
    const unsigned zero = 0u;
    const unsigned one = 1u;
    sprintf( buffer, "%s", "variable" );
    fout.write( ( char* ) buffer, sizeof( char ) * 8 );

    //BEGIN METIS PARTITIONING
    strcpy( buffer, "METIS_DD" );
    fout.write( ( char* ) buffer, sizeof( char ) * 8 );
    fout.write( ( char* ) &zero, sizeof( unsigned ) );

    vector1.resize( nel );
    unsigned icount = 0;
    for( int isdom = 0; isdom < _writer_one_level.n_processors(); isdom++ ) {
      for( unsigned ii = mesh->GetElementOffset(isdom); ii < mesh->GetElementOffset(isdom + 1); ii++ ) {
        vector1[icount] = isdom;
        icount++;
      }
    }
    fout.write( ( char* ) &vector1[0], nel * sizeof( double ) );
    //END METIS PARTITIONING

    //BEGIN SOLUTION
    if( solution != NULL )  {

      bool printAll = 0;
      for( unsigned ivar = 0; ivar < vars.size(); ivar++ ) {
        printAll += !( vars[ivar].compare( "All" ) ) + !( vars[ivar].compare( "all" ) ) + !( vars[ivar].compare( "ALL" ) );
      }

      for( unsigned ivar = 0; ivar < !printAll * vars.size() + printAll * solution->GetSolutionSize(); ivar++ ) {
        unsigned i = ( printAll == 0 ) ? solution->GetIndex( vars[ivar].c_str() ) : ivar;

        for( int name = 0; name < 4; name++ ) {
          
            const std::string printName = Writer_one_level::print_sol_bdc_res_eps_name( solution->GetSolName_from_index( i ) , name);

          
          if( name == Writer_one_level::_index_sol || ( _writer_one_level._debugOutput  && solution->is_unknown_of_system(i) ) ) {
            
            
              //BEGIN LAGRANGIAN Fem SOLUTION
            if( solution->GetSolutionType( i ) < NFE_FAMS_C_ZERO_LAGRANGE ) { // **********  on the nodes **********
              fout.write( ( char* ) printName.c_str(), sizeof( char ) * 8 );
              fout.write( ( char* ) &one, sizeof( unsigned ) );
              if( name == Writer_one_level::_index_sol ) {
                numVector->matrix_mult( *solution->_Sol[i], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( i ), * mesh ) );
              }
              else if( name == Writer_one_level::_index_bdc ) {
                numVector->matrix_mult( *solution->_Bdc[i], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( i ), * mesh ) );
              }
              else if( name == Writer_one_level::_index_res ) {
                numVector->matrix_mult( *solution->_Res[i], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( i ), * mesh ) );
              }
              else if( name == Writer_one_level::_index_eps ) {
                numVector->matrix_mult( *solution->_Eps[i], * _writer_one_level._fe_proj_matrices.GetQitoQjProjection( index, solution->GetSolutionType( i ), * mesh ) );
              }
              numVector->localize_to_one( vector1, 0 );
              fout.write( ( char* ) &vector1[0], nvt * sizeof( double ) );
            }
              //END LAGRANGIAN Fem SOLUTION

              //BEGIN DISCONTINUOUS Fem SOLUTION
            else if( solution->GetSolutionType( i ) < NFE_FAMS ) { // **********  on the elements **********
              fout.write( ( char* )  printName.c_str(), sizeof( char ) * 8 );
              fout.write( ( char* ) &zero, sizeof( unsigned ) );

              if( name == Writer_one_level::_index_sol ) {
                solution->_Sol[i]->localize_to_one( vector2, 0 );
              }
              else if( name == Writer_one_level::_index_bdc ) {
                solution->_Bdc[i]->localize_to_one( vector2, 0 );
              }
              else if( name == Writer_one_level::_index_res ) {
                solution->_Res[i]->localize_to_one( vector2, 0 );
              }
              else if( name == Writer_one_level::_index_eps ) {
                solution->_Eps[i]->localize_to_one( vector2, 0 );
              }
              vector1.resize( nel );
              for( unsigned ii = 0; ii < nel; ii++ ) {
                vector1[ii] = vector2[ mesh->GetSolutionDof( 0, ii, solution->GetSolutionType( i ) ) ];
              }

              fout.write( ( char* ) &vector1[0], nel * sizeof( double ) );
            }
              //END DISCONTINUOUS Fem SOLUTION

          }
          
        }
      }
    }
    //END SOLUTION

    sprintf( buffer, "%s", "endvars" );
    fout.write( ( char* ) buffer, sizeof( char ) * 8 );
    //END VARIABLES

    sprintf( buffer, "%s", "endgmv" );
    fout.write( ( char* ) buffer, sizeof( char ) * 8 );
    fout.close();
    //END GMV FILE PRINT

    delete numVector;
    delete [] buffer;

    return;
  }

} //end namespace femus




