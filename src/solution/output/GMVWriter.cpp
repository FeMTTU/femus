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
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include "Files.hpp"


namespace femus {

  GMVWriter::GMVWriter( MultiLevelSolution* ml_sol ) : Writer( ml_sol ) {
    _debugOutput = false;
  }

  GMVWriter::GMVWriter( MultiLevelMesh* ml_mesh ) : Writer( ml_mesh ) {
    _debugOutput = false;
  }

  GMVWriter::~GMVWriter() {

  }

  void GMVWriter::Write( const std::string output_path, const char order[], const std::vector<std::string>& vars, const unsigned time_step ) {

    // ********** linear -> index==0 *** quadratic -> index==1 **********
    unsigned index = ( strcmp( order, "linear" ) ) ? 1 : 0;

    std::string filename_prefix;
    if( _ml_sol != NULL )
      filename_prefix = "sol";
    else
      filename_prefix = "mesh";

    std::ostringstream filename;
    filename << output_path << "/" << filename_prefix << ".level" << _gridn << "." << time_step << "." << order << ".gmv";

    std::ofstream fout;

    if( _iproc != 0 ) {
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

    Mesh* mesh = _ml_mesh->GetLevel( _gridn - 1 );
    Solution* solution = _ml_sol->GetSolutionLevel( _gridn - 1 );
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
    numVector->init( mesh->_dofOffset[index][_nprocs], mesh->_ownSize[index][_iproc], true, AUTOMATIC );


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
      if( !_surface ) {
        numVector->matrix_mult( *mesh->_topology->_Sol[i],
                                *mesh->GetQitoQjProjection( index, 2 ) );
        if( _graph && i == 2 ) {
          unsigned indGraphVar = _ml_sol->GetIndex( _graphVariable.c_str() );
          numVector->matrix_mult( *solution->_Sol[indGraphVar],
                                  *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indGraphVar ) ) );
        }
      }
      else {
        unsigned indSurfVar = _ml_sol->GetIndex( _surfaceVariables[i].c_str() );
        numVector->matrix_mult( *solution->_Sol[indSurfVar],
                                *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indSurfVar ) ) );
      }

      numVector->localize_to_one( vector1, 0 );
      if( _ml_sol != NULL && _moving_mesh  && dim > i )  {
        unsigned indDXDYDZ = _ml_sol->GetIndex( _moving_vars[i].c_str() );
        numVector->matrix_mult( *solution->_Sol[indDXDYDZ],
                                *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( indDXDYDZ ) ) );
        numVector->localize_to_one( vector2, 0 );
        if( _iproc == 0 ) {
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

    //mesh->_topology->_Sol[mesh->GetTypeIndex()]->localize_to_one( vector1, 0 );

    for( unsigned isdom = 0; isdom < _nprocs; isdom++ ) {
      mesh->el->LocalizeElementDof( isdom );
      mesh->el->LocalizeElementQuantities(isdom);
      if( _iproc == 0 ) {
        for( unsigned ii = mesh->_elementOffset[isdom]; ii < mesh->_elementOffset[isdom + 1]; ii++ ) {
          //short unsigned ielt = static_cast < short unsigned >( vector1[ii] + 0.25 );
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
          fout.write( ( char* ) &NVE[ielt][index], sizeof( unsigned ) );
          for( unsigned j = 0; j < NVE[ielt][index]; j++ ) {
            unsigned jnode_Metis = mesh->GetSolutionDof( j, ii, index );
            topology[j] = jnode_Metis + 1;
          }
          fout.write( ( char* ) topology, sizeof( unsigned ) *NVE[ielt][index] );
        }
      }
      mesh->el->FreeLocalizedElementDof();
      mesh->el->FreeLocalizedElementQuantities();
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
    for( int isdom = 0; isdom < _nprocs; isdom++ ) {
      for( unsigned ii = mesh->_elementOffset[isdom]; ii < mesh->_elementOffset[isdom + 1]; ii++ ) {
        vector1[icount] = isdom;
        icount++;
      }
    }
    fout.write( ( char* ) &vector1[0], nel * sizeof( double ) );
    //END METIS PARTITIONING

    //BEGIN SOLUTION
    if( _ml_sol != NULL )  {
      bool printAll = 0;
      for( unsigned ivar = 0; ivar < vars.size(); ivar++ ) {
        printAll += !( vars[ivar].compare( "All" ) ) + !( vars[ivar].compare( "all" ) ) + !( vars[ivar].compare( "ALL" ) );
      }
      for( unsigned ivar = 0; ivar < !printAll * vars.size() + printAll * _ml_sol->GetSolutionSize(); ivar++ ) {
        unsigned i = ( printAll == 0 ) ? _ml_sol->GetIndex( vars[ivar].c_str() ) : ivar;

        for( int name = 0; name < 4; name++ ) {
          //BEGIN LAGRANGIAN Fem SOLUTION
          if( name == 0 ) {
            sprintf( buffer, "%s", _ml_sol->GetSolutionName( i ) );
          }
          else if( name == 1 ) {
            sprintf( buffer, "%s %s", "Bdc", _ml_sol->GetSolutionName( i ) );
          }
          else if( name == 2 ) {
            sprintf( buffer, "%s %s", "Res", _ml_sol->GetSolutionName( i ) );
          }
          else {
            sprintf( buffer, "%s %s", "Eps", _ml_sol->GetSolutionName( i ) );
          }
          if( name == 0 || ( _debugOutput  && solution->_ResEpsBdcFlag[i] ) ) {
            if( _ml_sol->GetSolutionType( i ) < 3 ) { // **********  on the nodes **********
              fout.write( ( char* ) buffer, sizeof( char ) * 8 );
              fout.write( ( char* ) &one, sizeof( unsigned ) );
              if( name == 0 ) {
                numVector->matrix_mult( *solution->_Sol[i], *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( i ) ) );
              }
              else if( name == 1 ) {
                numVector->matrix_mult( *solution->_Bdc[i], *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( i ) ) );
              }
              else if( name == 2 ) {
                numVector->matrix_mult( *solution->_Res[i], *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( i ) ) );
              }
              else {
                numVector->matrix_mult( *solution->_Eps[i], *mesh->GetQitoQjProjection( index, _ml_sol->GetSolutionType( i ) ) );
              }
              numVector->localize_to_one( vector1, 0 );
              fout.write( ( char* ) &vector1[0], nvt * sizeof( double ) );
              //END LAGRANGIAN Fem SOLUTION
            }
            else {
              //BEGIN DISCONTINUOUS Fem SOLUTION
              fout.write( ( char* ) buffer, sizeof( char ) * 8 );
              fout.write( ( char* ) &zero, sizeof( unsigned ) );

              if( name == 0 ) {
                solution->_Sol[i]->localize_to_one( vector2, 0 );
              }
              else if( name == 1 ) {
                solution->_Bdc[i]->localize_to_one( vector2, 0 );
              }
              else if( name == 2 ) {
                solution->_Res[i]->localize_to_one( vector2, 0 );
              }
              else {
                solution->_Eps[i]->localize_to_one( vector2, 0 );
              }
              vector1.resize( nel );
              for( unsigned ii = 0; ii < nel; ii++ ) {
                vector1[ii] = vector2[ mesh->GetSolutionDof( 0, ii, _ml_sol->GetSolutionType( i ) ) ];
              }

              fout.write( ( char* ) &vector1[0], nel * sizeof( double ) );
              //END DISCONTINUOUS Fem SOLUTION
            }
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




