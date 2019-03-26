/*=========================================================================

 Program: FEMUS
 Module: XDMFWriter
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
#include "FemusConfig.hpp"
#include "XDMFWriter.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <iomanip>

#ifdef HAVE_HDF5
#include "hdf5.h"
#endif

#include "MultiLevelMeshTwo.hpp"
#include "DofMap.hpp"
#include "SystemTwo.hpp"
#include "NumericVector.hpp"
#include "GeomElemBase.hpp"
#include "FETypeEnum.hpp"
#include "paral.hpp"

namespace femus {

  const std::string XDMFWriter::type_el[3][N_GEOM_ELS] = {{"Hexahedron", "Tetrahedron", "Wedge", "Quadrilateral", "Triangle", "Polyline"}, //linear
    {"Hexahedron_20", "Tetrahedron_10", "Not_implemented", "Quadrilateral_8", "Triangle_6", "Edge_3"},        //serendipity
    {"Hexahedron_27", "Tetrahedron_10", "Not_implemented", "Quadrilateral_9", "Triangle_6", "Edge_3"}
  };            //tensor-product quadratic (some real, some fake)
  /// @todo Tri6 and Tet10 are actually SERENDIPITY, not TENSOR-PRODUCT QUADRATIC.
  // The corresponding tensor-product ones should be Tri7 and Tet14

  const std::string XDMFWriter::_nodes_name = "/NODES";
  const std::string XDMFWriter::_elems_name = "/ELEMS";
  const std::string XDMFWriter::_conn = "CONN";

  XDMFWriter::XDMFWriter( MultiLevelSolution* ml_sol ) : Writer( ml_sol ) {
    _debugOutput = false;
  }

  XDMFWriter::XDMFWriter( MultiLevelMesh* ml_mesh ) : Writer( ml_mesh ) {
    _debugOutput = false;
  }

  XDMFWriter::~XDMFWriter() {}

  void XDMFWriter::Write( const std::string output_path, const char order[], const std::vector<std::string>& vars, const unsigned time_step ) {

#ifdef HAVE_HDF5

    bool print_all = 0;
    for( unsigned ivar = 0; ivar < vars.size(); ivar++ ) {
      print_all += !( vars[ivar].compare( "All" ) ) + !( vars[ivar].compare( "all" ) ) + !( vars[ivar].compare( "ALL" ) );
    }

    unsigned index_nd = 0;
    if( !strcmp( order, "linear" ) ) {   //linear
      index_nd = 0;
    }
    else if( !strcmp( order, "quadratic" ) ) {   //quadratic
      index_nd = 1;
    }
    else if( !strcmp( order, "biquadratic" ) ) {   //tensor-product quadratic (real and fake)
      index_nd = 2;
    }

    Mesh* mesh = _ml_mesh->GetLevel( _gridn - 1 );
    Solution* solution = _ml_sol->GetSolutionLevel( _gridn - 1 );

    /// @todo I assume that the mesh is not mixed
    std::string type_elem;
    unsigned iel0 = mesh->_elementOffset[_iproc];
    unsigned elemtype = mesh->GetElementType( iel0 );

    type_elem = XDMFWriter::type_el[index_nd][elemtype];

    if( type_elem.compare( "Not_implemented" ) == 0 )
    {
      std::cerr << "XDMF-Writer error: element type not supported!" << std::endl;
      abort();
    }

    unsigned nvt = mesh->_dofOffset[index_nd][_nprocs];
    unsigned nel = mesh->GetNumberOfElements();
    unsigned dim = mesh->GetDimension();
    unsigned maxDim = ( nvt > ( dim + 1 ) * nel ) ? nvt : ( dim + 1 ) * nel;
    unsigned ndofs = mesh->el->GetNVE( elemtype, index_nd );

    std::vector < int > var_conn( nel * ndofs );

    std::vector < double > vector1;
    vector1.reserve( nvt );

    std::vector < double > vector2;
    vector2.reserve( maxDim );

    NumericVector* numVector = NumericVector::build().release();
    numVector->init( nvt, mesh->_ownSize[index_nd][_iproc], true, AUTOMATIC );

    //BEGIN XMF FILE PRINT

    std::string filename_prefix;
    if( _ml_sol != NULL ) filename_prefix = "sol";
    else filename_prefix = "mesh";

    std::ostringstream xdmf_filename;
    xdmf_filename << output_path << "/" << filename_prefix << ".level" << _gridn << "." << time_step << "." << order << ".xmf";

    std::ofstream fout;

    if( _iproc != 0 ) {
      fout.rdbuf();   //redirect to dev_null
    }
    else {
      fout.open( xdmf_filename.str().c_str() );
      if( fout.is_open() ) {
        std::cout << std::endl << " The output is printed to file " << xdmf_filename.str() << " in XDMF-HDF5 format" << std::endl;
      }
      else {
        std::cout << std::endl << " The output file " << xdmf_filename.str() << " cannot be opened.\n";
        abort();
      }
    }

    // The HDF5 file name
    std::ostringstream hdf5_filename;
    std::ostringstream hdf5_filename2;
    hdf5_filename2 << filename_prefix << ".level" << _gridn << "." << time_step << "." << order << ".h5";
    hdf5_filename << output_path << "/" <<  hdf5_filename2.str();

    // head ************************************************
    fout << "<?xml version=\"1.0\" ?>" << std::endl;
    fout << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd []\">" << std::endl;
    fout << "<Xdmf>" << std::endl;
    fout << "<Domain>" << std::endl;
    fout << "<Grid Name=\"Mesh\">" << std::endl;
    fout << "<Time Value =\"" << time_step << "\" />" << std::endl;
    fout << "<Topology Type=\"" << type_elem << "\" Dimensions=\"" << nel << "\">" << std::endl;
    //Connectivity
    fout << "<DataStructure DataType=\"Int\" Dimensions=\"" << nel << " " << ndofs << "\"" << "  Format=\"HDF\">" << std::endl;
    fout << hdf5_filename2.str() << ":/CONNECTIVITY" << std::endl;
    fout << "</DataStructure>" << std::endl;
    fout << "</Topology>" << std::endl;
    fout << "<Geometry Type=\"X_Y_Z\">" << std::endl;
    //Node_X
    fout << "<DataStructure DataType=\"Double\" Precision=\"8\" Dimensions=\"" << nvt << "  1\"" << "  Format=\"HDF\">" << std::endl;
    fout << hdf5_filename2.str() << ":/NODES_X1" << std::endl;
    fout << "</DataStructure>" << std::endl;
    //Node_Y
    fout << "<DataStructure DataType=\"Double\" Precision=\"8\" Dimensions=\"" << nvt << "  1\"" << "  Format=\"HDF\">" << std::endl;
    fout << hdf5_filename2.str() << ":/NODES_X2" << std::endl;
    fout << "</DataStructure>" << std::endl;
    //Node_Z
    fout << "<DataStructure DataType=\"Double\" Precision=\"8\" Dimensions=\"" << nvt << "  1\"" << "  Format=\"HDF\">" << std::endl;
    fout << hdf5_filename2.str() << ":/NODES_X3" << std::endl;
    fout << "</DataStructure>" << std::endl;
    fout << "</Geometry>" << std::endl;
    //Metis partitions
    fout << "<Attribute Name=\"" << "Domain_partitions" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
    fout << "<DataItem DataType=\"Double\" Dimensions=\"" << nel << "  1\""  << "  Format=\"HDF\">" << std::endl;
    fout << hdf5_filename2.str() << ":/DOMAIN_PARTITIONS" << std::endl;
    fout << "</DataItem>" << std::endl;
    fout << "</Attribute>" << std::endl;

    if( _ml_sol != NULL ) {
      // Solution Variables
      for( unsigned i = 0; i < ( 1 - print_all ) * vars.size() + print_all * _ml_sol->GetSolutionSize(); i++ ) {
        unsigned indx = ( print_all == 0 ) ? _ml_sol->GetIndex( vars[i].c_str() ) : i;
        //Printing biquadratic solution on the nodes
        if( _ml_sol->GetSolutionType( indx ) < 3 ) {
          std::string solName =  _ml_sol->GetSolutionName( indx );
          for( int name = 0; name < 1 + 3 * _debugOutput * solution->_ResEpsBdcFlag[i]; name++ ) {
            std::string printName;
            if( name == 0 ) printName = solName;
            else if( name == 1 ) printName = "Bdc" + solName;
            else if( name == 2 ) printName = "Res" + solName;
            else printName = "Eps" + solName;
            fout << "<Attribute Name=\"" << printName << "\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
            fout << "<DataItem DataType=\"Double\" Precision=\"8\" Dimensions=\"" << nvt << "  1\"" << "  Format=\"HDF\">" << std::endl;
            fout << hdf5_filename2.str() << ":/" << printName << std::endl;
            fout << "</DataItem>" << std::endl;
            fout << "</Attribute>" << std::endl;
          }
        }
        else if( _ml_sol->GetSolutionType( indx ) >= 3 ) {   //Printing picewise constant solution on the element
          std::string solName =  _ml_sol->GetSolutionName( indx );
          for( int name = 0; name < 1 + 3 * _debugOutput * solution->_ResEpsBdcFlag[i]; name++ ) {
            std::string printName;
            if( name == 0 ) printName = solName;
            else if( name == 1 ) printName = "Bdc" + solName;
            else if( name == 2 ) printName = "Res" + solName;
            else printName = "Eps" + solName;
            fout << "<Attribute Name=\"" << printName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
            fout << "<DataItem DataType=\"Double\" Precision=\"8\" Dimensions=\"" << nel << "  1\"" << "  Format=\"HDF\">" << std::endl;
            fout << hdf5_filename2.str() << ":/" << printName << std::endl;
            fout << "</DataItem>" << std::endl;
            fout << "</Attribute>" << std::endl;
          }
        }
      }
    }
    //end ml_sol

    fout << "</Grid>" << std::endl;
    fout << "</Domain>" << std::endl;
    fout << "</Xdmf>" << std::endl;
    fout.close();
    //END XMF FILE PRINT

    //BEGIN HD5 FILE PRINT
    hid_t file_id;
    if( _iproc == 0 ) file_id = H5Fcreate( hdf5_filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    hsize_t dimsf[2];
    herr_t status;
    hid_t dataspace;
    hid_t dataset;

    //BEGIN COORDINATES
    for( int i = 0; i < 3; i++ ) {
      numVector->matrix_mult( *mesh->_topology->_Sol[i], *mesh->GetQitoQjProjection( index_nd, 2 ) );
      numVector->localize_to_one( vector1, 0 );

      if( _ml_sol != NULL && _moving_mesh && dim > i ) {
        unsigned varind_DXDYDZ = _ml_sol->GetIndex( _moving_vars[i].c_str() );
        numVector->matrix_mult( *solution->_Sol[varind_DXDYDZ],
                                *mesh->GetQitoQjProjection( index_nd, _ml_sol->GetSolutionType( varind_DXDYDZ ) ) );
        numVector->localize_to_one( vector2, 0 );
        if( _iproc == 0 ) {
          for( unsigned ii = 0; ii < nvt; ii++ ) vector1[ii] += vector2[ii];
        }
      }

      if( _iproc == 0 ) {
        dimsf[0] = nvt ;
        dimsf[1] = 1;
        std::ostringstream Name;
        Name << "/NODES_X" << i + 1;
        dataspace = H5Screate_simple( 2, dimsf, NULL );
        dataset   = H5Dcreate( file_id, Name.str().c_str(), H5T_NATIVE_DOUBLE,
                               dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vector1[0] );
        H5Sclose( dataspace );
        H5Dclose( dataset );
      }
    } //end 3d loop
    //END COORDINATES

    //BEGIN CONNETTIVITY
    unsigned icount = 0;
    for( unsigned isdom = 0; isdom < _nprocs; isdom++ ) {
      mesh->el->LocalizeElementDof( isdom );
      if( _iproc == 0 ) {
        for( unsigned iel = mesh->_elementOffset[isdom]; iel < mesh->_elementOffset[isdom + 1]; iel++ ) {
          for( unsigned j = 0; j < ndofs; j++ ) {
            unsigned vtk_loc_conn = FemusToVTKorToXDMFConn[j];
            var_conn[icount] = mesh->GetSolutionDof( vtk_loc_conn, iel, index_nd );
            icount++;
          }
        }
      }
      mesh->el->FreeLocalizedElementDof();
    }

    if(_iproc == 0) {
      dimsf[0] = nel * ndofs ;
      dimsf[1] = 1;
      dataspace = H5Screate_simple( 2, dimsf, NULL );
      dataset   = H5Dcreate( file_id, "/CONNECTIVITY", H5T_NATIVE_INT,
                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      status   = H5Dwrite( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var_conn[0] );
      H5Sclose( dataspace );
      H5Dclose( dataset );
    }

    //END CONNETTIVITY

    //BEGIN METIS PARTITIONING
    if( _iproc == 0 ) {
      unsigned icount = 0;
      for( int isdom = 0; isdom < _nprocs; isdom++ ) {
        for( unsigned ii = mesh->_elementOffset[isdom]; ii < mesh->_elementOffset[isdom + 1]; ii++ ) {
          vector1[icount] = isdom;
          icount++;
        }
      }
      dimsf[0] = nel;
      dimsf[1] = 1;
      dataspace = H5Screate_simple( 2, dimsf, NULL );
      dataset   = H5Dcreate( file_id, "/DOMAIN_PARTITIONS", H5T_NATIVE_DOUBLE,
                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      status   = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vector1[0] );
      H5Sclose( dataspace );
      H5Dclose( dataset );
    }
    //END METIS PARTITIONING

    //BEGIN SOLUTION
    if( _ml_sol != NULL )  {
      //BEGIN DISCONTINUOUS Fem SOLUTION
      for( unsigned i = 0; i < ( 1 - print_all ) *vars.size() + print_all * _ml_sol->GetSolutionSize(); i++ ) {
        unsigned indx = ( print_all == 0 ) ? _ml_sol->GetIndex( vars[i].c_str() ) : i;
        if( _ml_sol->GetSolutionType( indx ) >= 3 ) {
          for( int name = 0; name < 1 + 3 * _debugOutput * solution->_ResEpsBdcFlag[i]; name++ ) {

            std::string solName =  _ml_sol->GetSolutionName( indx );
            std::string printName;
            if( name == 0 ) {
              solution->_Sol[indx]->localize_to_one( vector2, 0 );
              printName = solName;
            }
            else if( name == 1 ) {
              solution->_Bdc[indx]->localize_to_one( vector2, 0 );
              printName = "Bdc" + solName;
            }
            else if( name == 2 ) {
              solution->_Res[indx]->localize_to_one( vector2, 0 );
              printName = "Res" + solName;
            }
            else {
              solution->_Eps[indx]->localize_to_one( vector2, 0 );
              printName = "Eps" + solName;
            }

            vector1.resize( nel );
            for( unsigned ii = 0; ii < nel; ii++ ) {
              vector1[ii] = vector2[ mesh->GetSolutionDof( 0, ii, _ml_sol->GetSolutionType( indx ) ) ];
            }

            if( _iproc == 0 ) {
              dimsf[0] = nel;
              dimsf[1] = 1;
              dataspace = H5Screate_simple( 2, dimsf, NULL );
              dataset   = H5Dcreate( file_id, printName.c_str(), H5T_NATIVE_DOUBLE,
                                     dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
              status   = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vector1[0] );
              H5Sclose( dataspace );
              H5Dclose( dataset );
            }
          }
        }
      }
      //END DISCONTINUOUS Fem SOLUTION

      //BEGIN LAGRANGIAN Fem SOLUTION
      for( unsigned i = 0; i < ( 1 - print_all ) *vars.size() + print_all * _ml_sol->GetSolutionSize(); i++ ) {
        unsigned indx = ( print_all == 0 ) ? _ml_sol->GetIndex( vars[i].c_str() ) : i;
        if( _ml_sol->GetSolutionType( indx ) < 3 ) {
          for( int name = 0; name < 1 + 3 * _debugOutput * solution->_ResEpsBdcFlag[i]; name++ ) {

            std::string solName =  _ml_sol->GetSolutionName( indx );
            std::string printName;
            if( name == 0 ) {
              numVector->matrix_mult( *solution->_Sol[indx],
                                      *mesh->GetQitoQjProjection( index_nd, _ml_sol->GetSolutionType( indx ) ) );
              printName = solName;
            }
            else if( name == 1 ) {
              numVector->matrix_mult( *solution->_Bdc[indx],
                                      *mesh->GetQitoQjProjection( index_nd, _ml_sol->GetSolutionType( indx ) ) );
              printName = "Bdc" + solName;
            }
            else if( name == 2 ) {
              numVector->matrix_mult( *solution->_Res[indx],
                                      *mesh->GetQitoQjProjection( index_nd, _ml_sol->GetSolutionType( indx ) ) );
              printName = "Res" + solName;
            }
            else {
              numVector->matrix_mult( *solution->_Eps[indx],
                                      *mesh->GetQitoQjProjection( index_nd, _ml_sol->GetSolutionType( indx ) ) );
              printName = "Eps" + solName;
            }

            numVector->localize_to_one( vector1, 0 );

            if( _iproc == 0 ) {
              dimsf[0] = nvt;
              dimsf[1] = 1;
              dataspace = H5Screate_simple( 2, dimsf, NULL );
              dataset   = H5Dcreate( file_id, printName.c_str(), H5T_NATIVE_DOUBLE,
                                     dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
              status   = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vector1[0] );
              H5Sclose( dataspace );
              H5Dclose( dataset );
            }
          }
        }
      }
      //END LAGRANGIAN Fem SOLUTION
    }
    //END SOLUTION

    if( _iproc == 0 ) H5Fclose( file_id );

    //END HD5 FILE PRINT

    delete numVector;

#endif

    return;
  }

  void XDMFWriter::write_solution_wrapper( const std::string output_path, const char type[] ) const {

#ifdef HAVE_HDF5

    // to add--> _time_step0, _ntime_steps
    int time_step0 = 0;
    int ntime_steps = 1;
    int print_step = 1;

    // Print The Xdmf transient wrapper
    std::ostringstream filename;
    filename << output_path << "/sol.level" << _gridn << "." << type << ".xmf";

    std::ofstream ftr_out;
    ftr_out.open( filename.str().c_str() );
    if( !ftr_out ) {
      std::cout << "Transient Output mesh file " << filename.str().c_str() << " cannot be opened.\n";
      exit( 0 );
    }

    ftr_out << "<?xml version=\"1.0\" ?> \n";
    ftr_out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd []\">" << std::endl;
    ftr_out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> " << std::endl;
    ftr_out << "<Domain> " << std::endl;
    ftr_out << "<Grid Name=\"Mesh\" GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
    // time loop for grid sequence
    for( unsigned time_step = time_step0; time_step < time_step0 + ntime_steps; time_step++ ) {
      if( !( time_step % print_step ) ) {
        filename << output_path << "/sol.level" << _gridn << "." << time_step << "." << type << ".xmf";
        ftr_out << "<xi:include href=\"" << filename.str().c_str() << "\" xpointer=\"xpointer(//Xdmf/Domain/Grid[" << 1 << "])\">\n";
        ftr_out << "<xi:include href=\"" << filename.str() << "\" xpointer=\"xpointer(//Xdmf/Domain/Grid[" << 1 << "])\">\n";
        ftr_out << "<xi:fallback/>\n";
        ftr_out << "</xi:include>\n";
      }
    }
    ftr_out << "</Grid> \n";
    ftr_out << "</Domain> \n";
    ftr_out << "</Xdmf> \n";
    ftr_out.close();
    ftr_out.close();
    //----------------------------------------------------------------------------------------------------------

#endif

  }



// =================================================================
//this function is ok here because it doesn't involve the map
//of the EQUATIONS, it is just a print of a time sequence to a .xmf file
/// Xdmf transient  print
//print from time t_idx_in to t_idx_final
//I will print a separate time sequence for each LEVEL
//TODO see if there is a way to read multiple sol collections, defined on different grids, in the SAME time file
  void XDMFWriter::transient_print_xmf( const std::string output_path, const uint t_idx_in, const uint t_idx_final, const int print_step, const uint nolevels_in ) {

    //multigrid
    uint NoLevels = nolevels_in;

    // time parameters
    const uint ndigits     = DEFAULT_NDIGITS;

    // dir names
    std::string    basetime     = DEFAULT_BASETIME;
    std::string    ext_xdmf     = DEFAULT_EXT_XDMF;
    std::string    basesol      = DEFAULT_BASESOL;
    std::string    aux_xdmf     = DEFAULT_AUX_XDMF;

// =================================
// ============= LEVELS ============
// =================================

    for( uint l = 0; l < NoLevels; l++ ) {

      // file
      std::ostringstream Name;
      Name << output_path << "/" << basetime << "."
           << std::setw( ndigits ) << std::setfill( '0' ) << t_idx_in << "-"
           << std::setw( ndigits ) << std::setfill( '0' ) << t_idx_final  << "_l" << l
           << ext_xdmf;
      std::ofstream out( Name.str().c_str() );
      if( out.fail() ) {
        std::cout << "transient_print_xmf: cannot print timeN.xmf" << std::endl;
        abort();
      }

      uint nprt = 1;
      std::string gname[3];
      gname[0] = basesol;

      //   Mesh -----------------------------------
      out << "<?xml version=\"1.0\" ?> \n";
      out << "<!DOCTYPE Xdmf SYSTEM "
          <<  "\"" << output_path << "/" << aux_xdmf << "\"" << "[]>\n";
      out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
      out << "<Domain> \n";
      for( uint kp = 0; kp < nprt; kp++ )    {
        out << "<Grid Name=\"" << gname[kp].c_str() << "\"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
        // time loop for grid sequence
        for( uint it = t_idx_in; it <= t_idx_final; it++ ) if( it % print_step == 0 )   {
            out << "<xi:include href=\""
                << basesol << "." << std::setw( ndigits ) << std::setfill( '0' ) << it << "_l" << l <<  ext_xdmf
                << "\"" << " xpointer=\"xpointer(//Xdmf/Domain/Grid[" << kp + 1 << "])\" >\n";
            out << "<xi:fallback />\n";
            out << " </xi:include>\n";
          }
        out << "</Grid> \n";
      }
      // Grid Collection end
      out << "</Domain> \n";
      out << "</Xdmf> \n";
      out.close();

    } //end levels

    return;
  }



// =============================================================
  hid_t XDMFWriter::read_Dhdf5( hid_t file, const std::string& name, double* data ) {

    hid_t  dataset = H5Dopen( file, name.c_str(), H5P_DEFAULT );
    hid_t status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
    H5Dclose( dataset );
    return status;
  }

// ===========================================================================
  hid_t XDMFWriter::read_Ihdf5( hid_t file, const std::string& name, int* data ) {

#ifdef HAVE_HDF5

    hid_t  dataset = H5Dopen( file, name.c_str(), H5P_DEFAULT );
    hid_t status = H5Dread( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
    H5Dclose( dataset );
    return status;

#endif

  }

// ===========================================================================
  hid_t XDMFWriter::read_UIhdf5( hid_t file, const std::string& name, uint* data ) {

#ifdef HAVE_HDF5

    hid_t  dataset = H5Dopen( file, name.c_str(), H5P_DEFAULT );
    hid_t status = H5Dread( dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
    H5Dclose( dataset );
    return status;

#endif

  }

//=========================================
  hid_t XDMFWriter::print_Dhdf5( hid_t file, const std::string& name, hsize_t dimsf[], double* data ) {

#ifdef HAVE_HDF5

    hid_t dataspace = H5Screate_simple( 2, dimsf, NULL );
    hid_t dataset = H5Dcreate( file, name.c_str(), H5T_NATIVE_DOUBLE,
                               dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    hid_t  status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
    H5Dclose( dataset );
    H5Sclose( dataspace );
    return status;

#endif

  }

/// Print int data into dhdf5 file
//TODO can we make a TEMPLATE function that takes either "double" or "int" and uses
//either H5T_NATIVE_DOUBLE or H5T_NATIVE_INT? Everything else is the same
  hid_t XDMFWriter::print_Ihdf5( hid_t file, const std::string& name, hsize_t dimsf[], int* data ) {

#ifdef HAVE_HDF5

    hid_t dataspace = H5Screate_simple( 2, dimsf, NULL );
    hid_t dataset = H5Dcreate( file, name.c_str(), H5T_NATIVE_INT,
                               dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    hid_t  status = H5Dwrite( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
    H5Dclose( dataset );
    H5Sclose( dataspace );
    return status;

#endif

  }

//H5T_NATIVE_UINT
  hid_t XDMFWriter::print_UIhdf5( hid_t file, const std::string& name, hsize_t dimsf[], uint* data ) {

#ifdef HAVE_HDF5

    hid_t dataspace = H5Screate_simple( 2, dimsf, NULL );
    hid_t dataset = H5Dcreate( file, name.c_str(), H5T_NATIVE_UINT,
                               dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    hid_t  status = H5Dwrite( dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
    H5Dclose( dataset );
    H5Sclose( dataspace );
    return status;

#endif

  }


  void XDMFWriter::PrintXDMFAttribute( std::ofstream& outstream,
                                       std::string hdf5_filename,
                                       std::string hdf5_field,
                                       std::string attr_name,
                                       std::string attr_type,
                                       std::string attr_center,
                                       std::string data_type,
                                       int data_dim_row,
                                       int data_dim_col
                                     ) {

    outstream << "<Attribute Name=\"" << attr_name << "\" "
              << " AttributeType=\"" << attr_type << "\" "
              << " Center=\"" <<  attr_center << "\">\n";
    outstream << "<DataItem  DataType=\"" << data_type << "\" "
              << " Precision=\"" << 8 << "\" "
              << " Dimensions=\"" << data_dim_row << " " << data_dim_col << "\" "
              << " Format=\"HDF\">  \n";
    outstream << hdf5_filename  << ":" << hdf5_field << "\n";
    outstream << "</DataItem>\n" << "</Attribute>\n";

    return;
  }


  void XDMFWriter::PrintXDMFTopology( std::ofstream& outfstream,
                                      std::string hdf5_file,
                                      std::string hdf5_field,
                                      std::string top_type,
                                      int top_dim,
                                      int datadim_n_elems,
                                      int datadim_el_nodes
                                    ) {

    outfstream << "<Topology Type="
               << "\"" <<  top_type  << "\" "
               << "Dimensions=\"" << top_dim << "\"> \n";
    outfstream << "<DataStructure DataType= \"Int\""
               << " Dimensions=\" "  <<  datadim_n_elems << " " <<  datadim_el_nodes
               << "\" Format=\"HDF\">  \n";
    outfstream << hdf5_file << ":" << hdf5_field << " \n";
    outfstream << "</DataStructure> \n" << "</Topology> \n";

    return;
  }


  void XDMFWriter::PrintXDMFGeometry( std::ofstream& outfstream,
                                      std::string hdf5_file,
                                      std::string hdf5_field,
                                      std::string coord_lev,
                                      std::string geom_type,
                                      std::string data_type,
                                      int data_dim_one,
                                      int data_dim_two ) {

    outfstream << "<Geometry Type=\"" << geom_type << "\"> \n";
    for( uint ix = 1; ix < 4; ix++ ) {
      outfstream << "<DataStructure DataType=\"" << data_type
                 << "\" Precision=\"8\" "
                 << "Dimensions=\"" << data_dim_one << " " << data_dim_two
                 << "\" Format=\"HDF\">  \n";
      outfstream << hdf5_file << ":" << hdf5_field << ix << coord_lev << "\n";
      outfstream << "</DataStructure> \n";
    }
    outfstream << " </Geometry>\n";

    return;
  }






// ============================================================
/// This function prints the solution: for quad and linear fem
/// the mesh is assumed to be quadratic, so we must print on quadratic nodes
/// first we print QUADRATIC variables over quadratic nodes, straightforward
/// then, LINEAR variables over quadratic nodes --> we interpolate
/// The NODE numbering of our mesh is so that LINEAR NODES COME FIRST
/// So, to interpolate we do the following:
/// First, we put the values of the LINEAR (coarse) MESH in the FIRST POSITIONS of the sol vector
/// Then, we loop over the elements:
///     we pick the element LINEAR values FROM THE FIRST POSITIONS of the SOL VECTOR itself,
///     we multiply them by the PROLONGATOR,
///     we put the result in the sol vector
/// Therefore the sol values are REWRITTEN every time and OVERWRITTEN because of adjacent elements.
// Ok, allora dobbiamo fare in modo che questa routine stampi Quadratici, Lineari e Costanti
// Ricordiamoci che tutto viene fatto in un mesh "RAFFINATO UNA VOLTA IN PIU' RISPETTO AL LIVELLO PIU' FINE"
// per quanto riguarda gli  elementi, noi abbiamo un valore per ogni elemento,
//e poi dovremo trasferirlo ai suoi figli originati dal raffinamento
// quindi dovremo avere una mappa che per ogni elemento ci da' i suoi FIGLI originati da un raffinamento,
// il tutto ovviamente seguendo gli ordinamenti delle connettivita' di FEMuS.
// Allora, per quanto riguarda gli elementi, viene fatta una mesh_conn_lin che si occupa di mostrare gli
// elementi del mesh linearizzato
//questa connettivita' linearizzata riguarda soltanto il livello FINE
// finora SIA le VARIABILI QUADRATICHE sia le VARIABILI LINEARI
//sono stampate sul "MESH FINE QUADRATICO, reso come MESH LINEARIZZATO"
// per le variabili COSTANTI, esse le stampero' sulla CELLA, che ovviamente sara' una CELLA LINEARIZZATA,
// la stessa che uso per il resto, il mio file soluzione .xmf credo che possa utilizzare UNA SOLA TIPOLOGIA.
// Quindi dobbiamo appoggiarci alla stessa topologia "raffinata linearizzata".
// Per fare questo penso che possiamo copiare quello che fa il file case.h5!

//The printing of the constant will be very similar to the printing of the PID for the case;
// so we could make a common routine which means "print_cell_property_on_linearized_mesh"

//By the way, observe that with multigrid and with the problems about printing in XDMF format with Paraview
//we are in practice having ONE COARSER LEVEL, for LINEAR COARSE NODES,
// and ONE FINER LEVEL for LINEAR NODES, given by this "auxiliary further mesh refinement"

//TODO Parallel: is this function called only by one proc I hope...
//TODO That information should stay INSIDE a FUNCTION, not outside

//TODO this function works only with Level=NoLevels - 1 !
//facciamola funzionare per ogni livello.
//ora che stampo le connettivita' a tutti i livelli per il mesh LINEARIZZATO, posso fare come voglio.
//Quindi stampo su MSH_CONN di un certo livello
//Ok, dobbiamo distinguere come si esplora la node_dof map e la lunghezza di ogni livello
//
// ===============================================
//An idea could be: can we make the projection of any FE solution onto this "REFINED LINEARIZED" MESH
// an automatic process, by using the Prol operators?

//Now, the behaviour of this function should be somewhat "parallel" to the construction of the node_dof.

//TODO the group must be REOPENED by EACH EQUATION.. forse Gopen anziche' Gcreate
// I need to have a group that is first created, then closed and reopened without being emptied...
// it should do the same as a File does

//ok, so far we will print the variables with _LEVEL... now we have to fix the wrapper as well

//Ok, now that we make it print things for ALL variables at ALL levels, we need to JUMP on the _node_dof map.
// TODO the NODE DOF NUMBERING is "CONSECUTIVE" ONLY ON THE FINE LEVEL... and only for QUADRATIC variables, 	for (uint i=0;i< n_nodes;i++)
// which is the same order as the mesh!
//So only on the FINE level you can loop over NODES, in all the other cases
// you need to loop over ELEMENTS
//For all the other levels we have to JUMP... so we will jump the same way as when we build the node dof!
//Now for all variables, loop over all subdomains, collect all the values, and print them...
// I cannot print them all together

//TODO TODO TODO here it is more complicated... pos_sol_ivar will sum up to the QUADRATIC NODES,
//but here we pick the quadratic nodes MORE THAN ONCE, so how can we count them EXACTLY?!!
//the problem is that when you are not on the fine level you have jumps, so, in order to count the LINEAR nodes,
//you need a flag for every quadratic node that tells you if it is also linear or not...
//we should either build a flag field in advance, so that we don't do it now, or maybe later is ok
//since we don't keep any extra info about each node (we do not have a Node class, or an Elem class),
//now it's time we have some flags...
//i guess i should loop over all the quadratic nodes, and let the flag for the linear only
// ok so when i pick the linear nodes i may set the flag is linear there
//flag = 1 means: it is linear, otherwise 0

// // //       for (int i = mesh->_off_nd[QQ][off_proc];
// // // 	       i < mesh->_off_nd[QQ][off_proc+Level+1];i++) {
// // //       if (i< mesh->_off_nd[QQ][off_proc]
// // // 	   + mesh->_off_nd[LL][off_proc + Level+1 ]
// // // 	   - mesh->_off_nd[LL][off_proc])
// // //
// // // 	 }
//now, in the quadratic nodes, what are the positions of the linear nodes?
//if I am not wrong, the first nodes in the COARSE LEVEL in every processor
//are the LINEAR ONES, that's why the offsets are what they are...
//in fact, we are looping on the GEomElObjects of each level.

//     int* flag_is_linear = new int[n_nodes]; // we are isolating MESH NODES correspoding to LINEAR DOFS
//what, but we already know from the node numbering what nodes are linear, because we distinguished them, right?!
//     for (uint i=0;i< n_nodes;i++) flag_is_linear[i] = 0;

//LINEAR VARIABLES
//these are re-dimensionalized //so you dont need to multiply below! //the offset for _node_dof is always quadratic, clearly, because the mesh is quadratic
//        //I am filling two QUADRATIC arrays first by the LINEAR POSITIONS
//so pay attention that if you are not setting to ZERO for the linear case, but exactly replacing at the required points

//Always remember that in order to pick the DofValues you have to provide the DofObjects in the correct manner

// QUADRATIC VARIABLES
// pos_in_mesh_obj gives me the position of the GEomElObject: in fact the quadratic dofs are built on the quadratic GeomEls exactly in this order
//here we are picking the NODES per subd and level, so we are sure don't pass MORE TIMES on the SAME NODE

// PRoblem with the linear variables in write and PrintBc.
//There is one line that brings to mistake, but TWO different mistakes.
// in PrintBc it seems to be related to HDF5;
// in write it seems to concern PETSC
//so that is the wrong line, if i comment it everything seems to work for any processor.
//with two and three processors it seems to give even different errors...
//when there is an error related to HDF5, the stack has the _start thing...
// with two procs there is an error related to Petsc,
// with three procs there is an error related to HDF5...

//If I only use PrintBc and not write,
//both with 2 and 3 processors the errors are related to HDF5... this is so absolutely weird...

//Now it seems like I am restricted to that line. That line is responsible for the error.
//Then, what is wrong with that? I would say: I am doing sthg wrong on the array positions.
//Let me put a control on the position. maybe the _el_map is somehow ruined

//Quando hai un segmentation fault, devi concentrarti non tanto sui VALORI ASSEGNATI
// quanto sugli INDICI DEGLI ARRAY !!!

// to do the interpolation over the fine mesh you loop over the ELEMENTS
//So of course you will pass through most nodes MORE THAN ONCE
//after setting correctly the linear nodes, then the quadratic ones are done by projection, without any problem
//for every element, I take the linear nodes.
//Then I loop over the quadratic nodes, and the value at each of them is the prolongation
// of the linear values.
//So, for every quadratic node, I compute the value and set the value
//now, the connectivity is that of the quadratic mesh, but with respect
//to the FINE node numbering
//now we have to convert from the fine node numbering to the node numbering at Level!!!
//do we already have something to do this in the MESH?!?
//given a quadratic node in FINE NUMBERING, obtained by an element AT LEVEL L,
//can we obtain the position of that node AT LEVEL L
//don't we have the connectivities AT LEVEL L,
//where the numbering goes from ZERO to n_nodes_level?
//I would say it is exactly the _node_dof FOR THE FIRST QUADRATIC VARIABLE!
//TODO OF COURSE THAT WOULD IMPLY THAT AT LEAST ONE QUADRATIC VARIABLE is BUILT
//TODO Also, I have to check that COARSER GEOMETRIES are ASSOCIATED to COARSER TOPOLOGIES!
// _node_dof in serial is the IDENTITY, BUT NOT IN PARALLEL!!
//Per passare dal Qnode in FINE NUMBERING al Qnode in SERIAL NUMBERING
//non basta passarlo alla _node_dof del livello, perche' quando sei in parallelo
//lui conta prima i dof quadratici, poi quelli lineari, poi quelli costanti, e bla bla bla...
//e quindi bisogna cambiare il modo di pigliarli!
//Si' bisogna usare qualcosa di SIMILE alla NODE_DOF, ma NON la node_dof,
//perche' quella, per dato livello, conta i dof QQ,LL,KK del proc0, poi QQ,LL,KK del proc1, and so on...
//OK OK OK! Direi che quello che dobbiamo usare e' proprio la node_map, quella che leggiamo dal mesh.h5!!!
//Quindi _node_dof e' per TUTTI i DOF,
// mentre _node_map e' solo per i NODI GEOMETRICI!
//mi sa pero' che la _node_map fa l'opposto di quello che vogliamo noi, cioe'
//dato il nodo di un certo livello ti restituisce il nodo FINE
//Noi invece abbiamo il NODO FINE, e vogliamo avere il nodo di un certo livello.
//Questo si otterrebbe leggendo TUTTA la map e non comprimendola come facciamo...
//oppure io la ricostruisco rifacendo il loop della node_dof ma solo per UNA variabile quadratica!
//praticamente, MA NON DEL TUTTO!, la _node_map e' l'inverso della _node_dof!!!
//TODO ma non c'e' nessun altro posto nel codice in cui devi passare
//da QNODI FINI a QNODI di LIVELLO?
// sembra di no, piu' che altro devi passare da QNODI FINI a DOF di LIVELLO!
//ora qui, essendo che dobbiamo stampare su griglie di diverso livello,
//e siccome le CONNETTIVITA' di TUTTI I LIVELLI SONO DATE RISPETTO all'unico FINE NUMBERING
// (TODO beh, volendo potrei utilizzare le connettivita' "di livello" che stampo nel file msh_conn_lin...
//Il fatto e' che quelle non sono variabile di classe (in Mesh) e quindi mi limito a calcolare, stampare e distruggere...)
// ALLORA DEVO FARE IL PASSAGGIO da QNODI FINI a QNODI DI LIVELLO.
//Questo me lo costruisco io domattina.
//Siccome e' gia' stato calcolato nel gencase, allora mi conviene evitare di fare questo calcolo.
//Quando leggo il vettore dal gencase, posso costruire un vettore di PAIRS...
//ovviamente e' piu' lento, perche' deve cercare un elemento con degli if...
//allora faccio un array 2x, per cui calcolo io l'indice da estrarre...

// We have to find the LINEAR NODES positions in the QUADRATIC LIST
// the QUADRATIC list at each level is based on
// the "LINEARIZED" CONNECTIVITIES of that level
//    plus the LIST OF COORDINATES of that level

//Now there is another mistake, in printing the LINEAR CONNECTIVITIES!
//The fact is always that _el_map gives the Qnodes in FINE NUMBERING,
//while I want the Qnode numbering AT EACH LEVEL.
//because at coarser levels I have fewer nodes, so the list of coords is shorter!

//Ok, the printing of the linearized connectivities is WRONG at NON-FINE LEVEL

//TODO AAAAAAAAAAAAAAAAA: Well, the thing is this: the connectivities of any level must be expressed
// in terms of the node numbering OF THAT LEVEL!
//if the connectivities are expressed with the FINE NODE NUMBERING,
//then we always have to convert to the LEVEL NODE NUMBERING!!!
//I guess the COORDINATES are WRONG...
//Ok, first of all we already have the connectivities at all levels
// of quadratic elements in mesh.h5, in FINE NODE NUMBERING.
//We only have to convert them to LEVEL NODE NUMBERING.

//Ok, I fixed the COORDINATES of the Qnodes at EVERY LEVEL.
//Plus, I have the "LINEARIZED" CONNECTIVITIES at EVERY LEVEL.
//Therefore, it seems like I can print any vector at EVERY LEVEL.
//Now, I still have a problem for the LINEAR variables.
//I guess I'm putting them in the WRONG PLACES before the interpolations.

//allora, quando faccio i loop con _off_nd, quei numeri li' corrispondono alle posizioni nel vettore sol?
//TODO chi da' l'ordine dei nodi del mesh A CIASCUN LIVELLO?
// Le COORDINATE DEL MESH A CIASCUN LIVELLO, quello e' l'ordine di riferimento!

//ok, the position i corresponds to the FINE MESH... then you must translate it for the LEVEL mesh

//For the linear variables we have TWO PROC LOOPS

//Now I guess I have to pick the position of the linear nodes on the mesh by using the ExtendedLevel on the map...
//In the following interpolation the loop is on the elements. From the elements you pick the connectivities,
//          which yield you the FINE QUADRATIC NODES, from which you pick the node numbers at level


//  PRINT OF SOL0000
// the first sol, and the case, are printed BEFORE the INITIAL CONDITIONS are set,
// or the initial conditions are set only at the FINE LEVEL?
// I would say the first one, otherwise I would expect different solutions...
// The first sol and the case are supposed to be equal...
//No wait, what i am saying is not true, because the first sol and the case
// at the FINE LEVEL are printed as they should,
// after the initial conditions,
//but not on COARSER LEVELS!
// Ok now the initial conditions are set only AT THE FINE LEVEL.
// Depending on the kind of multigrid cycle, you may want to set the initial conditions
// ONLY AT THE FINE LEVEL or AT ALL LEVELS...
// For us, let us just do the Initialize function in such a way that all the levels can be treated separately.
// then, if we need it, we call it for ALL LEVELS, or we call it for ONLY THE FINE, or ONLY THE COARSE, or whatever...

//TODO ok, we have to remember one basic principle about our multigrid algorithm:
//the true solution is contained only at the FINE LEVEL!
//at all the other levels we have the DELTAS
//so, the level where to pick the values is the FINE LEVEL.
//Then, of course, we will print at all levels, but the VALUES from x_old are taken from the FINE LEVEL
//no... but wait a second... i am printing at all levels, so that's fine! I wanna print the RESIDUAL for all levels,
//except for the fine level where i print the true solution

// This prints All Variables of One Equation
  void XDMFWriter::write( const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn ) {

    std::vector<GeomElemBase*> fe_in( QL );
    for( int fe = 0; fe < QL; fe++ )    fe_in[fe] = GeomElemBase::build( mesh->_geomelem_id[mesh->get_dim() - 1 - VV].c_str(), fe );


    hid_t file_id = H5Fopen( namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

    // ==========================================
    // =========== FOR ALL LEVELS ===============
    // ==========================================
    for( uint Level = 0; Level < mesh->_NoLevels; Level++ )  {

      std::ostringstream grname;
      grname << "LEVEL" << Level;
//      hid_t group_id = H5Gcreate(file_id, grname.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//      hid_t group_id = H5Gopen(file_id, grname.str().c_str(),H5P_DEFAULT);

      int NGeomObjOnWhichToPrint[QL];
      NGeomObjOnWhichToPrint[QQ] = mesh->_NoNodesXLev[Level];
      NGeomObjOnWhichToPrint[LL] = mesh->_NoNodesXLev[Level];
      NGeomObjOnWhichToPrint[KK] = mesh->_n_elements_vb_lev[VV][Level] * NRE[mesh->_eltype_flag[VV]];

      const uint n_nodes_lev = mesh->_NoNodesXLev[Level];
      double* sol_on_Qnodes  = new double[n_nodes_lev];  //TODO VALGRIND //this is QUADRATIC because it has to hold  either quadratic or linear variables and print them on a QUADRATIC mesh

      // ===================================
      // ========= QUADRATIC ===============
      // ===================================
      for( uint ivar = 0; ivar < dofmap->_nvars[QQ]; ivar++ )        {

        int pos_in_mesh_obj = 0;
        for( uint isubdom = 0; isubdom < mesh->_NoSubdom; isubdom++ ) {
          uint off_proc = isubdom * mesh->_NoLevels;

          for( int fine_node = mesh->_off_nd[QQ][off_proc];
               fine_node < mesh->_off_nd[QQ][off_proc + Level + 1]; fine_node++ ) {

            int pos_in_sol_vec_lev = dofmap->GetDof( Level, QQ, ivar, fine_node );
            int pos_on_Qnodes_lev  = mesh->_Qnode_fine_Qnode_lev[Level][ fine_node ];

#ifndef NDEBUG
            if( pos_on_Qnodes_lev >= ( int ) n_nodes_lev ) {
              std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl;
              abort();
            }
#endif
            sol_on_Qnodes[ pos_on_Qnodes_lev/* pos_in_mesh_obj*/ ] = ( * eqn->_LinSolver[Level]->_EPSC )( pos_in_sol_vec_lev ) * eqn->_refvalue[ ivar + dofmap->_VarOff[QQ] ];
            pos_in_mesh_obj++;
          }
        }  //end subd

#ifndef NDEBUG
        if( pos_in_mesh_obj != NGeomObjOnWhichToPrint[QQ] ) {
          std::cout << "Wrong counting of quadratic nodes" << std::endl;
          abort();
        }
#endif

        std::ostringstream var_name;
        var_name << eqn->_var_names[ ivar + dofmap->_VarOff[QQ] ] << "_" << grname.str(); 	 //         std::string var_name = grname.str() + "/" + _var_names[ivar];
        hsize_t  dimsf[2];
        dimsf[0] = NGeomObjOnWhichToPrint[QQ];
        dimsf[1] = 1;
        XDMFWriter::print_Dhdf5( file_id, var_name.str(), dimsf, sol_on_Qnodes );  //TODO VALGRIND

      }

      // =================================
      // ========= LINEAR ================
      // =================================
      uint elnds[QL_NODES];
      elnds[QQ] = mesh->_elnodes[VV][QQ];
      elnds[LL] = mesh->_elnodes[VV][LL];
      double* elsol_c = new double[elnds[LL]];

      for( uint ivar = 0; ivar < dofmap->_nvars[LL]; ivar++ )        {

//               for (uint i=0; i< n_nodes_lev; i++) { sol_on_Qnodes[i] = 0.; }

        for( uint isubdom = 0; isubdom < mesh->_NoSubdom; isubdom++ ) {
          uint off_proc = isubdom * mesh->_NoLevels;
          for( int fine_node = mesh->_off_nd[QQ][off_proc];
               fine_node < mesh->_off_nd[QQ][off_proc] +
               mesh->_off_nd[LL][off_proc + Level + 1 ]
               - mesh->_off_nd[LL][off_proc]; fine_node++ ) {

            int pos_in_sol_vec_lev = dofmap->GetDof( Level, LL, ivar, fine_node );
            int pos_on_Qnodes_lev = mesh->_Qnode_fine_Qnode_lev[Level][ fine_node ];

#ifndef NDEBUG
            if( pos_in_sol_vec_lev == -1 ) {
              std::cout << "Not correct DOF number at required level" << std::endl;
              abort();
            }
            if( pos_on_Qnodes_lev >= ( int ) n_nodes_lev ) {
              std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl;
              abort();
            }
#endif

            sol_on_Qnodes[ pos_on_Qnodes_lev ] = ( *eqn->_LinSolver[Level]->_EPSC )( pos_in_sol_vec_lev ) * eqn->_refvalue[ ivar + dofmap->_VarOff[LL] ];

          }
        }

        //  2bB element interpolation over the fine mesh -----------------------
        // the way you filled linear positions before completely affects what happens next, which is only geometric
        for( uint iproc = 0; iproc < mesh->_NoSubdom; iproc++ ) {
          uint off_proc = iproc * mesh->_NoLevels;
          int iel_b = mesh->_off_el[VV][off_proc + Level];
          int iel_e = mesh->_off_el[VV][off_proc + Level + 1];

          for( int iel = 0; iel < ( iel_e - iel_b ); iel++ ) {

            for( uint in = 0; in < elnds[LL]; in++ ) {
              int pos_Qnode_fine = mesh->_el_map[VV][( iel + iel_b ) * elnds[QQ] + in ];
              int pos_Qnode_lev  = mesh->_Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];
              elsol_c[in] = sol_on_Qnodes[ pos_Qnode_lev ];   /**_refvalue[ivar]*/ //Do not multiply here!
            }

            for( uint in = 0; in < elnds[QQ]; in++ ) {  //TODO this loop can be done from elnds[LL] instead of from 0
              double sum = 0.;
              for( uint jn = 0; jn < elnds[LL]; jn++ ) {
                sum += fe_in[LL]->get_prol( in * elnds[LL] + jn ) * elsol_c[jn];
              }

              int pos_Qnode_fine = mesh->_el_map[VV][( iel + iel_b ) * elnds[QQ] + in ];  //Qnode in FINE NUMBERING
              int pos_Qnode_lev  = mesh->_Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];  //Qnode in Level NUMBERING

#ifndef NDEBUG
              if( pos_Qnode_lev == -1 ) {
                std::cout << "Not correct node number at required level" << std::endl;
                abort();
              }
              if( pos_Qnode_lev >= ( int ) n_nodes_lev ) {
                std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl;
                abort();
              }
#endif

              sol_on_Qnodes[ pos_Qnode_lev ] = sum;
            }
          }
        } // 2bB end interpolation over the fine mesh --------

        std::ostringstream var_name;
        var_name << eqn->_var_names[ ivar + dofmap->_VarOff[LL] ] << "_" << grname.str();
        hsize_t  dimsf[2];
        dimsf[0] = NGeomObjOnWhichToPrint[LL];
        dimsf[1] = 1;
        XDMFWriter::print_Dhdf5( file_id, var_name.str(), dimsf, sol_on_Qnodes );

      } // ivar linear

      delete []elsol_c;
      delete []sol_on_Qnodes;

      // ===================================
      // ========= CONSTANT ================
      // ===================================
      double* sol_on_cells;
      sol_on_cells = new double[ NGeomObjOnWhichToPrint[KK] ];

      for( uint ivar = 0; ivar < dofmap->_nvars[KK]; ivar++ )        {

        int cel = 0;
        for( uint iproc = 0; iproc < mesh->_NoSubdom; iproc++ ) {
          uint off_proc = iproc * mesh->_NoLevels;

          int sum_elems_prev_sd_at_lev = 0;
          for( uint pr = 0; pr < iproc; pr++ ) {
            sum_elems_prev_sd_at_lev += mesh->_off_el[VV][pr * mesh->_NoLevels + Level + 1] - mesh->_off_el[VV][pr * mesh->_NoLevels + Level];
          }

          for( int iel = 0;
               iel <    mesh->_off_el[VV][off_proc + Level + 1]
               - mesh->_off_el[VV][off_proc + Level]; iel++ ) {
            int elem_lev = iel + sum_elems_prev_sd_at_lev;
            int dof_pos_lev = dofmap->GetDof( Level, KK, ivar, elem_lev );
            for( uint is = 0; is < NRE[mesh->_eltype_flag[VV]]; is++ ) {
              sol_on_cells[cel * NRE[mesh->_eltype_flag[VV]] + is] = ( * eqn->_LinSolver[Level]->_EPSC )( dof_pos_lev ) * eqn->_refvalue[ ivar + dofmap->_VarOff[KK] ];
            }
            cel++;
          }
        }

        std::ostringstream varname;
        varname << eqn->_var_names[ ivar + dofmap->_VarOff[KK] ] << "_" << grname.str();         //   std::string varname = grname.str() + "/" + _var_names[_nvars[QQ]+_nvars[LL]+ivar];
        hsize_t dimsf[2];
        dimsf[0] = NGeomObjOnWhichToPrint[KK];
        dimsf[1] = 1;
        XDMFWriter::print_Dhdf5( file_id, varname.str(), dimsf, sol_on_cells );

      } //end KK

      delete [] sol_on_cells;

//         H5Gclose(group_id);

    } //end Level

    H5Fclose( file_id );  //TODO VALGRIND

    for( int fe = 0; fe < QL; fe++ )  {
      delete fe_in[fe];
    }


    return;
  }


// ===================================================
/// This function reads the system solution from namefile.h5
//TODO this must be modified in order to take into account KK element dofs
  void XDMFWriter::read_system_solutions( const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn ) {
//this is done in parallel

    std::cout << "read_system_solutions still has to be written for CONSTANT elements, BEWARE!!! ==============================  " << std::endl;

    const uint Level = mesh->_NoLevels - 1;

    const uint offset   =       mesh->_NoNodesXLev[mesh->_NoLevels - 1];

    // file to read
    double* sol = new double[offset]; // temporary vector
    hid_t  file_id = H5Fopen( namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

    // reading loop over system varables
    for( uint ivar = 0; ivar < dofmap->_nvars[LL] + dofmap->_nvars[QQ]; ivar++ ) {
      uint el_nds = mesh->_elnodes[VV][QQ];
      if( ivar >= dofmap->_nvars[QQ] ) el_nds = mesh->_elnodes[VV][LL];
      // reading ivar param
      std::ostringstream grname;
      grname << eqn->_var_names[ivar] << "_" << "LEVEL" << Level;
      XDMFWriter::read_Dhdf5( file_id, grname.str(), sol );
      double Irefval = 1. / eqn->_refvalue[ivar]; // units

      // storing  ivar variables (in parallell)
      for( int iel = 0; iel <  mesh->_off_el[VV][mesh->_iproc * mesh->_NoLevels + mesh->_NoLevels]
           - mesh->_off_el[VV][mesh->_iproc * mesh->_NoLevels + mesh->_NoLevels - 1]; iel++ ) {
        uint elem_gidx = ( iel + mesh->_off_el[VV][mesh->_iproc * mesh->_NoLevels + mesh->_NoLevels - 1] ) * NVE[ mesh->_geomelem_flag[mesh->get_dim() - 1] ][BIQUADR_FE];
        for( uint i = 0; i < el_nds; i++ ) {  // linear and quad
          int k = mesh->_el_map[VV][elem_gidx + i]; // the global node
          eqn->_LinSolver[mesh->_NoLevels - 1]->_EPS->set( dofmap->GetDof( mesh->_NoLevels - 1, QQ, ivar, k ), sol[k]*Irefval );   // set the field
        }
      }
    }

    eqn->_LinSolver[mesh->_NoLevels - 1]->_EPS->localize( * eqn->_LinSolver[mesh->_NoLevels - 1]->_EPSC );
    // clean
    H5Fclose( file_id );
    delete []sol;

    return;
  }



// =====================================================================
/// This function  defines the boundary conditions for  DA systems:
// how does this function behave when e.g. both velocity and pressure are linear?
//this function runs over ALL the domain, so it's not only on a portion.
//ok we must think of this function in terms of the quantities
//we must also keep in mind that the DofMap is built by following a GIVEN ORDER;
// so, if you have 2 quadratic and 1 linear, you choose to put first the quadratic then the linear
// here, we loop over BOUNDARY ELEMENTS
//for every boundary element, loop over its nodes
// pick the coordinate of that node
//for every node, pick all its degrees of freedom from the dofmap
//at THAT COORDINATE, the DoF loop corresponds to a Quantity loop



// ================================================================
/// This function prints the boundary conditions: for quad and linear fem
//ok, this routine prints FLAGS, not VALUES.
//therefore, I'm not interested in INTERPOLATING FLAGS.
//Flags are either one or zero, no intermediate values are requested
//if all the nodes are 0, then it means do the pressure integral ---> put a zero
//if at least one node is 1 ---> put a 1. (anyway, the bc_p are not used in "volume" manner)
//the problem is always that we are printing LINEAR VARIABLES on a QUADRATIC mesh,
// so actually the interpolated values are there just for printing...
//

//TODO I must do in such a way as to print bc just like x_old,
// it is supposed to be exactly the same routine

//TODO I should print these bc flags, as well as the solution vectors, for ALL LEVELS!
// in this way I could see what happens at every level
// I already have the connectivities at all levels, thanks to gencase
// The only problem is that in 3D I do not see the connectivities, so I should convert all the levels to LINEARIZED REFINED

//also, apart from very little things, the routine is very similar to printing any field defined on each level
// And it's an integer and when you do the average you need to use the ceil() function...

//ok, when you print a field on a coarser grid, you need to specify a coarser topology,
// but also a coarser number of points...
// but, when I print only BOUNDARY MESH or VOLUME MESH at all levels,
// I put the FINE numbers and it's ok...
// so how does it pick the coordinates?


// A = NUMBER SEEN in PARAVIEW
// B = NUMBER in HDF5

//Ok, when i load coarser meshes, the number of cells changes but not the number of points...
// but the thing is: it seems that the NODES are ORDERED in SUCH A WAY THAT
// you have first COARSER then

//AAA: the equation A=B only holds with the WHOLE MESH!!!
// When you add some FILTER in paraview, for instance if you do a CLIP,
// Then the numbers DO NOT CORRESPOND ANYMORE!!!

//Ok, if i have to print a data attribute on a grid, i need to provide a SMALLER NUMBER of COORDINATES,
// or maybe make a LARGER DATA VECTOR...
//Due to the way we order the nodes, I guess we can say in the XDMF that the array CAN BE CUT! Let's see...
//You cannot put in the XDMF a SMALLER NUMBER than the dimension of the HDF5 ARRAY,
//otherwise it gives SEGMENTATION FAULT!
//So I need to have HDF5 fields for the coordinates of EACH LEVEL!
//I dont wanna put useless values, so I'll just print the coordinates for each level...

//Now there is a problem with the TIME printing... I am including all the solution files but it is picking up
// the COARSE LEVEL instead of the FINEST ONE

//TODO the problem is that now every solution file contains DIFFERENT GRIDS,
// and when I load the TIME file then it loads THE FIRST ONE IT MEETS in the SOL FILES!!!

// AAA, ecco l'errore! la connettivita' che percorriamo per interpolare i valori lineari
// e' quella FINE, ma noi ora dobbiamo prendere quella DI CIASCUN LIVELLO SEPARATAMENTE!


  void XDMFWriter::write_bc( const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn, const int* bc, int** bc_fe_kk ) {

    std::vector<GeomElemBase*> fe_in( QL );
    for( int fe = 0; fe < QL; fe++ )    fe_in[fe] = GeomElemBase::build( mesh->_geomelem_id[mesh->get_dim() - 1 - VV].c_str(), fe );

    hid_t file_id = H5Fopen( namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

    std::string  bdry_suffix = DEFAULT_BDRY_SUFFIX;

    const uint Lev_pick_bc_NODE_dof = mesh->_NoLevels - 1; //we use the FINE Level as reference

    // ==========================================
    // =========== FOR ALL LEVELS ===============
    // ==========================================
    for( uint Level = 0; Level < mesh->_NoLevels; Level++ )  {

      std::ostringstream grname;
      grname << "LEVEL" << Level;

      int NGeomObjOnWhichToPrint[QL];
      NGeomObjOnWhichToPrint[QQ] = mesh->_NoNodesXLev[Level];
      NGeomObjOnWhichToPrint[LL] = mesh->_NoNodesXLev[Level];
      NGeomObjOnWhichToPrint[KK] = mesh->_n_elements_vb_lev[VV][Level] * NRE[mesh->_eltype_flag[VV]];

      const uint n_nodes_lev = mesh->_NoNodesXLev[Level];
      int* sol_on_Qnodes = new int[n_nodes_lev];  //this vector will contain the values of ONE variable on ALL the QUADRATIC nodes
      //for the quadratic variables it'll be just a copy, for the linear also interpolation

      // ===================================
      // ========= QUADRATIC ================
      // ===================================
      for( uint ivar = 0; ivar < eqn->_dofmap._nvars[QQ]; ivar++ )        {

        for( uint isubdom = 0; isubdom < mesh->_NoSubdom; isubdom++ ) {
          uint off_proc = isubdom * mesh->_NoLevels;

          for( int fine_node = mesh->_off_nd[QQ][off_proc];
               fine_node < mesh->_off_nd[QQ][off_proc + Level + 1]; fine_node ++ ) {

            int pos_in_sol_vec_lev = eqn->_dofmap.GetDof( Lev_pick_bc_NODE_dof, QQ, ivar, fine_node );
            int pos_on_Qnodes_lev  = mesh->_Qnode_fine_Qnode_lev[Level][ fine_node ];

            sol_on_Qnodes[ pos_on_Qnodes_lev ] = bc[pos_in_sol_vec_lev];
          }
        }  //end subd

        std::ostringstream var_name;
        var_name << eqn->_var_names[ ivar + eqn->_dofmap._VarOff[QQ] ] << "_" << grname.str() << bdry_suffix;
        hsize_t dimsf[2];
        dimsf[0] = NGeomObjOnWhichToPrint[QQ];
        dimsf[1] = 1;
        XDMFWriter::print_Ihdf5( file_id, var_name.str(), dimsf, sol_on_Qnodes );
      }

      // ===================================
      // ========= LINEAR ==================
      // ===================================
      uint elnds[QL_NODES];
      elnds[QQ] = mesh->_elnodes[VV][QQ];
      elnds[LL] = mesh->_elnodes[VV][LL];
      double* elsol_c = new double[elnds[LL]];

      for( uint ivar = 0; ivar < eqn->_dofmap._nvars[LL]; ivar++ )   {

        for( uint isubdom = 0; isubdom < mesh->_NoSubdom; isubdom++ ) {
          uint off_proc = isubdom * mesh->_NoLevels;

          for( int fine_node =   mesh->_off_nd[QQ][off_proc];
               fine_node <   mesh->_off_nd[QQ][off_proc]
               + mesh->_off_nd[LL][off_proc + Level + 1]
               - mesh->_off_nd[LL][off_proc]; fine_node++ ) {

            int pos_in_sol_vec_lev = eqn->_dofmap.GetDof( Lev_pick_bc_NODE_dof, LL, ivar, fine_node );
            int pos_on_Qnodes_lev  = mesh->_Qnode_fine_Qnode_lev[Level][ fine_node ];

            sol_on_Qnodes[ pos_on_Qnodes_lev ] =  bc[pos_in_sol_vec_lev];
          }
        }

        //  2bB element interpolation over the fine mesh
        for( uint iproc = 0; iproc < mesh->_NoSubdom; iproc++ ) {
          uint off_proc = iproc * mesh->_NoLevels;
          int iel_b = mesh->_off_el[VV][off_proc + Level];
          int iel_e = mesh->_off_el[VV][off_proc + Level + 1];
          for( int iel = 0; iel < ( iel_e - iel_b ); iel++ ) {

            for( uint in = 0; in < elnds[LL]; in++ )  {
              int pos_Qnode_fine = mesh->_el_map[VV][( iel + iel_b ) * elnds[QQ] + in ];
              int pos_Qnode_lev  = mesh->_Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];
              elsol_c[in] = sol_on_Qnodes[ pos_Qnode_lev ];
            }

            for( uint in = 0; in < elnds[QQ]; in++ ) {  // mid-points
              double sum = 0;
              for( uint jn = 0; jn < elnds[LL]; jn++ ) {
                sum += fe_in[LL]->get_prol( in * elnds[LL] + jn ) * elsol_c[jn];
              }

              int pos_Qnode_fine = mesh->_el_map[VV][( iel + iel_b ) * elnds[QQ] + in ];
              int pos_Qnode_lev  = mesh->_Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];

              sol_on_Qnodes[ pos_Qnode_lev ] = ceil( sum );  //the ceiling, because you're putting double over int!
            }
          }
        } // 2bB end interpolation over the fine mesh

        std::ostringstream var_name;
        var_name << eqn->_var_names[ ivar + eqn->_dofmap._VarOff[LL] ] << "_" << grname.str() << bdry_suffix;
        hsize_t  dimsf[2];
        dimsf[0] = NGeomObjOnWhichToPrint[LL];
        dimsf[1] = 1;
        XDMFWriter::print_Ihdf5( file_id, var_name.str(), dimsf, sol_on_Qnodes );
      } // ivar

      delete [] elsol_c;
      delete [] sol_on_Qnodes;

      // ===================================
      // ========= CONSTANT ================
      // ===================================
      for( uint ivar = 0; ivar < eqn->_dofmap._nvars[KK]; ivar++ )        {

        int* sol_on_cells;
        sol_on_cells = new int[ NGeomObjOnWhichToPrint[KK] ];

        int cel = 0;
        for( uint iproc = 0; iproc < mesh->_NoSubdom; iproc++ ) {
          uint off_proc = iproc * mesh->_NoLevels;

          int sum_elems_prev_sd_at_lev = 0;
          for( uint pr = 0; pr < iproc; pr++ ) {
            sum_elems_prev_sd_at_lev += mesh->_off_el[VV][pr * mesh->_NoLevels + Level + 1] - mesh->_off_el[VV][pr * mesh->_NoLevels + Level];
          }

          for( int iel = 0;
               iel <    mesh->_off_el[VV][off_proc + Level + 1]
               - mesh->_off_el[VV][off_proc + Level]; iel++ ) {
            for( uint is = 0; is < NRE[mesh->_eltype_flag[VV]]; is++ ) {
              sol_on_cells[cel * NRE[mesh->_eltype_flag[VV]] + is] = bc_fe_kk[Level][iel + sum_elems_prev_sd_at_lev + ivar * mesh->_n_elements_vb_lev[VV][Level]]; //this depends on level!
            }
            cel++;
          }
        }

        std::ostringstream var_name;
        var_name << eqn->_var_names[ ivar + eqn->_dofmap._VarOff[KK] ] << "_" << grname.str() << bdry_suffix;
        hsize_t dimsf[2];
        dimsf[0] = NGeomObjOnWhichToPrint[KK];
        dimsf[1] = 1;
        XDMFWriter::print_Ihdf5( file_id, var_name.str(), dimsf, sol_on_cells );

        delete [] sol_on_cells;

      } //end KK

    } //end Level


    H5Fclose( file_id );

    for( int fe = 0; fe < QL; fe++ )  {
      delete fe_in[fe];
    }

    return;
  }



// ========================================================
/// It prints the volume/boundary Mesh (connectivity) in Xdmf format
//this is where the file mesh.xmf goes: it is a good rule that the
//file .xmf and the .h5 it needs should be in the same folder (so that READER and DATA are ALWAYS TOGETHER)
//well, there can be MANY READERS for ONE DATA file (e.g. sol.xmf, case.xmf)
// or MANY DATA files for  ONE READER,
//but if you put ALL THE DATA in THE SAME FOLDER as the READER(s)
//then they are always independent and the data will always be readable

  void XDMFWriter::PrintMeshXDMF( const std::string output_path, const MultiLevelMeshTwo& mesh, const uint order_fe ) {

    if( mesh._iproc == 0 ) {

      std::string meshname[VB];
      meshname[VV] = "Volume";
      meshname[BB] = "Boundary";

      std::string ord_string;
      if( order_fe == BIQUADR_FE )     {
        ord_string = "_biquadratic";
      }
      else if( order_fe == LINEAR_FE ) {
        ord_string = "_linear";
      }

      std::string     ext_xdmf = DEFAULT_EXT_XDMF;
      std::string       ext_h5 = DEFAULT_EXT_H5;
      std::string     basemesh = DEFAULT_BASEMESH;
      std::string      connlin = DEFAULT_BASEMESH_LIN;

      std::ostringstream top_file;
      top_file  << "mesh" << ord_string << ext_h5;
      std::ostringstream geom_file;
      geom_file << "mesh" << DEFAULT_BASEMESH_BIQ << ext_h5;

      std::ostringstream namefile;
      namefile << output_path << "/" << "mesh" << ord_string << ext_xdmf;

      std::ofstream out( namefile.str().c_str() );

      if( out.fail() ) {
        std::cout << " The file is not open" << std::endl;
        abort();
      }

      out << "<?xml version=\"1.0\" ?> \n";
      out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n";
      out << " [ <!ENTITY HeavyData \"mesh.h5 \"> ] ";
      out << ">\n";
      out << " \n";
      out << "<Xdmf> \n" << "<Domain> \n";


      for( uint vb = 0; vb < VB; vb++ )  {

        uint n_children;
        if( order_fe == BIQUADR_FE )     {
          n_children = 1;
        }
        else if( order_fe == LINEAR_FE ) {
          n_children = NRE[mesh._eltype_flag[vb]];
        }


        for( uint l = 0; l < mesh._NoLevels; l++ ) {

          out << "<Grid Name=\"" << meshname[vb].c_str() << "_L" << l << "\"> \n";

          PrintXDMFTopGeom( out, top_file, geom_file, l, vb, mesh, order_fe );

          std::ostringstream  pid_field;
          pid_field << "/PID/PID_VB" << vb << "_L" << l;
          PrintXDMFAttribute( out, top_file.str(), pid_field.str(), "PID", "Scalar", "Cell", "Int", mesh._n_elements_vb_lev[vb][l]*n_children, 1 );

          out << "</Grid> \n";

        }
      }

      out << "</Domain> \n";
      out << "</Xdmf> \n";
      out.close();


    } //end iproc


    return;
  }



// ========================================================
//print topology and geometry, useful for both case.xmf and sol.xmf
  void XDMFWriter::PrintXDMFTopGeom( std::ofstream& out,
                                     std::ostringstream& top_file,
                                     std::ostringstream& geom_file, const uint Level, const uint vb, const MultiLevelMeshTwo& mesh, const uint order_fe ) {

#ifdef HAVE_HDF5

    uint n_children;
    if( order_fe == BIQUADR_FE )     {
      n_children = 1;
    }
    else if( order_fe == LINEAR_FE ) {
      n_children = NRE[mesh._eltype_flag[vb]];
    }
    else {
      std::cout << "Mesh Not supported" << std::endl;
      abort();
    }

    uint nel = mesh._n_elements_vb_lev[vb][Level];

    std::ostringstream hdf_field;
    hdf_field << _elems_name << "/VB" << vb << "/" << _conn  << "_L" << Level;

    std::ostringstream coord_lev;
    coord_lev << "_L" << Level;

    PrintXDMFTopology( out, top_file.str(), hdf_field.str(), type_el[order_fe][mesh._eltype_flag[vb]], nel * n_children, nel * n_children, NVE[mesh._eltype_flag[vb]][order_fe] );

    PrintXDMFGeometry( out, geom_file.str(), _nodes_name + "/COORD/X", coord_lev.str(), "X_Y_Z", "Double", mesh._NoNodesXLev[Level], 1 );

#endif

    return;
  }



  void XDMFWriter::PrintSubdomFlagOnCellsAllVBAllLevHDF5( hid_t& file, std::string filename, const MultiLevelMeshTwo& mesh, const uint order ) {

    if( mesh._iproc == 0 )   {

      hid_t group_id = H5Gcreate( file, "/PID", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

      for( int  vb = 0; vb < VB; vb++ ) {
        for( int  Level = 0; Level < mesh._NoLevels; Level++ )  {
          PrintSubdomFlagOnCellsHDF5( vb, Level, filename, mesh, order );
        }
      }

      H5Gclose( group_id );

    }  //end iproc

    return;
  }


// ===============================================================
/// this function is done only by _iproc == 0!
/// it prints the PID index on the cells of the linear mesh
  void XDMFWriter::PrintSubdomFlagOnCellsHDF5( const uint vb, const int l, std::string filename, const MultiLevelMeshTwo& mesh, const uint order ) {

    if( mesh._iproc == 0 )   {

      uint n_children;
      if( order == BIQUADR_FE )      n_children = 1;
      else if( order == LINEAR_FE ) n_children = NRE[mesh._eltype_flag[vb]];
      else {
        std::cout << "Mesh Not supported" << std::endl;
        abort();
      }

      uint n_elements =  mesh._n_elements_vb_lev[vb][l];
      int* ucoord;
      ucoord = new int[n_elements * n_children];
      int cel = 0;
      for( uint iproc = 0; iproc <  mesh._NoSubdom; iproc++ ) {
        for( int iel =   mesh._off_el[vb][ l   + iproc * mesh._NoLevels];
             iel <  mesh._off_el[vb][ l + 1 + iproc * mesh._NoLevels]; iel++ ) {
          for( uint is = 0; is < n_children; is++ )
            ucoord[cel * n_children + is] = iproc;
          cel++;
        }
      }

      hid_t file_id = H5Fopen( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
      hsize_t dimsf[2];
      dimsf[0] = n_elements * n_children;
      dimsf[1] = 1;
      std::ostringstream pidname;
      pidname << "/PID/PID_VB" << vb << "_L" << l;

      XDMFWriter::print_Ihdf5( file_id, pidname.str(), dimsf, ucoord );
      H5Fclose( file_id );


    }

    return;
  }



// ========================================================
/// It manages the printing in Xdmf format
  void XDMFWriter::PrintMeshLinear( const std::string output_path, const MultiLevelMeshTwo& mesh ) {

    const uint iproc = mesh._iproc;
    if( iproc == 0 ) {
      PrintConnAllLEVAllVBLinearHDF5( output_path, mesh );
      PrintMeshXDMF( output_path, mesh, LINEAR_FE );
    }

    return;
  }



// ========================================================
/// It prints the connectivity in hdf5 format
  void XDMFWriter::PrintConnAllLEVAllVBLinearHDF5( const std::string output_path, const MultiLevelMeshTwo& mesh ) {

    std::string auxvb[VB];
    auxvb[0] = "0";
    auxvb[1] = "1";

    std::string    basemesh = DEFAULT_BASEMESH;
    std::string    ext_h5   = DEFAULT_EXT_H5;
    std::string    connlin  = DEFAULT_BASEMESH_LIN;

    std::ostringstream namefile;
    namefile << output_path << "/" <<  connlin << ext_h5;

    std::cout << namefile.str() << std::endl;
    hid_t file = H5Fcreate( namefile.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

//================================
    std::string elems_fem = _elems_name;
    hid_t subgroup_id_0 = H5Gcreate( file, elems_fem.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    for( uint vb = 0; vb < VB; vb++ )  {

      std::string elems_fem_vb = _elems_name + "/VB" + auxvb[vb];

      hid_t subgroup_id = H5Gcreate( file, elems_fem_vb.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

      for( uint l = 0; l < mesh._NoLevels; l++ )   PrintConnLinearHDF5( file, l, vb, mesh );

      H5Gclose( subgroup_id );

    }

    H5Gclose( subgroup_id_0 );

    PrintSubdomFlagOnCellsAllVBAllLevHDF5( file, namefile.str(), mesh, LINEAR_FE );

//================================

    H5Fclose( file );
    return;
  }



  void XDMFWriter::PrintConnLinearHDF5( hid_t file, const uint Level, const uint vb, const MultiLevelMeshTwo& mesh ) {

    int conn[8][8];  //TODO this is the largest dimension, bad programming
    uint* gl_conn;

    uint icount = 0;
    uint mode = NVE[ mesh._geomelem_flag[mesh._dim - 1 - vb] ][BIQUADR_FE];
    uint n_elements = mesh._n_elements_vb_lev[vb][Level];
    uint nsubel, nnodes;

    switch( mesh._dim )   {
      case 2:  {
        switch( mode ) {
            // -----------------------------------
          case 9:	// Quad 9  0-4-1-5-2-6-3-7-8
            gl_conn = new uint[n_elements * 4 * 4];
            conn[0][0] = 0;
            conn[0][1] = 4;
            conn[0][2] = 8;
            conn[0][3] = 7;// quad4  0-4-8-7
            conn[1][0] = 4;
            conn[1][1] = 1;
            conn[1][2] = 5;
            conn[1][3] = 8;// quad4  4-1-5-8
            conn[2][0] = 8;
            conn[2][1] = 5;
            conn[2][2] = 2;
            conn[2][3] = 6;// quad4  8-5-2-6
            conn[3][0] = 7;
            conn[3][1] = 8;
            conn[3][2] = 6;
            conn[3][3] = 3;// quad4  7-8-6-3
            nsubel = 4;
            nnodes = 4;
            break;
//=====================
          case 6:	// Quad 9  0-4-1-5-2-6-3-7-8
            gl_conn = new uint[n_elements * 4 * 3];
            conn[0][0] = 0;
            conn[0][1] = 3;
            conn[0][2] = 5; // quad4  0-4-8-7
            conn[1][0] = 3;
            conn[1][1] = 4;
            conn[1][2] = 5; // quad4  4-1-5-8
            conn[2][0] = 3;
            conn[2][1] = 1;
            conn[2][2] = 4; // quad4  8-5-2-6
            conn[3][0] = 4;
            conn[3][1] = 2;
            conn[3][2] = 5; // quad4  7-8-6-3
            nsubel = 4;
            nnodes = 3;
            break;
//===================
          case 3: // boundary edge 2 linear 0-2-1
            gl_conn = new uint[n_elements * 2 * 2];
            conn[0][0] = 0;
            conn[0][1] = 2;		// element 0-2
            conn[1][0] = 2;
            conn[1][1] = 1;		// element 1-2
            nsubel = 2;
            nnodes = 2;
            break;
            // -----------------------------------------
          default:   // interior 3D
            gl_conn = new uint[n_elements * mode];
            for( uint n = 0; n < mode; n++ ) conn[0][n] = n;
            nsubel = 1;
            nnodes = mode;
            break;
        } //end switch mode 1

        break; //end case 2
      }
      case 3:  {
        switch( mode ) {
            // ----------------------
          case  27: //  Hex 27 (8 Hex8)
            gl_conn = new uint[n_elements * 8 * 8];
            conn[0][0] = 0;
            conn[0][1] = 8;
            conn[0][2] = 20;
            conn[0][3] = 11;
            conn[0][4] = 12;
            conn[0][5] = 21;
            conn[0][6] = 26;
            conn[0][7] = 24;
            conn[1][0] = 8;
            conn[1][1] = 1;
            conn[1][2] = 9;
            conn[1][3] = 20;
            conn[1][4] = 21;
            conn[1][5] = 13;
            conn[1][6] = 22;
            conn[1][7] = 26;
            conn[2][0] = 11;
            conn[2][1] = 20;
            conn[2][2] = 10;
            conn[2][3] = 3;
            conn[2][4] = 24;
            conn[2][5] = 26;
            conn[2][6] = 23;
            conn[2][7] = 15;
            conn[3][0] = 20;
            conn[3][1] = 9;
            conn[3][2] = 2;
            conn[3][3] = 10;
            conn[3][4] = 26;
            conn[3][5] = 22;
            conn[3][6] = 14;
            conn[3][7] = 23;
            conn[4][0] = 12;
            conn[4][1] = 21;
            conn[4][2] = 26;
            conn[4][3] = 24;
            conn[4][4] = 4;
            conn[4][5] = 16;
            conn[4][6] = 25;
            conn[4][7] = 19;
            conn[5][0] = 21;
            conn[5][1] = 13;
            conn[5][2] = 22;
            conn[5][3] = 26;
            conn[5][4] = 16;
            conn[5][5] = 5;
            conn[5][6] = 17;
            conn[5][7] = 25;
            conn[6][0] = 24;
            conn[6][1] = 26;
            conn[6][2] = 23;
            conn[6][3] = 15;
            conn[6][4] = 19;
            conn[6][5] = 25;
            conn[6][6] = 18;
            conn[6][7] = 7;
            conn[7][0] = 26;
            conn[7][1] = 22;
            conn[7][2] = 14;
            conn[7][3] = 23;
            conn[7][4] = 25;
            conn[7][5] = 17;
            conn[7][6] = 6;
            conn[7][7] = 18;
            nsubel = 8;
            nnodes = 8;
            break;
            // ---------------------------------------
          case  10: // Tet10
            gl_conn = new uint[n_elements * 8 * 4];
            conn[0][0] = 0;
            conn[0][1] = 4;
            conn[0][2] = 6;
            conn[0][3] = 7;
            conn[1][0] = 4;
            conn[1][1] = 1;
            conn[1][2] = 5;
            conn[1][3] = 8;
            conn[2][0] = 5;
            conn[2][1] = 2;
            conn[2][2] = 6;
            conn[2][3] = 9;
            conn[3][0] = 7;
            conn[3][1] = 8;
            conn[3][2] = 9;
            conn[3][3] = 3;
            conn[4][0] = 4;
            conn[4][1] = 8;
            conn[4][2] = 6;
            conn[4][3] = 7;
            conn[5][0] = 4;
            conn[5][1] = 5;
            conn[5][2] = 6;
            conn[5][3] = 8;
            conn[6][0] = 5;
            conn[6][1] = 9;
            conn[6][2] = 6;
            conn[6][3] = 8;
            conn[7][0] = 7;
            conn[7][1] = 6;
            conn[7][2] = 9;
            conn[7][3] = 8;
            nsubel = 8;
            nnodes = 4;
            break;
            // ---------------------------------------
          case 6:	// Tri6  0-3-1-4-2-5
            gl_conn = new uint[n_elements * 4 * 3];
            conn[0][0] = 0;
            conn[0][1] = 3;
            conn[0][2] = 5; // quad4  0-4-8-7
            conn[1][0] = 3;
            conn[1][1] = 4;
            conn[1][2] = 5; // quad4  4-1-5-8
            conn[2][0] = 3;
            conn[2][1] = 1;
            conn[2][2] = 4; // quad4  8-5-2-6
            conn[3][0] = 4;
            conn[3][1] = 2;
            conn[3][2] = 5; // quad4  7-8-6-3
            nsubel = 4;
            nnodes = 3;
            break;
            // ---------------------------------------
          case 9:  // Quad9 elements ( 4 Quad4)
            gl_conn = new uint[n_elements * 4 * 4];
            conn[0][0] = 0;
            conn[0][1] = 4;
            conn[0][2] = 8;
            conn[0][3] = 7;
            conn[1][0] = 4;
            conn[1][1] = 1;
            conn[1][2] = 5;
            conn[1][3] = 8;
            conn[2][0] = 8;
            conn[2][1] = 5;
            conn[2][2] = 2;
            conn[2][3] = 6;
            conn[3][0] = 7;
            conn[3][1] = 8;
            conn[3][2] = 6;
            conn[3][3] = 3;
            nsubel = 4;
            nnodes = 4;
            break;
            // -----------------------------------------
          default:   // interior 3D
            gl_conn = new uint[n_elements * mode];
            for( uint n = 0; n < mode; n++ ) conn[0][n] = n;
            nsubel = 1;
            nnodes = mode;
            break;
        } //end switch mode 2

        break;
      } //end case 3

      default:
        std::cout << "PrintConnLin not working" << std::endl;
        abort();
        break;
    }

    // mapping
    for( uint iproc = 0; iproc < mesh._NoSubdom; iproc++ ) {
      for( int el = mesh._off_el[vb][iproc * mesh._NoLevels + Level];
           el < mesh._off_el[vb][iproc * mesh._NoLevels + Level + 1]; el++ ) {
        for( uint se = 0; se < nsubel; se++ ) {
          for( uint i = 0; i < nnodes; i++ ) {
            const uint pos = el * mode + conn[se][i];
            uint Qnode_fine = mesh._el_map[vb][pos];
            uint Qnode_lev = mesh._Qnode_fine_Qnode_lev[Level][Qnode_fine];
            gl_conn[icount] = Qnode_lev;
            icount++;
          }
        }
      }
    }

    // Print mesh in hdf files
    hsize_t dimsf[2];
    dimsf[0] = icount;
    dimsf[1] = 1;
    hid_t dtsp = H5Screate_simple( 2, dimsf, NULL );
    std::ostringstream Name;
    Name << _elems_name << "/VB" << vb << "/" << _conn << "_L" << Level;

    hid_t dtset = H5Dcreate( file, Name.str().c_str(), H5T_NATIVE_INT, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );  //TODO VALGRIND
    H5Dwrite( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gl_conn );

    H5Dclose( dtset );

    H5Sclose( dtsp );

    delete[] gl_conn;

    return;
  }



// ==========================================================
//prints conn and stuff for either vol or bdry mesh
//When you have to construct the connectivity,
//you go back to the libmesh elem ordering,
//then you pick the nodes of that element in LIBMESH numbering,
//then you pick the nodes in femus NUMBERING,
//and that's it
  void XDMFWriter::PrintElemVBBiquadraticHDF5( hid_t file,
      const uint vb ,
      const std::vector<int>& nd_libm_fm,
      ElemStoBase** el_sto_in,
      const std::vector<std::pair<int, int> >  el_fm_libm_in, const MultiLevelMeshTwo& mesh ) {

// const unsigned from_libmesh_to_xdmf[27] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};  //id
// const unsigned from_libmesh_to_xdmf[27] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,21,22,23,24,20,25,26};  //from libmesh to femus
// const unsigned from_libmesh_to_xdmf[27]    = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,24,22,21,23,20,25,26};  //from libmesh to xdmf

    std::ostringstream name;

    std::string auxvb[VB];
    auxvb[0] = "0";
    auxvb[1] = "1";
    std::string elems_fem = _elems_name;
    std::string elems_fem_vb = elems_fem + "/VB" + auxvb[vb];  //VV later

    hid_t subgroup_id = H5Gcreate( file, elems_fem_vb.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    hsize_t dimsf[2];
    // NoElements ------------------------------------
    dimsf[0] = mesh._NoLevels;
    dimsf[1] = 1;
    XDMFWriter::print_UIhdf5( file, ( elems_fem_vb + "/NExLEV" ), dimsf, mesh._n_elements_vb_lev[vb] );
    // offset
    dimsf[0] = mesh._NoSubdom * mesh._NoLevels + 1;
    dimsf[1] = 1;
    XDMFWriter::print_Ihdf5( file, ( elems_fem_vb + "/OFF_EL" ), dimsf, mesh._off_el[vb] );

    //here you pick all the elements at all levels,
    //and you print their connectivities according to the libmesh ordering
    int* tempconn;
    tempconn = new int[mesh._n_elements_sum_levs[vb]*mesh._elnodes[vb][QQ]]; //connectivity of all levels

// // //     if (_elnodes[vb][QQ] == 27)  { //HEX27
// // //
// // //        for (int ielem=0;ielem<_n_elements_sum_levs[vb];ielem++) {
// // //         for (uint inode=0;inode<_elnodes[vb][QQ];inode++) {
// // //             int el_libm =   el_fm_libm_in[ielem].second;
// // //             int nd_libm = el_sto_in[el_libm]->_elnds[inode];
// // //             tempconn[ from_libmesh_to_xdmf[inode] + ielem*_elnodes[vb][QQ] ] = nd_libm_fm[nd_libm];
// // //         }
// // //     }
// // //
// // //     }//HEX27
// // //     else {

    for( int ielem = 0; ielem < mesh._n_elements_sum_levs[vb]; ielem++ ) {
      for( uint inode = 0; inode < mesh._elnodes[vb][QQ]; inode++ ) {
        int el_libm =   el_fm_libm_in[ielem].second;
        int nd_libm = el_sto_in[el_libm]->_elnds[inode];
        tempconn[inode + ielem * mesh._elnodes[vb][QQ]] = nd_libm_fm[nd_libm];
      }
    }

// // //   }

    dimsf[0] = mesh._n_elements_sum_levs[vb] * mesh._elnodes[vb][QQ];
    dimsf[1] = 1;
    XDMFWriter::print_Ihdf5( file, ( elems_fem_vb + "/CONN" ), dimsf, tempconn );

    // level connectivity ---------------------------------
    for( int ilev = 0; ilev < mesh._NoLevels; ilev++ ) {

      int* conn_lev = new int[mesh._n_elements_vb_lev[vb][ilev]*mesh._elnodes[vb][QQ]]; //connectivity of ilev


      int ltot = 0;
      for( int iproc = 0; iproc < mesh._NoSubdom; iproc++ ) {
        for( int iel = mesh._off_el[vb][iproc * mesh._NoLevels + ilev];
             iel < mesh._off_el[vb][iproc * mesh._NoLevels + ilev + 1]; iel++ ) {
          for( uint inode = 0; inode < mesh._elnodes[vb][QQ]; inode++ ) {
            conn_lev[ltot * mesh._elnodes[vb][QQ] + inode ] =
              tempconn[  iel * mesh._elnodes[vb][QQ] + inode ];
          }
          ltot++;
        }
      }


      dimsf[0] = mesh._n_elements_vb_lev[vb][ilev] * mesh._elnodes[vb][QQ];
      dimsf[1] = 1;

      name.str( "" );
      name << elems_fem_vb << "/CONN" << "_L" << ilev ;
      XDMFWriter::print_Ihdf5( file, name.str(), dimsf, conn_lev );
      //clean
      delete []conn_lev;

    }


    delete []tempconn;

    H5Gclose( subgroup_id );

    return;


  }



// ========================================================
/// Read mesh from hdf5 file (namefile)
//is this function for ALL PROCESSORS
// or only for PROC==0? Seems to be for all processors
// TODO do we need the leading "/" for opening a dataset?
// This routine reads the mesh file and also makes it NONDIMENSIONAL, so that everything is solved on a nondimensional mesh
  void XDMFWriter::ReadMeshAndNondimensionalizeBiquadraticHDF5( const std::string output_path, MultiLevelMeshTwo& mesh )   {

    std::string    basemesh = DEFAULT_BASEMESH;
    std::string      ext_h5 = DEFAULT_EXT_H5;

    std::ostringstream meshname;
    meshname << output_path << "/" << basemesh  << ext_h5;

//==================================
// OPEN FILE
//==================================
    std::cout << " Reading mesh from= " <<  meshname.str() <<  std::endl;
    hid_t  file_id = H5Fopen( meshname.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
//   if (file_id < 0) {std::cout << "MultiLevelMeshTwo::read_c(): File Mesh input in data_in is missing"; abort();}
// he does not do this check, if things are wrong the H5Fopen function detects the error

//==================================
// DLS (Dimension, Levels, Subdomains)
// =====================
    uint topdata[3];
    XDMFWriter::read_UIhdf5( file_id, "/DIM_LEVS_PROCS", topdata );

//==================================
// CHECKS
// =====================
    if( mesh._NoLevels !=  topdata[1] )  {
      std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the number of mesh levels is " <<
                "different in the mesh file and in the configuration file" << std::endl;
      abort();
    }


    if( mesh._NoSubdom != topdata[2] )  {
      std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the number of mesh subdomains is " << mesh._NoSubdom
                << " while the processor size of this run is " << paral::get_size()
                << ". Re-run gencase and your application with the same number of processors" << std::endl;
      abort();
    }

//alright, there is a check to make: if you run gencase in 3D, and then you run main in 2D then things go wrong...
//this is because we read the dimension from 'dimension', we should read it from the mesh file in principle,
//in fact it is that file that sets the space in which we are simulating...
//I'll put a check

    if( mesh._dim != topdata[0] ) {
      std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the mesh dimension is " << mesh._dim
                << " while the dimension in the configuration file is " << mesh.get_dim()
                << ". Recompile either gencase or your application appropriately" << std::endl;
      abort();
    }
//it seems like it doesn't print to file if I don't put the endline "<< std::endl".
//Also, "\n" seems to have no effect, "<< std::endl" must be used
//This fact doesn't seem to be related to PARALLEL processes that abort sooner than the others

//==================================
// FEM element DoF number
// =====================
//Reading this is not very useful... well, it may be a check
    XDMFWriter::read_UIhdf5( file_id, "/ELNODES_VB", &mesh._type_FEM[0] );

    for( int vb = 0; vb < VB; vb++ ) {
      if( mesh._type_FEM[vb] !=  NVE[ mesh._geomelem_flag[mesh._dim - 1 - vb] ][BIQUADR_FE] )  {
        std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the element type of the mesh is" <<
                  "different from the element type as given by the GeomEl" << std::endl;
        abort();
      }
    }

// ===========================================
// ===========================================
//  NODES
// ===========================================
// ===========================================

// ++++++++++++++++++++++++++++++++++++++++++++++++++
// nodes X lev
// ++++++++++++++++++++++++++++++++++++++++++++++++++
    mesh._NoNodesXLev = new uint[mesh._NoLevels + 1];
    XDMFWriter::read_UIhdf5( file_id, "/NODES/MAP/NDxLEV", mesh._NoNodesXLev );

// ===========================================
//  COORDINATES  (COORD)
// ===========================================
    //in the mesh file now I have the coordinates for all levels, but let me read only those of the FINEST level
    double ILref = 1. / mesh.get_Lref();
    int lev_for_coords = mesh._NoLevels - 1;
    uint n_nodes = mesh._NoNodesXLev[lev_for_coords];
    mesh._xyz = new double[mesh._dim * n_nodes];
    double* coord;
    coord = new double[n_nodes];
    for( uint kc = 0; kc < mesh._dim; kc++ ) {
      std::ostringstream Name;
      Name << "NODES/COORD/X" << kc + 1 << "_L" << lev_for_coords;
      XDMFWriter::read_Dhdf5( file_id, Name.str().c_str(), coord );
      for( uint inode = 0; inode < n_nodes; inode++ ) mesh._xyz[inode + kc * n_nodes] = coord[inode] * ILref;  //NONdimensionalization!!!
    }
    delete []coord;

// ===================================================
//  OFF_ND
// ===================================================
    mesh._off_nd = new int*[QL_NODES];

    for( int fe = 0; fe < QL_NODES; fe++ )    {
      mesh._off_nd[fe] = new int[mesh._NoSubdom * mesh._NoLevels + 1];
      std::ostringstream namefe;
      namefe << "/NODES/MAP/OFF_ND" << "_F" << fe;
      XDMFWriter::read_Ihdf5( file_id, namefe.str().c_str(), mesh._off_nd[fe] );
    }


    // ====================================================
// NODE MAP
// ====================================================
    // node mapping
    //this map "mixes" linear and quadratic, in the sense that
    // if you have levels 0 1 2, you can define
    // 3 levels for the QUADRATIC DOFS, which are A=0 B=1 C=2,
    //and similarly you can define 3 levels for the linear dofs.
    //In order to be consistent with the levels you can define
    // another auxiliary level 3 which is
    // the coarse linear level
    //in each _node_map there is no "subdomain" stuff, only LEVEL
    // Each node_map contains the list of NODES of THAT LEVEL,
    //of course in the FINE NODE NUMBERING
    //in practice gives us, for every level, the FINE NODE INDICES corresponding to the DOFS
    //so it is like an INVERSE DOF MAP: from DOF TO NODE.
    //if you use levels 0 1 2, you do from QUADRATIC dof to FINE NODE
    //if you use levels 3 0 1, you do from LINEAR dof to FINE NODE

    uint n_nodes_top = mesh._NoNodesXLev[mesh._NoLevels - 1];
    mesh._Qnode_lev_Qnode_fine = new uint *[mesh._NoLevels + 1];
    mesh._Qnode_fine_Qnode_lev = new int *[mesh._NoLevels + 1];

    for( uint ilev = 0; ilev <= mesh._NoLevels; ilev++ ) {  //loop on Extended levels
      mesh._Qnode_lev_Qnode_fine[ilev] = new uint [mesh._NoNodesXLev[ilev]];
      mesh._Qnode_fine_Qnode_lev[ilev] = new int [n_nodes_top];  //THIS HAS TO BE INT because it has -1!!!
      std::ostringstream Name;
      Name << "/NODES/MAP/MAP" << "_XL" << ilev;
      XDMFWriter::read_Ihdf5( file_id, Name.str().c_str(), mesh._Qnode_fine_Qnode_lev[ilev] );
      for( uint inode = 0; inode < n_nodes_top; inode++ ) {
        int val_lev = mesh._Qnode_fine_Qnode_lev[ilev][inode];
        if( val_lev != -1 ) mesh._Qnode_lev_Qnode_fine[ilev][ val_lev ] = inode;  //this doesnt have -1 numbers
      }
      //This is how you read the _node_map:
      //- you remove the "-1" that come from the gencase
      //- and you dont put the content (we dont need it!) but the position,
      //   i.e. the node index in the fine numbering
    }

// ===========================================
//   /ELEMS
// ===========================================

// ===========================================
//   NUMBER EL
// ===========================================
    mesh._n_elements_vb_lev = new uint*[VB];
    for( uint vb = 0; vb < VB; vb++ ) {
      mesh._n_elements_vb_lev[vb] = new uint[mesh._NoLevels];
      std::ostringstream Name;
      Name << "/ELEMS/VB" << vb  << "/NExLEV";
      XDMFWriter::read_UIhdf5( file_id, Name.str().c_str(), mesh._n_elements_vb_lev[vb] );
    }

// ===========================================
//   OFF_EL
// ===========================================
    mesh._off_el = new int*[VB];

    for( int vb = 0; vb < VB; vb++ )    {
      mesh._off_el[vb] = new int [mesh._NoSubdom * mesh._NoLevels + 1];
      std::ostringstream offname;
      offname << "/ELEMS/VB" << vb << "/OFF_EL";
      XDMFWriter::read_Ihdf5( file_id, offname.str().c_str(), mesh._off_el[vb] );
    }

// ===========================================
//   CONNECTIVITY
// ===========================================
    mesh._el_map = new uint*[VB];
    for( int vb = 0; vb < VB; vb++ )    {
      mesh._el_map[vb] = new uint [mesh._off_el[vb][mesh._NoSubdom * mesh._NoLevels]*NVE[ mesh._geomelem_flag[mesh._dim - 1 - vb] ][BIQUADR_FE]];
      std::ostringstream elName;
      elName << "/ELEMS/VB" << vb  << "/CONN";
      XDMFWriter::read_UIhdf5( file_id, elName.str().c_str(), mesh._el_map[vb] );
    }

// ===========================================
//   BDRY EL TO VOL EL
// ===========================================
    mesh._el_bdry_to_vol = new int*[mesh._NoLevels];
    for( uint lev = 0; lev < mesh._NoLevels; lev++ )    {
      mesh._el_bdry_to_vol[lev] = new int[mesh._n_elements_vb_lev[BB][lev]];
      std::ostringstream btov;
      btov << "/ELEMS/BDRY_TO_VOL_L" << lev;
      XDMFWriter::read_Ihdf5( file_id, btov.str().c_str(), mesh._el_bdry_to_vol[lev] );
    }

// ===========================================
//  CLOSE FILE
// ===========================================
    H5Fclose( file_id );

    return;

  }


// ===============================================================
  void XDMFWriter::PrintMeshBiquadraticHDF5( const std::string output_path, const MultiLevelMeshTwo& mesh )  {

    std::ostringstream name;

    std::string basemesh  = DEFAULT_BASEMESH;
    std::string ext_h5    = DEFAULT_EXT_H5;

    std::ostringstream inmesh;
    inmesh << output_path << "/" << basemesh << ext_h5;

//==================================
// OPEN FILE
//==================================
    hid_t file = H5Fcreate( inmesh.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

//==================================
// DLS (Dimension, Levels, Subdomains)
// =====================
    int* tdata;

    tdata    = new int[3];
    tdata[0] = mesh._dim;
    tdata[1] = mesh._NoLevels;
    tdata[2] = mesh._NoSubdom;

    hsize_t dimsf[2];
    dimsf[0] = 3;
    dimsf[1] = 1;
    XDMFWriter::print_Ihdf5( file, "DIM_LEVS_PROCS", dimsf, tdata );
    delete [] tdata;

//==================================
// FEM element DoF number
// =====================
    int* ttype_FEM;
    ttype_FEM = new int[VB];

    for( uint vb = 0; vb < VB; vb++ )   ttype_FEM[vb] = NVE[ mesh._geomelem_flag[mesh.get_dim() - 1 - vb] ][BIQUADR_FE];

    dimsf[0] = VB;
    dimsf[1] = 1;
    XDMFWriter::print_Ihdf5( file, "ELNODES_VB", dimsf, ttype_FEM );

// ===========================================
// ===========================================
//  NODES
// ===========================================
// ===========================================
    hid_t group_id = H5Gcreate( file, XDMFWriter::_nodes_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

// ++++ NODES/MAP ++++++++++++++++++++++++++++++++++++++++++++++
    std::string ndmap = XDMFWriter::_nodes_name + "/MAP";
    hid_t subgroup_id = H5Gcreate( file, ndmap.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// nodes X lev
    std::string ndxlev =  ndmap + "/NDxLEV";
    dimsf[0] = mesh._NoLevels + 1;
    dimsf[1] = 1;
    XDMFWriter::print_UIhdf5( file, ndxlev.c_str(), dimsf, mesh._NoNodesXLev );

// ++++++++++++++++++++++++++++++++++++++++++++++++++
// node map (XL, extended levels)
// ++++++++++++++++++++++++++++++++++++++++++++++++++
    dimsf[0] = mesh._n_nodes;
    dimsf[1] = 1;

    for( int ilev = 0; ilev < mesh._NoLevels + 1; ilev++ ) {
      name.str( "" );
      name << ndmap << "/MAP" << "_XL" << ilev;
      XDMFWriter::print_Ihdf5( file, name.str(), dimsf, mesh._Qnode_fine_Qnode_lev[ilev] );
    }

// ++++++++++++++++++++++++++++++++++++++++++++++++++
//   OFF_ND: node offset quadratic and linear
// ++++++++++++++++++++++++++++++++++++++++++++++++++
    dimsf[0] = mesh._NoSubdom * mesh._NoLevels + 1;
    dimsf[1] = 1;
    for( int fe = 0; fe < QL_NODES; fe++ ) {
      std::ostringstream namefe;
      namefe <<  ndmap << "/OFF_ND" << "_F" << fe;
      XDMFWriter::print_Ihdf5( file, namefe.str(), dimsf, mesh._off_nd[fe] );
    }

    H5Gclose( subgroup_id );

// ===========================================
//  COORDINATES  (COORD)
// ===========================================
    //ok, we need to print the coordinates of the nodes for each LEVEL
    //The array _nod_coords holds the coordinates of the FINE Qnodes
    //how do we take the Qnodes of each level?
    //Well, I'd say we need to use  _off_nd for the quadratics
    //Ok, the map _nd_fm_libm goes from FINE FEMUS NODE ORDERING to FINE LIBMESH NODE ORDERING
    //I need to go from the QNODE NUMBER at LEVEL
    //to the QNODE NUMBER at FINE LEVEL according to FEMUS,
    //and finally to the number at fine level according to LIBMESH.
    //The fine level is the only one where the Qnode numbering is CONTIGUOUS



    std::string ndcoords = XDMFWriter::_nodes_name + "/COORD";
    subgroup_id = H5Gcreate( file, ndcoords.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // nodes are PRINTED ACCORDING to FEMUS ordering, which is inode, i.e. the INVERSE of v[inode].second
    //you use this because you do DIRECTLY a NODE LOOP
    //now let us loop coordinates over ALL LEVELS
    for( int l = 0; l < mesh._NoLevels; l++ )  {
      double* xcoord = new double[mesh._NoNodesXLev[l]];
      for( int kc = 0; kc < 3; kc++ ) {

        int Qnode_lev = 0;
        for( uint isubdom = 0; isubdom < mesh._NoSubdom; isubdom++ ) {
          uint off_proc = isubdom * mesh._NoLevels;
          for( int k1 = mesh._off_nd[QQ][off_proc];
               k1 < mesh._off_nd[QQ][off_proc + l + 1 ]; k1++ ) {
            int Qnode_fine_fm = mesh._Qnode_lev_Qnode_fine[l][Qnode_lev];
            xcoord[Qnode_lev] = mesh._nd_coords_libm[ mesh._nd_fm_libm[Qnode_fine_fm].second + kc * mesh._n_nodes ];
            Qnode_lev++;
          }
        } //end subdomain

// old        for (int inode=0; inode <_n_nodes_lev[l];inode++)  xcoord[inode] = _nod_coords[_nd_fm_libm[inode].second+kc*_n_nodes]; //the offset is fine
        dimsf[0] = mesh._NoNodesXLev[l];
        dimsf[1] = 1;
        name.str( "" );
        name << ndcoords << "/X" << kc + 1 << "_L" << l;
        XDMFWriter::print_Dhdf5( file, name.str(), dimsf, xcoord );
      }

      delete [] xcoord;
    }  //levels

    H5Gclose( subgroup_id );

    H5Gclose( group_id );

// ===========================================
// ===========================================
//   /ELEMS
// ===========================================
// ===========================================

    group_id = H5Gcreate( file, XDMFWriter::_elems_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    ElemStoBase** elsto_out;   //TODO delete it!
    elsto_out = new ElemStoBase*[mesh._n_elements_sum_levs[VV]];
    for( int i = 0; i < mesh._n_elements_sum_levs[VV]; i++ ) {
      elsto_out[i] = static_cast<ElemStoBase*>( mesh._el_sto[i] );
    }

    ElemStoBase** elstob_out;
    elstob_out = new ElemStoBase*[mesh._n_elements_sum_levs[BB]];
    for( int i = 0; i < mesh._n_elements_sum_levs[BB]; i++ ) {
      elstob_out[i] = static_cast<ElemStoBase*>( mesh._el_sto_b[i] );
    }

    XDMFWriter::PrintElemVBBiquadraticHDF5( file, VV, mesh._nd_libm_fm, elsto_out, mesh._el_fm_libm, mesh );
    XDMFWriter::PrintElemVBBiquadraticHDF5( file, BB, mesh._nd_libm_fm, elstob_out, mesh._el_fm_libm_b, mesh );

    // ===============
    // print child to father map for all levels for BOUNDARY ELEMENTS
    // ===============

    for( int lev = 0; lev < mesh._NoLevels; lev++ )  {
      std::ostringstream   bname;
      bname << XDMFWriter::_elems_name << "/BDRY_TO_VOL_L" << lev;
      dimsf[0] = mesh._n_elements_vb_lev[BB][lev];
      dimsf[1] = 1;
      XDMFWriter::print_Ihdf5( file, bname.str(), dimsf, mesh._el_child_to_fath[lev] );
    }




//             std::cout <<  "==================" << std::endl;
//          for (int i=0;i<_n_elements_vb_lev[BB][lev];i++)  {
//              std::cout <<  _el_child_to_fath[lev][i] << std::endl;
//        }
// 	 }



// ok, so, i had to do two cast vectors so i could pass them both to the print_elem_vb routine
//now i have to be careful in destroying these, in relation with their fathers...

    //delete temp  //I AM HERE TODO
//     for (int i=0;i< _n_elements_sum_levs_vb[BB];i++) {   delete /*[]*/ elstob_out[i]; } // delete [] el_sto[i]; with [] it doesnt work
//    delete [] elstob_out;
//     for (int i=0;i< _n_elements_sum_levs_vb[VV];i++) {   delete /*[]*/ elsto_out[i]; } // delete [] el_sto[i]; with [] it doesnt work
//    delete [] elsto_out;
    H5Gclose( group_id );

// ===========================================
//   PID
// ===========================================
    PrintSubdomFlagOnCellsAllVBAllLevHDF5( file, inmesh.str().c_str(), mesh, BIQUADR_FE );

// ===========================================
//  CLOSE FILE
// ===========================================
    H5Fclose( file );

    //============
    delete [] ttype_FEM;

    return;
  }



// =======================================================================
// For every matrix, compute dimension, position, len, offlen
// this function concerns a COUPLE of finite element families, that's it
// fe_row = rows
// fe_col = columns
// the arrays that I print with HDF5 should already be filled at this point and have their dimension explicit
// the dimension of POS (Mat)        is "count"
// the dimension of LEN (len)        is "n_dofs_lev_fe[fe_row]+1"
// the dimension of OFFLEN (len_off) is "n_dofs_lev_fe[fe_row]+1"

  void XDMFWriter::PrintOneVarMatrixHDF5( const std::string& name, const std::string& groupname, uint** n_dofs_lev_fe, int count,
                                          int* Mat, int* len, int* len_off,
                                          int fe_row, int fe_col, int* FELevel ) {


    hid_t file = H5Fopen( name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

    hsize_t dimsf[2];
    dimsf[1] = 1;  //for all the cases

    std::ostringstream fe_couple;
    fe_couple <<  "_F" << fe_row << "_F" << fe_col;

    //==== DIM ========
    std::ostringstream name0;
    name0  << groupname << "/" << "DIM" << fe_couple.str();
    dimsf[0] = 2;
    int rowcln[2];
    rowcln[0] = n_dofs_lev_fe[fe_row][FELevel[fe_row]]; //row dimension
    rowcln[1] = n_dofs_lev_fe[fe_col][FELevel[fe_col]]; //column dimension
    XDMFWriter::print_Ihdf5( file, name0.str().c_str(), dimsf, rowcln );

    //===== POS =======
    std::ostringstream name1;
    name1 << groupname << "/"  << "POS" << fe_couple.str();
    dimsf[0] = count;
    XDMFWriter::print_Ihdf5( file, name1.str().c_str(), dimsf, Mat );

    //===== LEN =======
    std::ostringstream name2;
    name2 << groupname << "/"  << "LEN" << fe_couple.str();
    dimsf[0] = rowcln[0] + 1;
    XDMFWriter::print_Ihdf5( file, name2.str().c_str(), dimsf, len );

    //==== OFFLEN ========
    std::ostringstream name3;
    name3 << groupname << "/"  << "OFFLEN" << fe_couple.str();
    dimsf[0] = rowcln[0] + 1;
    XDMFWriter::print_Ihdf5( file, name3.str().c_str(), dimsf, len_off );

    H5Fclose( file );  //TODO this file seems to be closed TWICE, here and in the calling function
    //you open and close the file for every (FE_ROW,FE_COL) couple

    return;
  }

// ===============================================================
// here pay attention, the EXTENDED LEVELS are not used for distinguishing
//
  void XDMFWriter::PrintOneVarMGOperatorHDF5( const std::string& filename, const std::string& groupname, uint* n_dofs_lev, int count, int* Op_pos, double* Op_val, int* len, int* len_off, int FELevel_row, int FELevel_col, int fe ) {

    hid_t file = H5Fopen( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );  //TODO questo apri interno e' per assicurarsi che il file sia aperto... quello fuori credo che non serva...
    // e invece credo che quello serva per CREARE il file, altrimenti non esiste

    hsize_t dimsf[2];
    dimsf[1] = 1;  //for all the cases

    std::ostringstream fe_family;
    fe_family <<  "_F" << fe;
    //==== DIM ========
    std::ostringstream name0;
    name0 << groupname << "/" << "DIM" << fe_family.str();
    dimsf[0] = 2;
    int rowcln[2]; //for restrictor row is coarse, column is fine
    rowcln[0] = n_dofs_lev[FELevel_row];
    rowcln[1] = n_dofs_lev[FELevel_col];

    XDMFWriter::print_Ihdf5( file, name0.str().c_str(), dimsf, rowcln );
    //===== POS =======
    std::ostringstream name1;
    name1 << groupname << "/" << "POS" << fe_family.str();
    dimsf[0] = count;
    XDMFWriter::print_Ihdf5( file, name1.str().c_str(), dimsf, Op_pos );
    //===== VAL =======
    std::ostringstream name1b;
    name1b << groupname << "/" << "VAL" << fe_family.str();
    dimsf[0] = count;
    XDMFWriter::print_Dhdf5( file, name1b.str().c_str(), dimsf, Op_val );
    //===== LEN =======
    std::ostringstream name2;
    name2 << groupname << "/" << "LEN" << fe_family.str();
    dimsf[0] = rowcln[0] + 1;
    XDMFWriter::print_Ihdf5( file, name2.str().c_str(), dimsf, len );
    //==== OFFLEN ========
    std::ostringstream name3;
    name3 << groupname << "/" << "OFFLEN" << fe_family.str();
    dimsf[0] = rowcln[0] + 1;
    XDMFWriter::print_Ihdf5( file, name3.str().c_str(), dimsf, len_off );


    //TODO Attenzione!!! Per scrivere una matrice HDF5 la vuole TUTTA INSIEME!!!
    //QUINDI DEVI FARE IN MODO DA ALLOCARLA SU UNO SPAZIO DI MEMORIA CONTIGUO!!!
    // per questo devi fare l'allocazione non con i new separati per ogni riga,
    //perche' in tal modo ogni riga puo' essere allocata in uno spazio di memoria
    //differente, rompendo la contiguita'!

    //TODO I must do two template functions out of these

    //ok let me print also in matrix form, more easily readable
    //i have a one dimensional array, i want to break it into two dimensional arrays
    //ok if I just create it it puts only zeros
    //now i want to extract vectors and in each row put the vector i want
    //Ok, non lo so , lo faccio alla brutta, traduco il mio vettore lungo
    // in una matrice e stampo la matrice
    int n_cols = 40; //facciamo una roba in grande ora
    int** mat_op = new int*[rowcln[0]];
    mat_op[0] = new int[ rowcln[0]*n_cols ];

    int  sum_prev_rows = 0;
    for( uint i = 0; i < rowcln[0]; i++ ) {
      mat_op[i] = mat_op[0] + i * n_cols;  //sum of pointers, to keep contiguous in memory!!!
      for( uint j = 0; j < n_cols; j++ ) {
        mat_op[i][j]  = 0;
        if( j < ( len[i + 1] - len[i] ) ) {
          mat_op[i][j] = Op_pos[ sum_prev_rows + j ];
        }
      }
      sum_prev_rows += ( len[i + 1] - len[i] );
    }

    dimsf[0] = rowcln[0];
    dimsf[1] = n_cols;
    std::ostringstream name4;
    name4 << groupname << "/" << "POS_MAT" << fe_family.str();
    XDMFWriter::print_Ihdf5( file, name4.str().c_str(), dimsf, &mat_op[0][0] );

    delete [] mat_op[0];
    delete [] mat_op;

//===========================
    n_cols = 40;
    double** mat = new double*[rowcln[0]];
    mat[0] = new double[ rowcln[0]*n_cols ];

    sum_prev_rows = 0;
    for( uint i = 0; i < rowcln[0]; i++ ) {
      mat[i] = mat[0] + i * n_cols;  //sum of pointers, to keep contiguous in memory!!!
      for( uint j = 0; j < n_cols; j++ ) {
        mat[i][j]  = 0;
        if( j < ( len[i + 1] - len[i] ) ) {
          mat[i][j] = Op_val[ sum_prev_rows + j ];
        }
      }
      sum_prev_rows += ( len[i + 1] - len[i] );
    }

    dimsf[0] = rowcln[0];
    dimsf[1] = n_cols;
    std::ostringstream name5;
    name5 << groupname << "/" << "VAL_MAT" << fe_family.str();
    XDMFWriter::print_Dhdf5( file, name5.str().c_str(), dimsf, &mat[0][0] );

    delete [] mat[0];
    delete [] mat;


    H5Fclose( file );

    return;
  }


// MultiLevelProblem
// =================================================================
/// This function prints the attributes into the corresponding hdf5 file
  void XDMFWriter::PrintSolHDF5Linear( const std::string output_path, const uint t_flag, const MultiLevelProblem& ml_prob ) {

    const uint    iproc  = ml_prob.GetMeshTwo()._iproc;
    if( iproc == 0 ) {

      const uint     ndigits  = DEFAULT_NDIGITS;
      std::string    basesol  = DEFAULT_BASESOL;
      std::string     ext_h5  = DEFAULT_EXT_H5;
      std::ostringstream filename;
      filename << output_path << "/" << basesol << "." << std::setw( ndigits ) << std::setfill( '0' ) << t_flag << ext_h5;

      hid_t   file = H5Fcreate( filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
      H5Fclose( file );

      MultiLevelProblem::const_system_iterator pos = ml_prob.begin();
      MultiLevelProblem::const_system_iterator pos_e = ml_prob.end();
      for( ; pos != pos_e; pos++ )    {
        SystemTwo* eqn = static_cast<SystemTwo*>( pos->second );
        XDMFWriter::write( filename.str(), & ml_prob.GetMeshTwo(), & ( eqn->_dofmap ), eqn );
      }

    } //end print iproc

    return;
  }

// ===================================================================
/// It prints the attributes in  Xdmf format for one time step
// Qui dobbiamo separare le NODE variables dalle CELL variables!
//le node variables sono le _nvars[QQ] e _nvars[LL], le Cell sono _nvars[KK]
//The SCALAR variables are ordered FIRST QUADRATIC, then LINEAR, then CONSTANT

//Ok, there is a problem when I want to print on coarser grids,
// it seems that it needs the coarse coordinates
//Ok i think i solved it

//Now the thing is: i need to print from FINEST LEVEL, so that i can see in paraview the FINEST sequence!!!
//Ok, now I want to print an XDMF for every level, because I've seen that paraview has quite a little bit
// of problems in reading files with multiple grids, so I'll put everything in separate files per level.
//The file name depends on the level so i will put it inside the level loop

// If you want to put all the grids in the SAME XDMF file, then it is better to put the FINEST first, and then all the others...
// This is due to the Paraview TIME Files that would read the FIRST GRID THEY MEET in the XDMF FILE!
//So you should do  in reverse order:  for (int l=NoLevels-1; l >= 0; l--)
// I DO NOT WANT TO SEPARATE the HDF5 file, because it works fine as a single file with no problems

  void XDMFWriter::PrintSolXDMFLinear( const std::string output_path, const uint t_step, const double curr_time, const MultiLevelProblem& ml_prob ) {

    const uint    iproc = ml_prob.GetMeshTwo()._iproc;
    if( iproc == 0 ) {

      const uint ndigits  = DEFAULT_NDIGITS;
      const uint NoLevels = ml_prob.GetMeshTwo()._NoLevels;

      //FE print
      std::string DofType[QL];
      DofType[QQ] = "Node";
      DofType[LL] = "Node";
      DofType[KK] = "Cell";

      std::string basesol     = DEFAULT_BASESOL;
      std::string basemesh    = DEFAULT_BASEMESH;
      std::string aux_xdmf    = DEFAULT_AUX_XDMF;
      std::string connlin     = DEFAULT_BASEMESH_LIN;
      std::string     ext_h5  = DEFAULT_EXT_H5;
      std::string    ext_xdmf = DEFAULT_EXT_XDMF;

      std::ostringstream top_file;
      top_file << connlin << ext_h5;
      std::ostringstream geom_file;
      geom_file << basemesh << ext_h5;

// =================================
// ============= LEVELS ============
// =================================

      for( uint l = 0; l < NoLevels; l++ ) {

        std::ostringstream filename_xdmf;
        filename_xdmf << output_path << "/" << basesol << "." << std::setw( ndigits ) << std::setfill( '0' ) << t_step << "_l" << l << ext_xdmf;
        std::ostringstream hdf_file;
        hdf_file << basesol << "." << std::setw( ndigits ) << std::setfill( '0' ) << t_step << ext_h5;

        std::ofstream out( filename_xdmf.str().c_str() );
        if( out.fail() ) {
          std::cout << "MultiLevelProblemTwo::print_soln_xmf: cannot print " << filename_xdmf.str().c_str() << std::endl;
          abort();
        }

        //  ++++++++++++ Header ++++++++++++++
        out << "<?xml version=\"1.0\" ?> \n";
        out << "<!DOCTYPE Xdmf SYSTEM ";
        out <<  "\"" << output_path << "/" << aux_xdmf << "\" \n";
//   out << " [ <!ENTITY HeavyData \"\"> ] ";
        out << ">\n";
        out << "<Xdmf> \n" << "<Domain> \n";


        int NGeomObjOnWhichToPrint[QL];
        NGeomObjOnWhichToPrint[QQ] = ml_prob.GetMeshTwo()._NoNodesXLev[l];
        NGeomObjOnWhichToPrint[LL] = ml_prob.GetMeshTwo()._NoNodesXLev[l];
        NGeomObjOnWhichToPrint[KK] = ml_prob.GetMeshTwo()._n_elements_vb_lev[VV][l] * NRE[ml_prob.GetMeshTwo()._eltype_flag[VV]];

        out << "<Grid Name=\"Volume_L" << l << "\"> \n";

        out << "<Time Value =\"" << curr_time << "\" /> \n";

        PrintXDMFTopGeom( out, top_file, geom_file, l, VV, ml_prob.GetMeshTwo(), LINEAR_FE );

        MultiLevelProblem::const_system_iterator pos1   = ml_prob.begin();
        MultiLevelProblem::const_system_iterator pos1_e = ml_prob.end();
        for( ; pos1 != pos1_e; pos1++ )   {
          SystemTwo* mgsol = static_cast<SystemTwo*>( pos1->second );
          int OffVarNames[QL];
          OffVarNames[QQ] = 0;
          OffVarNames[LL] = mgsol->_dofmap._nvars[QQ];
          OffVarNames[KK] = mgsol->_dofmap._nvars[QQ] + mgsol->_dofmap._nvars[LL];
          for( int fe = 0; fe < QL; fe++ )  {
            for( uint ivar = 0; ivar < mgsol->_dofmap._nvars[fe]; ivar++ )   {
              std::ostringstream var_name;
              var_name << mgsol->_var_names[ OffVarNames[fe] + ivar] << "_LEVEL" << l;
              PrintXDMFAttribute( out, hdf_file.str(), var_name.str(), var_name.str(), "Scalar", DofType[fe], "Double", NGeomObjOnWhichToPrint[fe], 1 );
            }
          } // end fe
        }

        out << "</Grid>\n";

        // ============= END GRID ===========

        out << "</Domain> \n" << "</Xdmf> \n";
        out.close();

      }//end levels

    } //end iproc==0

    return;
  }

// =================================================================
/// This function prints xdmf and hdf5 file
  void XDMFWriter::PrintSolLinear( const std::string output_path, const uint t_step, const double curr_time, const MultiLevelProblem& ml_prob ) {

    XDMFWriter::PrintSolHDF5Linear( output_path, t_step, ml_prob );
    XDMFWriter::PrintSolXDMFLinear( output_path, t_step, curr_time, ml_prob );

    return;
  }


// ====================================================
/// This function prints initial and boundary data in xdmf+hdf5 format
// of course whenever you change the fields printed in the case h5 file
// then you need to change also the xdmf file
// we should do a routine that for a given field prints both the hdf5 dataset
// and the  xdmf tag... well it's not so automatic, because you need to know
// what is the grid on which to print, bla bla bla
  void XDMFWriter::PrintCaseLinear( const std::string output_path, const uint t_init, const MultiLevelProblem& ml_prob ) {


    XDMFWriter::PrintCaseHDF5Linear( output_path, t_init, ml_prob );
    XDMFWriter::PrintCaseXDMFLinear( output_path, t_init, ml_prob );

    return;
  }


// =============================================================================
/// This function prints initial and boundary data in hdf5 fromat
/// in the file case.h5
  void XDMFWriter::PrintCaseHDF5Linear( const std::string output_path, const uint t_init, const MultiLevelProblem& ml_prob ) {

    const uint    iproc = ml_prob.GetMeshTwo()._iproc;
    if( iproc == 0 ) {

      const uint ndigits      = DEFAULT_NDIGITS;
      std::string    basecase = DEFAULT_BASECASE;
      std::string   ext_xdmf  = DEFAULT_EXT_XDMF;
      std::string     ext_h5  = DEFAULT_EXT_H5;

      std::ostringstream filename;
      filename << output_path << "/" << basecase << "." << std::setw( ndigits ) << std::setfill( '0' ) << t_init << ext_h5;

      hid_t file = H5Fcreate( filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

      H5Fclose( file );

      MultiLevelProblem::const_system_iterator pos   = ml_prob.begin();
      MultiLevelProblem::const_system_iterator pos_e = ml_prob.end();
      for( ; pos != pos_e; pos++ ) {
        SystemTwo* eqn = static_cast<SystemTwo*>( pos->second );
        XDMFWriter::write( filename.str(), & ml_prob.GetMeshTwo(), & ( eqn->_dofmap ), eqn );  // initial solution
        XDMFWriter::write_bc( filename.str(), & ml_prob.GetMeshTwo(), & ( eqn->_dofmap ), eqn, eqn->_bcond._bc, NULL );  // boundary condition
      }

    } //end iproc

    return;
  }



// ====================================================================
/// It prints the Xdmf file to read the initial and boundary conditions
//now, this function does exactly the same as print sol;
//more over, it prints the PID, and the BOUNDARY CONDITIONS
//so it prints, INITIAL CONDITIONS, BOUNDARY CONDITIONS, PID
//let us split so that we can have a unique function

//clearly, you need to know WHERE to print this file.
//so you need the absolute paths
//but, inside the lines of this file, you dont need to put the absolute paths,
//because you already know you'll not separate .xmf and .h5

  void XDMFWriter::PrintCaseXDMFLinear( const std::string output_path, const uint t_init, const MultiLevelProblem& ml_prob ) {

    const uint    iproc = ml_prob.GetMeshTwo()._iproc;
    if( iproc == 0 ) {

      const uint NoLevels = ml_prob.GetMeshTwo()._NoLevels;
      const uint ndigits  = DEFAULT_NDIGITS;

      std::string     basecase = DEFAULT_BASECASE;
      std::string     basemesh = DEFAULT_BASEMESH;
      std::string       ext_h5 = DEFAULT_EXT_H5;
      std::string     ext_xdmf = DEFAULT_EXT_XDMF;
      std::string     aux_xdmf = DEFAULT_AUX_XDMF;
      std::string      connlin = DEFAULT_BASEMESH_LIN;
      std::string  bdry_suffix = DEFAULT_BDRY_SUFFIX;

      std::ostringstream top_file;
      top_file << connlin << ext_h5;
      std::ostringstream geom_file;
      geom_file << basemesh << ext_h5;

      //FE print
      std::string DofType[QL];
      DofType[QQ] = "Node";
      DofType[LL] = "Node";
      DofType[KK] = "Cell";

      std::string var_name[VB];
      std::string var_type[VB];

// =================================
// ============= LEVELS ============
// =================================

      for( uint l = 0; l < NoLevels; l++ ) {

        std::ostringstream filename_xdmf;
        filename_xdmf << output_path << "/"
                      << basecase  << "." << std::setw( ndigits ) << std::setfill( '0' ) << t_init << "_l" << l << ext_xdmf;
        std::ostringstream hdf_file;
        hdf_file <<  basecase << "." << std::setw( ndigits ) << std::setfill( '0' ) << t_init << ext_h5;

        std::ofstream out( filename_xdmf.str().c_str() );
        if( out.fail() ) {
          std::cout << "MultiLevelProblemTwo::print_case_xmf: cannot print " << filename_xdmf.str().c_str() << std::endl;
          abort();
        }

        // BEGIN XDMF =======
        out << "<?xml version=\"1.0\" ?> \n";
        out << "<!DOCTYPE Xdmf SYSTEM ";
        out <<  "\"" << output_path << "/" << aux_xdmf << "\" \n";
//    out << " [ <!ENTITY HeavyData \"\"> ] ";
        out << ">\n";
        out << "<Xdmf> \n" << "<Domain> \n";


        int NGeomObjOnWhichToPrint[QL];
        NGeomObjOnWhichToPrint[QQ] = ml_prob.GetMeshTwo()._NoNodesXLev[l];
        NGeomObjOnWhichToPrint[LL] = ml_prob.GetMeshTwo()._NoNodesXLev[l];
        NGeomObjOnWhichToPrint[KK] = ml_prob.GetMeshTwo()._n_elements_vb_lev[VV][l] * NRE[ml_prob.GetMeshTwo()._eltype_flag[VV]];

        out << "<Grid Name=\"Volume_L" << l << "\"> \n";

        // TOPOLOGY GEOMETRY ===========
        PrintXDMFTopGeom( out, top_file, geom_file, l, VV, ml_prob.GetMeshTwo(), LINEAR_FE );

        // ATTRIBUTES FOR EACH SYSTEM ===========
        MultiLevelProblem::const_system_iterator pos1 = ml_prob.begin();
        MultiLevelProblem::const_system_iterator pos1_e = ml_prob.end();
        for( ; pos1 != pos1_e; pos1++ )   {
          SystemTwo* mgsol = static_cast<SystemTwo*>( pos1->second );
          int OffVarNames[QL];
          OffVarNames[QQ] = 0;
          OffVarNames[LL] = mgsol->_dofmap._nvars[QQ];
          OffVarNames[KK] = mgsol->_dofmap._nvars[QQ] + mgsol->_dofmap._nvars[LL];
          for( int fe = 0; fe < QL; fe++ )  {
            for( uint ivar = 0; ivar < mgsol->_dofmap._nvars[fe]; ivar++ )     {
              std::ostringstream  varstream;
              varstream << mgsol->_var_names[OffVarNames[fe] + ivar] << "_LEVEL" << l;
              var_name[VV] = varstream.str();
              var_type[VV] = "Double";
              var_name[BB] = var_name[VV] + bdry_suffix;
              var_type[BB] = "Int";
              for( int vb = 0; vb < VB; vb++ ) {
                PrintXDMFAttribute( out, hdf_file.str(), var_name[vb], var_name[vb], "Scalar", DofType[fe], var_type[vb], NGeomObjOnWhichToPrint[fe], 1 );
              }
            }
          } //end fe

        }  //end eqn

        out << "</Grid>\n";

        out << "</Domain> \n" << "</Xdmf> \n";
        out.close();

      } //end levels

    } //if iproc==0

    return;
  }

// ========================================================================
/// This function read the solution form all the system (restart)
  void XDMFWriter::ReadSol( const std::string output_path, const uint t_step, double& time_out, const MultiLevelProblem& ml_prob ) {

    const uint ndigits      = DEFAULT_NDIGITS;
    std::string    basesol  = DEFAULT_BASESOL;
    std::string   ext_xdmf  = DEFAULT_EXT_XDMF;
    std::string     ext_h5  = DEFAULT_EXT_H5;
// ---------------------------------------------------
    // reading time from from sol.N.xmf file
    // ---------------------------------------------------
    // open file -----------------------------
    std::ostringstream namefile;
    namefile << output_path << "/"
             << basesol << "." << std::setw( ndigits ) << std::setfill( '0' ) << t_step << "_l" << ( ml_prob.GetMeshTwo()._NoLevels - 1 ) << ext_xdmf;   //TODO here we should avoid doing this process TWICE because we already do it in the TransientSetup calling function

#ifdef DEFAULT_PRINT_INFO // --------  info ------------------
    std::cout << "\n MultiLevelProblemTwo::read_soln: Reading time  from "
              << namefile.str().c_str();
#endif  // -------------------------------------------
    std::ifstream in ;
    in.open( namefile.str().c_str() );  //associate the file stream with the name of the file
    if( !in.is_open() ) {
      std::cout << " read_soln: restart .xmf file not found "  << std::endl;
      abort();
    }

    // reading time from xmf file --------------
    std::string buf = "";
    while( buf != "<Time" ) in >> buf;
    in >> buf >> buf;
    buf = buf.substr( 2, buf.size() - 3 );
//create an istringstream from a string
    std::istringstream buffer( buf );
    double restart_time;
    buffer >> restart_time;

    //pass  the time value to the calling function
    time_out = restart_time;

//add parameter to system dont need that now
//   _utils.set_par("restartime",restart_time);

    // ---------------------------------------------------
    // reading data from  sol.N.h5
    // ---------------------------------------------------
    // file name -----------------------------------------
    namefile.str( "" );  //empty string
    namefile << output_path << "/"
             << basesol << "." << std::setw( ndigits ) << std::setfill( '0' ) << t_step << ext_h5;
    //if i put the path of this file to be relative, will the read depend on where I launched the executable...
    // or where the executable is I think... no, the path is given by where the executable is LAUNCHED

#ifdef DEFAULT_PRINT_INFO  // --------------- info ---------------
    std::cout << "\n MultiLevelProblemTwo::read_soln: Reading from file "
              << namefile.str().c_str() << std::endl;
#endif // ---------------------------------------------
    // loop reading over the variables ---------------------
    for( MultiLevelProblem::const_system_iterator eqn = ml_prob.begin(); eqn != ml_prob.end(); eqn++ ) {
      SystemTwo* mgsol = static_cast<SystemTwo*>( eqn->second );
    } //  loop --------------------------------------------------------

    return;
  }

} //end namespace femus


