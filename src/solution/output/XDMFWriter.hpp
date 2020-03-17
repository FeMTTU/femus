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

#ifndef __femus_solution_XDMFWriter_hpp__
#define __femus_solution_XDMFWriter_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "MultiLevelProblem.hpp"

#ifdef HAVE_HDF5

namespace femus {

  class DofMap;
  class MultiLevelMeshTwo;
  class SystemTwo;



  class XDMFWriter : public Writer {

    public:

      /** Constructor. */
      XDMFWriter( MultiLevelSolution* ml_sol );

      /** Constructor. */
      XDMFWriter( MultiLevelMesh* ml_mesh );

      /** Destructor */
      virtual ~XDMFWriter();

      /** write output function */
      void Write( const std::string output_path, const char order[], const std::vector < std::string >& vars = std::vector < std::string > (), const unsigned time_step = 0 ) ;

      /** write a wrapper file for paraview to open all the files of an history together */
      void write_solution_wrapper( const std::string output_path, const char type[] ) const;

      //==================
      static void transient_print_xmf( const std::string output_path, const uint t_idx_in, const uint t_idx_final, const int print_step, const uint nolevels_in );

      static void write_bc( const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn, const int* bc, int** bc_fe_kk );
      static void write( const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn ); ///prints on a "Quadratic-Linearized" Mesh //TODO this should be PrintNumericVector of the equation //Writer//
      static void  read_system_solutions( const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn );                     ///read from a "Quadratic-Linearized" Mesh                                      //Writer/Reader//

      //hdf5 ------------------------------------
      static hid_t print_Dhdf5( hid_t file, const std::string& name, hsize_t* dimsf, double* data );
      static hid_t print_Ihdf5( hid_t file, const std::string& name, hsize_t* dimsf, int* data );
      static hid_t print_UIhdf5( hid_t file, const std::string& name, hsize_t* dimsf, uint* data );
      static hid_t read_Dhdf5( hid_t file, const std::string& name, double* data );
      static hid_t read_Ihdf5( hid_t file, const std::string& name, int* data );
      static hid_t read_UIhdf5( hid_t file, const std::string& name, uint* data );

      /** MESH PRINTING */
      static void PrintXDMFAttribute( std::ofstream& outstream,
                                      std::string hdf5_filename,
                                      std::string hdf5_field,
                                      std::string attr_name,
                                      std::string attr_type,
                                      std::string attr_center,
                                      std::string data_type,
                                      int data_dim_row,
                                      int data_dim_col
                                    );

      static void PrintXDMFTopology( std::ofstream& outfstream,
                                     std::string hdf5_file,
                                     std::string hdf5_field,
                                     std::string top_type,
                                     int top_dim,
                                     int datadim_n_elems,
                                     int datadim_el_nodes
                                   );

      static void PrintXDMFGeometry( std::ofstream& outfstream,
                                     std::string hdf5_file,
                                     std::string hdf5_field,
                                     std::string coord_lev,
                                     std::string geom_type,
                                     std::string data_type,
                                     int data_dim_one,
                                     int data_dim_two );

      static void PrintMeshXDMF( const std::string output_path, const MultiLevelMeshTwo& mesh, const uint order_fe );

      static void PrintXDMFTopGeom( std::ofstream& out, std::ostringstream& top_file,
                                    std::ostringstream& geom_file,
                                    const uint Level,
                                    const uint vb,
                                    const MultiLevelMeshTwo& mesh, const uint order_fe );

      static void PrintSubdomFlagOnCellsAllVBAllLevHDF5( hid_t& file, std::string filename, const MultiLevelMeshTwo& mesh, const uint order );

      static void PrintSubdomFlagOnCellsHDF5( const uint vb, const int Level, std::string filename, const MultiLevelMeshTwo& mesh, const uint order );

      static void PrintMeshLinear( const std::string output_path, const MultiLevelMeshTwo& mesh );

      static void PrintConnAllLEVAllVBLinearHDF5( const std::string output_path, const MultiLevelMeshTwo& mesh );

      static void PrintConnLinearHDF5( hid_t file, const uint Level, const uint vb, const MultiLevelMeshTwo& mesh );

      static void PrintElemVBBiquadraticHDF5( hid_t file, const uint vb, const std::vector<int>& nd_libm_fm, ElemStoBase** el_sto_in, const std::vector<std::pair<int, int> >  el_fm_libm_in, const MultiLevelMeshTwo& mesh );

      static void ReadMeshAndNondimensionalizeBiquadraticHDF5( const std::string output_path, MultiLevelMeshTwo& mesh );

      static void PrintMeshBiquadraticHDF5( const std::string output_path, const MultiLevelMeshTwo& mesh );

      /** Matrices */
      static void PrintOneVarMatrixHDF5( const std::string& name, const std::string& groupname, uint** n_nodes_all, int count, int* Mat, int* len, int* len_off, int type1, int type2, int* FELevel );

      static void PrintOneVarMGOperatorHDF5( const std::string& filename, const std::string& groupname, uint* n_dofs_lev, int count, int* Op_pos, double* Op_val, int* len, int* len_off, int FELevel_row, int FELevel_col, int fe );

      /** MultiLevelProblem */
      static void PrintSolXDMFLinear( const std::string output_path, const uint t_step, const double curr_time, const MultiLevelProblem& ml_prob );

      static void PrintSolHDF5Linear( const std::string output_path, const uint t_flag, const MultiLevelProblem& ml_prob );

      static void PrintSolLinear( const std::string output_path, const uint t_step, const double curr_time, const MultiLevelProblem& ml_prob );

      static void PrintCaseXDMFLinear( const std::string output_path, const uint t_init, const MultiLevelProblem& ml_prob );

      static void PrintCaseHDF5Linear( const std::string output_path, const uint t_init, const MultiLevelProblem& ml_prob );

      static void PrintCaseLinear( const std::string output_path, const uint t_init, const MultiLevelProblem& ml_prob ); ///< Print ic and bc

      static void ReadSol( const std::string output_path, const uint t_step, double& time_out, const MultiLevelProblem& ml_prob ); ///< Read solution //TODO must be updated, not implemented

      /** Set if to print or not to prind the debugging variables */
      void SetDebugOutput( bool value ) {
        _debugOutput = value;
      }

    private:

      bool _debugOutput;

      static const std::string type_el[3][N_GEOM_ELS];

      static const std::string _nodes_name;
      static const std::string _elems_name;
      static const std::string _conn;
      //     std::string _nd_coord_folder;
//     std::string _el_pid_name;
//     std::string _nd_map_FineToLev;



  };


} //end namespace femus

#endif

#endif
