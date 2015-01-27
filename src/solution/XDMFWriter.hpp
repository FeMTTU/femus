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

#ifndef __xdmfwriter_h_
#define __xdmfwriter_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer.hpp"
#include "MultiLevelMeshTwo.hpp"

namespace femus {

class DofMap;
class MultiLevelMeshTwo;
class SystemTwo;



class XDMFWriter : public Writer {

public:

    /** Constructor. */
    XDMFWriter(MultiLevelSolution& ml_sol);
    
    /** Destructor */
    virtual ~XDMFWriter();

    /** write output function */
    virtual void write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step = 0);

    /** write a wrapper file for paraview to open all the files of an history toghether */
    void write_solution_wrapper(const char type[]) const;

  //==================    
   static void write_system_solutions_bc(const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn, const int* bc, int** bc_fe_kk);      
   static void write_system_solutions(const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn);   ///prints on a "Quadratic-Linearized" Mesh //TODO this should be PrintNumericVector of the equation //Writer//
   static void  read_system_solutions(const std::string namefile, const MultiLevelMeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn);                       ///read from a "Quadratic-Linearized" Mesh                                      //Writer/Reader// 
    
  //hdf5 ------------------------------------
   static hid_t print_Dhdf5(hid_t file,const std::string & name, hsize_t* dimsf,double* data);
   static hid_t print_Ihdf5(hid_t file,const std::string & name, hsize_t* dimsf,int* data);
   static hid_t print_UIhdf5(hid_t file,const std::string & name, hsize_t* dimsf,uint* data);
   static hid_t read_Dhdf5(hid_t file,const std::string & name,double* data);
   static hid_t read_Ihdf5(hid_t file,const std::string & name,int* data);
   static hid_t read_UIhdf5(hid_t file,const std::string & name,uint* data);

  //XDMF  
   static void PrintXDMFAttribute(std::ofstream& outstream, 
				      std::string hdf5_filename, 
				      std::string hdf5_field,
				      std::string attr_name,
				      std::string attr_type,
				      std::string attr_center,
				      std::string data_type,
				      int data_dim_row,
				      int data_dim_col
 				    );
  
  static void PrintXDMFTopology(std::ofstream& outfstream,
				     std::string hdf5_file,
				     std::string hdf5_field,
				     std::string top_type,
				     int top_dim,
				     int datadim_n_elems,
				     int datadim_el_nodes
				    );  
  
  static void PrintXDMFGeometry(std::ofstream& outfstream,
				     std::string hdf5_file,
				     std::string hdf5_field,
				     std::string coord_lev,
				     std::string geom_type,
				     std::string data_type,
				     int data_dim_one,
				     int data_dim_two); 
  
  static void PrintMultimeshXdmfBiquadratic(const std::string output_path, const MultiLevelMeshTwo & mesh);
  
  static void PrintXDMFAllLEVAllVB(const std::string output_path, const MultiLevelMeshTwo & mesh);
  
  static void PrintXDMFGridVBLinear(std::ofstream& out, std::ostringstream& top_file,
			      std::ostringstream& geom_file,
			      const uint Level,
			      const uint vb,
			      const MultiLevelMeshTwo & mesh);
  
   /** */
  static void PrintXDMFTopologyGeometryLinear(std::ofstream& out,const unsigned Level, const unsigned vb, const MultiLevelMeshTwo& mesh);

  static void PrintSubdomFlagOnQuadrCells(const int vb, const int Level, std::string filename, const MultiLevelMeshTwo & mesh);
  
  static void PrintSubdomFlagOnLinCells(std::string filename, const MultiLevelMeshTwo & mesh);
  
private:
  
   static const std::string type_el[4][6];
    

};


} //end namespace femus



#endif
