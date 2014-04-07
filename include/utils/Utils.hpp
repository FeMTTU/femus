#ifndef __io_hpp__
#define __io_hpp__

#include <map>
#include <string>
#include "hdf5.h"

#include "Typedefs.hpp"


namespace IO {
  
  //hdf5 ------------------------------------
   hid_t print_Dhdf5(hid_t file,const std::string & name, hsize_t dimsf[],double data[]);
   hid_t print_Ihdf5(hid_t file,const std::string & name, hsize_t dimsf[],int data[]);
   hid_t read_Dhdf5(hid_t file,const std::string & name,double data[]);
   hid_t read_Ihdf5(hid_t file,const std::string & name,int data[]);

  //XDMF  
   void PrintXDMFAttribute(std::ofstream& outstream, 
				      std::string hdf5_filename, 
				      std::string hdf5_field,
				      std::string attr_name,
				      std::string attr_type,
				      std::string attr_center,
				      std::string data_type,
				      int data_dim_row,
				      int data_dim_col
 				    );
  
  void PrintXDMFTopology(std::ofstream& outfstream,
				     std::string hdf5_file,
				     std::string hdf5_field,
				     std::string top_type,
				     int top_dim,
				     int datadim_n_elems,
				     int datadim_el_nodes
				    );  
  
  void PrintXDMFGeometry(std::ofstream& outfstream,
				     std::string hdf5_file,
				     std::string hdf5_field,
				     std::string coord_lev,
				     std::string geom_type,
				     std::string data_type,
				     int data_dim_one,
				     int data_dim_two); 
  
  
}//end namespace IO



#endif