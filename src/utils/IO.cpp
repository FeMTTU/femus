// Class includes
#include "IO.hpp"

// std libraries 
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>


namespace femus {


namespace IO {
// =============================================================
hid_t read_Dhdf5(hid_t file,const std::string & name,double* data) {

  hid_t  dataset = H5Dopen(file,name.c_str(), H5P_DEFAULT);
  hid_t status=H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

// ===========================================================================
hid_t read_Ihdf5(hid_t file,const std::string & name,int* data) {

  hid_t  dataset = H5Dopen(file,name.c_str(), H5P_DEFAULT);
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

// ===========================================================================
hid_t read_UIhdf5(hid_t file,const std::string & name,uint* data) {

  hid_t  dataset = H5Dopen(file,name.c_str(), H5P_DEFAULT);
  hid_t status=H5Dread(dataset,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

//=========================================
hid_t print_Dhdf5(hid_t file,const std::string & name, hsize_t dimsf[],double* data) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_DOUBLE,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return status;
}

/// Print int data into dhdf5 file
//TODO can we make a TEMPLATE function that takes either "double" or "int" and uses
//either H5T_NATIVE_DOUBLE or H5T_NATIVE_INT? Everything else is the same
hid_t print_Ihdf5(hid_t file,const std::string & name, hsize_t dimsf[],int* data) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return status;
}

//H5T_NATIVE_UINT
hid_t print_UIhdf5(hid_t file,const std::string & name, hsize_t dimsf[],uint* data) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_UINT,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return status;
}


void PrintXDMFAttribute(std::ofstream& outstream, 
				      std::string hdf5_filename, 
				      std::string hdf5_field,
				      std::string attr_name,
				      std::string attr_type,
				      std::string attr_center,
				      std::string data_type,
				      int data_dim_row,
				      int data_dim_col
 				    ) {

                    outstream << "<Attribute Name=\""<< attr_name << "\" "
                              << " AttributeType=\"" << attr_type << "\" "
			      << " Center=\"" <<  attr_center << "\">\n";
                    outstream << "<DataItem  DataType=\""<< data_type << "\" "
                              << " Precision=\"" << 8 << "\" " 
			      << " Dimensions=\"" << data_dim_row << " " << data_dim_col << "\" "
			      << " Format=\"HDF\">  \n";
                    outstream << hdf5_filename  << ":" << hdf5_field << "\n";
                    outstream << "</DataItem>\n" << "</Attribute>\n";
  
  return;
}


void PrintXDMFTopology(std::ofstream& outfstream,
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


void PrintXDMFGeometry(std::ofstream& outfstream,
				     std::string hdf5_file,
				     std::string hdf5_field,
			             std::string coord_lev,
				     std::string geom_type,
				     std::string data_type,
				     int data_dim_one,
				     int data_dim_two) {

    outfstream << "<Geometry Type=\"" << geom_type << "\"> \n";
    for (uint ix=1; ix<4; ix++) {
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

} //end namespace IO


} //end namespace femus


