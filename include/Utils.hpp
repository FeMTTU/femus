#ifndef __utils_h__
#define __utils_h__

// external libs
#include <map>
#include <string>
#include "hdf5.h"

// Femus
#include "Typedefs_conf.hpp"
#include "RunTimeMap.hpp"


class Files;

// =========================================
//              Utils
// =========================================

class Utils  {

public:
    Files&                 _files;     ///< Files class pointer
    RunTimeMap<double>     _urtmap;
  
  // Constructor-Destructor --------------------
  Utils(Files& mgfiles_in);  ///< Constructor
  ~Utils() {};                 ///< Destructor


 // Operations ---------------------------------
  inline void zeroN(double* x,const uint N) const; 
  inline double dotN(const double* x,const double* y,const uint N) const;  //this will be deleted 
  inline double dot(const double* x,const double* y, const uint spacedim) const;
  inline  void cross(const double* a,const double* b, double* res) const;
  inline void extend(const double* a, double* a3D, const uint spacedim) const;
  inline void extend_nds(const uint,const double*, double*, const uint spacedim) const;
  inline void normalize(double* x,const double fac, const uint spacedim) const;

  //hdf5 ------------------------------------
  hid_t print_Dhdf5(hid_t file,const std::string & name, hsize_t dimsf[],double data[]);
  hid_t print_Ihdf5(hid_t file,const std::string & name, hsize_t dimsf[],int data[]);
  hid_t read_Dhdf5(hid_t file,const std::string & name,double data[]);
  hid_t read_Ihdf5(hid_t file,const std::string & name,int data[]);
  
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
};



/// set to zero - n components
inline void Utils::zeroN(double* x,const uint N) const {

  for (uint i=0; i< N; i++)  x[i]=0.;
  return;
}


//  useful  functions --------------------------------------------
/// dot product
inline void Utils::normalize(double* x,const double fac, const uint spacedim) const {
  for (uint idim=0; idim< spacedim; idim++)  x[idim] /= fac;
  return;
}

/// dot product - n components
inline double Utils::dotN(const double* x,const double* y,const uint N) const {
  double dotprod=0.;
  for (uint idim=0; idim< N; idim++)  dotprod += x[idim]*y[idim];
  return dotprod;
}

/// dot product
inline double Utils::dot(const double* x,const double* y, const uint spacedim) const {
  double dotprod=0.;
  for (uint idim=0; idim < spacedim; idim++)  dotprod += x[idim]*y[idim];
  return dotprod;
}

/// Cross product
inline void Utils::cross(const double* a,const double* b, double* res) const {
//a,b,res are 3D vectors
//clean then fill
  for (uint i=0; i<3; i++) res[i]=0.;
  for (uint i=0; i<3; i++) res[i] = (a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3]);
  return;
}

/// extend to a 3D vector a vector with dimension 
inline void Utils::extend(const double* a, double* a3D, const uint spacedim) const {
  for (uint i=0; i<3; i++) a3D[i]=0.;
  for (uint i=0; i< spacedim; i++)  a3D[i]=a[i];
  return;
}

/// extend to 3D an element dof vector
void Utils::extend_nds(const uint el_ndofs,const double* a_nds, double* a_nds3D, const uint spacedim) const {

//AAA: valid from spacedim to 3

//set to zero
  for (uint eln=0; eln<el_ndofs; eln++)  {
    for (uint i=0; i<3; i++) {
      a_nds3D[eln+i*el_ndofs]=0.;
    }
  }
//extend
  for (uint eln=0; eln<el_ndofs; eln++)    {
    for (uint idim=0; idim<spacedim; idim++) {
      a_nds3D[eln+idim*el_ndofs] = a_nds[eln+idim*el_ndofs];
    }
  }

  return;
}



#endif


//  std::map<std::string,double>& get_utils_map() const   {return _utilsmap;}
// error: invalid initialization of reference of type ‘std::__debug::map<std::basic_string<char>, double>&’ from expression of type ‘const std::__debug::map<std::basic_string<char>, double>’
//this function must return a REFERENCE to the object _utilsmap in this class, NOT A COPY...
//libmesh does like that:
//   std::map<std::string,double>& get_utils_map()    {return _utilsmap;}
// in this case you give a reference,so if something changes in the origin then you can get it
