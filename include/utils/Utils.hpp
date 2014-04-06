#ifndef __utils_h__
#define __utils_h__

// external libs
#include <map>
#include <string>
#include "hdf5.h"

// Femus
#include "Typedefs.hpp"
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


};


namespace Math {
 // Operations ---------------------------------
  static inline void zeroN(double* x,const uint N); 
  static inline double dotN(const double* x,const double* y,const uint N);  //TODO this will be deleted 
  static inline double dot(const double* x,const double* y, const uint spacedim);
  static inline  void cross(const double* a,const double* b, double* res);
  static inline void extend(const double* a, double* a3D, const uint spacedim);
  static inline void extend_nds(const uint,const double*, double*, const uint spacedim);
  static inline void normalize(double* x,const double fac, const uint spacedim);

/// set to zero - n components
inline void zeroN(double* x,const uint N)  {

  for (uint i=0; i< N; i++)  x[i]=0.;
  return;
}


//  useful  functions --------------------------------------------
/// dot product
inline void normalize(double* x,const double fac, const uint spacedim)  {
  for (uint idim=0; idim< spacedim; idim++)  x[idim] /= fac;
  return;
}

/// dot product - n components
inline double dotN(const double* x,const double* y,const uint N)  {
  double dotprod=0.;
  for (uint idim=0; idim< N; idim++)  dotprod += x[idim]*y[idim];
  return dotprod;
}

/// dot product
inline double dot(const double* x,const double* y, const uint spacedim) {
  double dotprod=0.;
  for (uint idim=0; idim < spacedim; idim++)  dotprod += x[idim]*y[idim];
  return dotprod;
}

/// Cross product
inline void cross(const double* a,const double* b, double* res) {
//a,b,res are 3D vectors
//clean then fill
  for (uint i=0; i<3; i++) res[i]=0.;
  for (uint i=0; i<3; i++) res[i] = (a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3]);
  return;
}

/// extend to a 3D vector a vector with dimension 
inline void extend(const double* a, double* a3D, const uint spacedim)  {
  for (uint i=0; i<3; i++) a3D[i]=0.;
  for (uint i=0; i< spacedim; i++)  a3D[i]=a[i];
  return;
}

/// extend to 3D an element dof vector
inline void extend_nds(const uint el_ndofs,const double* a_nds, double* a_nds3D, const uint spacedim)  {

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



} //end namespace Math



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


//  std::map<std::string,double>& get_utils_map() const   {return _utilsmap;}
// error: invalid initialization of reference of type ‘std::__debug::map<std::basic_string<char>, double>&’ from expression of type ‘const std::__debug::map<std::basic_string<char>, double>’
//this function must return a REFERENCE to the object _utilsmap in this class, NOT A COPY...
//libmesh does like that:
//   std::map<std::string,double>& get_utils_map()    {return _utilsmap;}
// in this case you give a reference,so if something changes in the origin then you can get it
