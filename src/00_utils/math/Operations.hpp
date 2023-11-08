#ifndef __femus_utils_Math_hpp__
#define __femus_utils_Math_hpp__


#include "TypeDefsAndTypeCasts.hpp"

#include <vector>



namespace femus {

// remember that you have to declare all these functions "inline", otherwise you get "multiple definition" in linking

 // Operations ---------------------------------

namespace Math {
  
   inline void zeroN(double* x,const uint N); 
   inline double dotN(const double* x,const double* y,const uint N);  //TODO this will be deleted 
   inline double dot(const double* x,const double* y, const uint spacedim);
   inline double dot(const std::vector< double > x, const std::vector< double > y, const uint spacedim);
   inline std::vector< double > tangent_vector_from_normal(const std::vector< double > orig_vector, const std::vector< double > normal, const uint spacedim);
   inline  void cross(const double* a,const double* b, double* res);
   inline void extend(const double* a, double* a3D, const uint spacedim);
   inline void extend_nds(const uint,const double*, double*, const uint spacedim);
   inline void normalize(double* x,const double fac, const uint spacedim);



// ============== inline functions - BEGIN ==========================   
   
/// set to zero - n components
inline void zeroN(double* x,const uint N)  {

  for (uint i=0; i< N; i++)  x[i]=0.;
  return;
}


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

inline double dot(const std::vector< double > x,const std::vector< double > y, const uint spacedim) {
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



   inline std::vector< double > tangent_vector_from_normal(const std::vector< double > orig_vector, const std::vector< double > normal, const uint spacedim) {
     
      const double  orig_vector_dot_normal = Math::dot(orig_vector, normal, orig_vector.size());

      std::vector < double > orig_vector_dot_tangent(orig_vector.size(), 0.);
            for (unsigned d = 0; d < orig_vector_dot_tangent.size() ; d++) {
              orig_vector_dot_tangent[d] = orig_vector[d] - orig_vector_dot_normal * normal[d];   // a - (a.n) n
            }
     
    return orig_vector_dot_tangent;
   }

// ============== inline functions - END ==========================


 
} //end namespace Math




} //end namespace femus



#endif
