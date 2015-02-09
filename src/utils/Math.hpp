#ifndef __math_hpp__
#define __math_hpp__


#include "Typedefs.hpp"

#include "CurrentQuantity.hpp"
#include "ElemType.hpp"
#include "MultiLevelProblem.hpp"

namespace femus {



 // Operations ---------------------------------

namespace Math {
  
   inline void zeroN(double* x,const uint N); 
   inline double dotN(const double* x,const double* y,const uint N);  //TODO this will be deleted 
   inline double dot(const double* x,const double* y, const uint spacedim);
   inline  void cross(const double* a,const double* b, double* res);
   inline void extend(const double* a, double* a3D, const uint spacedim);
   inline void extend_nds(const uint,const double*, double*, const uint spacedim);
   inline void normalize(double* x,const double fac, const uint spacedim);



// ============== inline functions ==========================   
   
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


//=================================================================

inline double FunctionIntegral (const uint vb, MultiLevelProblem & eqnmap_in, double (*pt2func)(double, const std::vector<double> ) )  {

  const uint mesh_vb = vb;
  
const uint Level  = eqnmap_in.GetMeshTwo()._NoLevels - 1;
const uint myproc = eqnmap_in.GetMeshTwo()._iproc;
  double time = 0.;
  
    CurrentElem       currelem(Level,vb,NULL,eqnmap_in.GetMeshTwo(),eqnmap_in.GetElemType()); //element without equation
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,eqnmap_in.GetQrule(currelem.GetDim()));

  //======== ELEMENT MAPPING =======
  const uint meshql = (int) eqnmap_in.GetMeshTwo().GetRuntimeMap().get("meshql");  
 
//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = eqnmap_in.GetMeshTwo().get_dim();
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

  double integral = 0.;
  
//loop over the geom el types
      const uint el_ngauss = eqnmap_in.GetQrule(currelem.GetDim()).GetGaussPointsNumber();

//parallel sum
    const uint nel_e = eqnmap_in.GetMeshTwo()._off_el[mesh_vb][eqnmap_in.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = eqnmap_in.GetMeshTwo()._off_el[mesh_vb][eqnmap_in.GetMeshTwo()._NoLevels*myproc+Level];
  
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {
  
    currelem.SetDofobjConnCoords(myproc,iel);
    currelem.SetMidpoint(); 
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);

     
    for (uint qp = 0; qp < el_ngauss; qp++) {

for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  }  
     
double  Jac_g=0.;
          if (vb==0)   Jac_g = currgp.JacVectVV_g(xyz);  //not xyz_refbox!      
     else if (vb==1)   Jac_g = currgp.JacVectBB_g(xyz);  //not xyz_refbox!      

   const double  wgt_g = eqnmap_in.GetQrule(currelem.GetDim()).GetGaussWeight(qp);

     for (uint fe = 0; fe < QL; fe++)     {          currgp.SetPhiElDofsFEVB_g (fe,qp);  }

 xyz.val_g();
double myval_g = pt2func(time,xyz._val_g); 

 
  integral += wgt_g*Jac_g*myval_g;

   
    }//gauss loop
     
    }//element loop
    
         std::cout << std::endl  << " ^^^^^^^^^^^^^^^^^L'integrale sul processore "<< myproc << " vale: " << integral << std::endl;

    double weights_sum = 0.;
    for (uint qp = 0; qp < el_ngauss; qp++)  weights_sum += eqnmap_in.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
       std::cout << std::endl << " ^^^^^^^^^^^^^^^^^ La somma dei pesi  vale: " << weights_sum << std::endl;

       double J=0.;
   #ifdef HAVE_MPI
//       MPI_Reduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );  //This one gives J only to processor 0 !
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
   #else
   J = integral;
   #endif
    
     std::cout << std::endl << " ^^^^^^^^^^^^^^^^^L'integrale totale vale: " << J << std::endl;

  return J;  
  
}



} //end namespace Math


} //end namespace femus



#endif