#ifndef __currgausspointBASE_h__
#define __currgausspointBASE_h__

#include "Typedefs.hpp"
#include "FETypeEnum.hpp"
#include "VBTypeEnum.hpp"
#include "GaussPoints.hpp"
#include "CurrentElem.hpp"

namespace femus {


class CurrentElem;
class elem_type;
class CurrentQuantity;
 

//TODO maybe this gausspoint need the CurrentElement also
//Basically we will have 
// a Current Element, 
// an Abstract Mathematical Element (FEElem),
// an Abstract Geometrical Element (GeomEl)


  class CurrentGaussPointBase {
  
  public:
    
    CurrentGaussPointBase(const CurrentElem & curr_el_in, const Gauss & qrule_in );
   ~CurrentGaussPointBase();
     
   inline const CurrentElem & GetCurrentElem() const {
     return _current_elem;
   }
   
   double**  get_tangent_ptr();   //TODO should be only for BOUNDARY 
   double*   get_normal_ptr();   //TODO should be only for BOUNDARY 
virtual double        JacVectVV_g(CurrentQuantity& xyz )/*const*/ = 0;  //TODO should be only for VOLUME
virtual double        JacVectBB_g(CurrentQuantity& xyz )/* const*/ = 0;  //TODO should be only for BOUNDARY 

virtual void         SetPhiElDofsFEVB_g(const uint qlflag, const uint qp) = 0;
virtual void SetDPhiDxezetaElDofsFEVB_g(const uint qlflag, const uint qp) = 0;
virtual void    SetDPhiDxyzElDofsFEVB_g(const uint qlflag, const uint qp) = 0;
virtual void ExtendDphiDxyzElDofsFEVB_g(const uint qlflag/*, const uint qp*/) = 0;

// TODO encapsulation for phi
inline double Phi(const uint ql,const uint dof) const { 
  return _phi_ndsQLVB_g[ql][dof];
}

  static CurrentGaussPointBase & build(const CurrentElem & elem_in,  const Gauss & qrule_in);

    double* _dphidxezeta_ndsQLVB_g[QL];  //canonical derivatives
    double*    _dphidxyz_ndsQLVB_g[QL];  //physical derivatives
    double*  _dphidxyz_ndsQLVB_g3D[QL];  //physical derivatives in 3D
    double*         _phi_ndsQLVB_g[QL];  //canonical functions
  
  protected:
    
   uint                   _IntDim[VB];   // = {dimension,dimension-1};  //  the dimension of the domain where you integrate based on vb  //TODO is here the correct place?!?
   const CurrentElem & _current_elem;
   const std::vector<const elem_type*>  &  _elem_type;  //[QL]
   const Gauss & _qrule;
   double**  _InvJac_g;/*[FM_DIM][FM_DIM]*/   //gauss Current Geometric Element //no new is needed here //TODO this is only for VOLUME elements
   double*  _normal_g;/*[FM_DIM]*/              //gauss Current Geometric Element //no new is needed here //TODO this should be only for a BOUNDARY GAUSS POINT
   double** _tangent_g;/*[FM_DIM - 1][FM_DIM]*/   //gauss  Current Geometric Element //no new is needed here //TODO this should be only for a BOUNDARY GAUSS POINT
 
  };
  


inline double** CurrentGaussPointBase::get_tangent_ptr() { return _tangent_g;}

inline double*  CurrentGaussPointBase::get_normal_ptr() { return _normal_g;}



} //end namespace femus



#endif