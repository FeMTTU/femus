#ifndef __currgausspointBASE_h__
#define __currgausspointBASE_h__


#include "Typedefs_conf.hpp"
#include "FEType_enum.hpp"
#include "VBType_enum.hpp"

class FEElemBase;
class QRule;
class QuantityLocal;
class EquationsMap;
 

//TODO maybe this gausspoint need the CurrElement also
//Basically we will have 
// a Current Element, 
// an Abstract Mathematical Element (FEElem),
// an Abstract Geometrical Element (GeomEl)


  class CurrGaussPointBase {
  
  protected:
    
   double**  _InvJac_g;/*[FM_DIM][FM_DIM]*/   //gauss Current Geometric Element //no new is needed here //TODO this is only for VOLUME elements
   double*  _normal_g;/*[FM_DIM]*/              //gauss Current Geometric Element //no new is needed here //TODO this should be only for a BOUNDARY GAUSS POINT
   double** _tangent_g;/*[FM_DIM - 1][FM_DIM]*/   //gauss  Current Geometric Element //no new is needed here //TODO this should be only for a BOUNDARY GAUSS POINT
 
  public:
    
    CurrGaussPointBase( EquationsMap& e_map_in );
   ~CurrGaussPointBase();
 
    EquationsMap         & _eqnmap;
    FEElemBase*         _AbsFEVect[QL];
    QRule&   _qrule;
    uint                   _IntDim[VB];   // = {dimension,dimension-1};  //  the dimension of the domain where you integrate based on vb  //TODO is here the correct place?!?
    double*         _phi_ndsQLVB_g[VB][QL];  //canonical functions  //TODO here it seems to contain GAUSS x ELDOFS
    double* _dphidxezeta_ndsQLVB_g[VB][QL];  //canonical derivatives
    double*    _dphidxyz_ndsQLVB_g[VB][QL];  //physical derivatives
    double*  _dphidxyz_ndsQLVB_g3D[VB][QL];  //physical derivatives in 3D

 double**  get_tangent_ptr();
 double*   get_normal_ptr();
virtual double        JacVectVV_g(const uint vbflag, QuantityLocal& xyz )/*const*/ = 0;  //TODO should be only for VOLUME
virtual double        JacVectBB_g(const uint vbflag, QuantityLocal& xyz )/* const*/ = 0;  //TODO should be only for BOUNDARY 

virtual void         SetPhiElDofsFEVB_g(const uint vbflag,const uint qlflag, const uint qp) = 0;
virtual void SetDPhiDxezetaElDofsFEVB_g(const uint vbflag,const uint qlflag, const uint qp) = 0;
virtual void    SetDPhiDxyzElDofsFEVB_g(const uint vbflag,const uint qlflag, const uint qp) = 0;
virtual void ExtendDphiDxyzElDofsFEVB_g(const uint vbflag,const uint qlflag/*, const uint qp*/) = 0;


  static CurrGaussPointBase& build(EquationsMap& e_map_in, const uint dim);  //Let us try with REFERENCE instead of POINTER

  };
  


inline double** CurrGaussPointBase::get_tangent_ptr() { return _tangent_g;}

inline double*  CurrGaussPointBase::get_normal_ptr() { return _normal_g;}


#endif