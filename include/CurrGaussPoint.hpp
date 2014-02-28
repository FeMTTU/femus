#ifndef __currgausspoint_h__
#define __currgausspoint_h__

#include "Typedefs_conf.hpp"
#include "FEType_enum.hpp"
#include "VBType_enum.hpp"
#include "CurrGaussPointBase.hpp"



//TODO maybe this gausspoint needs the CurrElement also
//Basically we will have 
// a Current Element, 
// an Abstract Mathematical Element (FEElem),
// an Abstract Geometrical Element (GeomEl)
   
template <unsigned int FM_DIM>
  class CurrGaussPoint : public CurrGaussPointBase {
    
  public:
    
     CurrGaussPoint( EquationsMap& e_map_in );
    ~CurrGaussPoint();
 
double        JacVectVV_g(const uint vbflag, QuantityLocal& xyz )/*const*/;  //TODO should be only for VOLUME
double        JacVectBB_g(const uint vbflag, QuantityLocal& xyz )/* const*/;  //TODO should be only for BOUNDARY 


void         SetPhiElDofsFEVB_g(const uint vbflag,const uint qlflag, const uint qp);
void SetDPhiDxezetaElDofsFEVB_g(const uint vbflag,const uint qlflag, const uint qp);
void    SetDPhiDxyzElDofsFEVB_g(const uint vbflag,const uint qlflag, const uint qp);
void ExtendDphiDxyzElDofsFEVB_g(const uint vbflag,const uint qlflag/*, const uint qp*/);
   
  protected: 
    
        double  _dxyzdxieta_g[FM_DIM - 1][FM_DIM];  //TODO TODO TODO this should stay TOGETHER with ALL THE OTHER DATA that now are in the BASE gauss point CLASS!!!                                                   //IF YOU PUT all the allocations HERE as STATIC it should be faster (maybe...)
        double _dxyzdxezeta_g[FM_DIM][FM_DIM]; //first index: x, second index: csi

      double ComputeJacVV();
        void ComputeJacBB();
	
  };
  
  //se metto alcune funzioni inline QUI queste potrebbero avere una IMPLICIT INSTANTIATION... ma se ne metto almeno una nel SOURCE 
  //allora devo per forza fare una EXPLICIT INSTANTIATION...

#endif