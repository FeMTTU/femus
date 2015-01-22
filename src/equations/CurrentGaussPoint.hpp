#ifndef __currgausspoint_h__
#define __currgausspoint_h__

#include "Typedefs.hpp"
#include "FETypeEnum.hpp"
#include "VBTypeEnum.hpp"
#include "CurrentGaussPointBase.hpp"


namespace femus {





//TODO maybe this gausspoint needs the CurrentElement also
//Basically we will have 
// a Current Element, 
// an Abstract Mathematical Element (FEElem),
// an Abstract Geometrical Element (GeomEl)

// This gauss point receives the Current Element, which can be a volume or boundary (or even less) element.
// So far this class has the VV and BB parts altogether 
// Questa classe ha potenzialmente sia la parte VV sia la parte BB, ma all'atto della costruzione  
// ne sceglie solo uno, e alloca le strutture solo per quello. In questo modo solo le funzioni VV o BB dovranno essere usate.  


  
template <unsigned int FM_DIM>
  class CurrentGaussPoint : public CurrentGaussPointBase {
    
  public:
    
     CurrentGaussPoint(const CurrentElem & curr_el_in, MultiLevelProblemTwo& e_map_in );
    ~CurrentGaussPoint();
 
double        JacVectVV_g(CurrentQuantity& xyz )/*const*/;  //TODO should be only for VOLUME
double        JacVectBB_g(CurrentQuantity& xyz )/* const*/;  //TODO should be only for BOUNDARY 


void         SetPhiElDofsFEVB_g(const uint qlflag, const uint qp);
void SetDPhiDxezetaElDofsFEVB_g(const uint qlflag, const uint qp);
void    SetDPhiDxyzElDofsFEVB_g(const uint qlflag, const uint qp);
void ExtendDphiDxyzElDofsFEVB_g(const uint qlflag/*, const uint qp*/);
   
  protected: 
    
        double  _dxyzdxieta_g[FM_DIM - 1][FM_DIM];  //only for BOUNDARY //TODO TODO TODO this should stay TOGETHER with ALL THE OTHER DATA that now are in the BASE gauss point CLASS!!!                                                   //IF YOU PUT all the allocations HERE as STATIC it should be faster (maybe...)
        double _dxyzdxezeta_g[FM_DIM][FM_DIM]; //only for VOLUME //first index: x, second index: csi

      double ComputeJacVV();     //TODO should be only for VOLUME
        void ComputeJacBB();     //TODO should be only for BOUNDARY 
	
  };
  
  //se metto alcune funzioni inline QUI queste potrebbero avere una IMPLICIT INSTANTIATION... ma se ne metto almeno una nel SOURCE 
  //allora devo per forza fare una EXPLICIT INSTANTIATION...


} //end namespace femus



#endif