#ifndef __mgsolver_h__
#define __mgsolver_h__


//Inherited classes
#include "SystemTwo.hpp"

namespace femus {

// Forwarded classes
class MultiLevelProblem;


class EqnNS : public SystemTwo {

  public:

   const uint   _AdvPic_fl;
   const uint   _AdvNew_fl;
   const uint   _Stab_fl;
   const double _Komp_fac;
  
  EqnNS(MultiLevelProblem & mg_equations_map,
        const std::string & eqname_in,const unsigned int number, const MgSmoother & smoother_type);

  ~EqnNS();

 void GenMatRhs(const uint Level);
 
 
};


} //end namespace femus





#endif
