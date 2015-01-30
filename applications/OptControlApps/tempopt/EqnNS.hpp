#ifndef __mgsolver_h__
#define __mgsolver_h__


//Inherited classes
#include "SystemTwo.hpp"

namespace femus {

// Forwarded classes
class MultiLevelProblemTwo;


class EqnNS : public SystemTwo {

  public:

   const uint   _AdvPic_fl;
   const uint   _AdvNew_fl;
   const uint   _Stab_fl;
   const double _Komp_fac;
  
  EqnNS(MultiLevelProblemTwo& mg_equations_map,
        std::string eqname_in,const unsigned int number, const MgSmoother & smoother_type);

  ~EqnNS();

 void elem_bc_read(const double * xp, int& surf_id, double normal[], int bc_flag[]) const;

 void GenMatRhs(const uint Level);
 
 
};


} //end namespace femus





#endif
