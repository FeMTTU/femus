#ifndef __mgsolvermhdadj0_h__
#define __mgsolvermhdadj0_h__


//inherited classes
#include "SystemTwo.hpp"

namespace femus {

// Forwarded classes
class MultiLevelProblemTwo;

class EqnMHDAD : public SystemTwo {

  public:
    
EqnMHDAD(   MultiLevelProblemTwo& mg_equations_map_in,
            const std::string & eqname_in, const unsigned int number, const MgSmoother & smoother_type);

  ~EqnMHDAD();

 void elem_bc_read(const double xp[],int& surf_id,double normal[],int bc_flag[]) const {};

 void GenMatRhs(const uint Level);
  


};



} //end namespace femus


#endif
