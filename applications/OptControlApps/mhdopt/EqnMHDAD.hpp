#ifndef __mgsolvermhdadj0_h__
#define __mgsolvermhdadj0_h__


//inherited classes
#include "SystemTwo.hpp"

namespace femus {

// Forwarded classes
class MultiLevelProblem;

class EqnMHDAD : public SystemTwo {

  public:
    
EqnMHDAD(   MultiLevelProblem& mg_equations_map_in,
            const std::string & eqname_in, const unsigned int number, const MgSmoother & smoother_type);

  ~EqnMHDAD();

 void GenMatRhs(const uint Level);

};



} //end namespace femus


#endif
