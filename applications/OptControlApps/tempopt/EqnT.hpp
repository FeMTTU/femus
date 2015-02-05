#ifndef __equationt_h__
#define __equationt_h__


// Local Includes -----------------
#include "SystemTwo.hpp"

namespace femus {

// Forward declarations ----------
class MultiLevelProblem;


class EqnT : public SystemTwo {

public:

  
    EqnT( MultiLevelProblem & mg_equations_map,
          const std::string & eqname_in, const unsigned int number, const MgSmoother & smoother_type);

    ~EqnT();  
     
 void GenMatRhs(const uint Level);
 
};


} //end namespace femus



#endif