#ifndef __equationfetest_h__
#define __equationfetest_h__


// Local Includes -----------------
#include "SystemTwo.hpp"


namespace femus {

// Forward declarations ----------
class MultiLevelProblem;


class EqnT : public SystemTwo {

public:

    EqnT( MultiLevelProblem& mg_equations_map,
          const std::string & eqname_in, const unsigned int number, const MgSmoother & smoother_type);

    ~EqnT();  
     

 void GenMatRhs(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assembe_matrix);
  
 
};


} //end namespace femus



#endif