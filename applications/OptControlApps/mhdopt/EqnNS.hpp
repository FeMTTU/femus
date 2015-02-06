#ifndef __mgsolver_h__
#define __mgsolver_h__


//Inherited classes
#include "SystemTwo.hpp"


namespace femus {


class EqnNS : public SystemTwo {

  public:
  
  EqnNS(   MultiLevelProblem& mg_equations_map,
           const std::string & eqname_in, const unsigned int number, const MgSmoother & smoother_type);
   

  ~EqnNS();

 void GenMatRhs(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assembe_matrix);
 
};


} //end namespace femus


#endif
