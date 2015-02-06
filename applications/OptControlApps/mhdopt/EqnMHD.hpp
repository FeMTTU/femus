#ifndef __mgsolvermhd0_h__
#define __mgsolvermhd0_h__


#include "SystemTwo.hpp"

namespace femus {

//Forward declarations
class MultiLevelProblem;


class EqnMHD : public SystemTwo {

  public:
  
  EqnMHD( MultiLevelProblem& mg_equations_map_in,
          const std::string & eqname_in,const unsigned int number, const MgSmoother & smoother_type);
	   
  ~EqnMHD();

  void GenMatRhs(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assembe_matrix);
  
};



} //end namespace femus



#endif
