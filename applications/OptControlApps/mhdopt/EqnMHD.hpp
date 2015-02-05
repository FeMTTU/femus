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

  void GenMatRhs(const uint Level);
  
};



} //end namespace femus



#endif
