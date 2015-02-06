#ifndef __mgsolvermhdcont0_h__
#define __mgsolvermhdcont0_h__

//C++ includes
#include <vector>

// Local Includes
#include "SystemTwo.hpp"


namespace femus {
  
// Forwarded classes
class NumericVector;

class EqnMHDCONT : public SystemTwo {

  public:
    
   
EqnMHDCONT( MultiLevelProblem& mg_equations_map_in,
            const std::string & eqname_in, const unsigned int number, const MgSmoother & smoother_type);

  ~EqnMHDCONT();

  void GenMatRhs(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assembe_matrix);

};



} //end namespace femus


#endif
