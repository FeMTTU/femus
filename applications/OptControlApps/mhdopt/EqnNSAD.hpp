#ifndef __mgsosadj0_h__
#define __mgsosadj0_h__


// Local Includes
#include "SystemTwo.hpp"

namespace femus {


// Forwarded classes
class MultiLevelProblem;



class EqnNSAD : public SystemTwo {

  public:

     EqnNSAD(  MultiLevelProblem& mg_equations_map_in,
               const std::string & eqname_in, const unsigned int number, const MgSmoother & smoother_type);

  ~EqnNSAD();

  void GenMatRhs(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assembe_matrix);



};



} //end namespace femus


#endif
