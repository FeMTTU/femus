#ifndef __mgsosadj0_h__
#define __mgsosadj0_h__


// Local Includes
#include "SystemTwo.hpp"

namespace femus {


// Forwarded classes
class MultiLevelProblemTwo;



class EqnNSAD : public SystemTwo {

  public:

     EqnNSAD(  MultiLevelProblemTwo& mg_equations_map_in,
               std::string eqname_in, const unsigned int number, const MgSmoother & smoother_type);

  ~EqnNSAD();

  void elem_bc_read(const double xp[],int& surf_id,double normal[],int bc_flag[]) const {};
  
 void GenMatRhs(const uint Level);



};



} //end namespace femus


#endif
