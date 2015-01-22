#ifndef __mgsosadj0_h__
#define __mgsosadj0_h__


// Local Includes
#include "SystemTwo.hpp"

namespace femus {


// Forwarded classes
class MultiLevelProblemTwo;



class EqnNSAD : public SystemTwo {

  public:

     EqnNSAD(  std::vector<Quantity*> int_map_in,
	           MultiLevelProblemTwo& mg_equations_map_in,
                   std::string eqname_in="Eqn_NSAD",
                   std::string varname_in="lambda");

  ~EqnNSAD();

  void ic_read(const double * xp, double * u_value, const double * el_xm) const;

  void  bc_read(const double * xp,const double * normal, int * bc) const;

  void elem_bc_read(const double xp[],int& surf_id,double normal[],int bc_flag[]) const;

  
 void GenMatRhs(const uint Level);



};



} //end namespace femus


#endif
