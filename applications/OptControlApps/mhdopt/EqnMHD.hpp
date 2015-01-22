#ifndef __mgsolvermhd0_h__
#define __mgsolvermhd0_h__


#include "SystemTwo.hpp"

namespace femus {

//Forward declarations
class MultiLevelProblemTwo;


class EqnMHD : public SystemTwo {

  public:
  
  EqnMHD(
    std::vector<Quantity*> int_map_in,
    MultiLevelProblemTwo& mg_equations_map_in,
    std::string eqname_in="Eqn_MHD",
    std::string varname_in="B");
	   
  ~EqnMHD();

  void ic_read(const double * xp, double * u_value, const double * el_xm) const;

  void bc_read(const double * xp,const double * normal, int * bc) const;
  
  void elem_bc_read(const double el_xm[],int& surf_id, double value[],int el_flag[]) const;
  
  void GenMatRhs(const uint Level);
  
   
  
};



} //end namespace femus



#endif
