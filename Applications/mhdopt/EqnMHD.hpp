#ifndef __mgsolvermhd0_h__
#define __mgsolvermhd0_h__


#include "EqnBase.hpp"

//Forward declarations
class EquationsMap;


class EqnMHD : public EqnBase {

  public:
  
  EqnMHD(
    std::vector<Quantity*> int_map_in,
    EquationsMap& mg_equations_map_in,
    std::string eqname_in="Eqn_MHD",
    std::string varname_in="B");
	   
  ~EqnMHD();

  void ic_read(double xp[],double u_value[], double el_xm[]);

  void bc_read(double xp[],double normal[],int u_value[]);
  
  void elem_bc_read(double el_xm[],int& surf_id, double value[],int el_flag[]);
  
  void GenMatRhsVB(const uint vb,const double time,const uint Level);
  
   
  
};


#endif
