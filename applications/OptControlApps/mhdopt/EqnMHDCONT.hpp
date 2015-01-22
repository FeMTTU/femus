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
    
   
EqnMHDCONT(  std::vector<Quantity*> int_map_in,
	     MultiLevelProblemTwo& mg_equations_map_in,
                   std::string eqname_in="Eqn_MHDCONT",
                   std::string varname_in="Becont");

  ~EqnMHDCONT();

  void ic_read(const double * xp, double * u_value, const double * el_xm) const;
  
  void bc_read(const double * xp,const double * normal, int * bc) const;

  void elem_bc_read(const double el_xm[],int& surf_id, double value[],int el_flag[]) const;
  
  void GenMatRhs(const uint Level);

};



} //end namespace femus


#endif
