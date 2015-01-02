#ifndef __mgsolvermhdadj0_h__
#define __mgsolvermhdadj0_h__


//inherited classes
#include "EqnBase.hpp"

namespace femus {

// Forwarded classes
class EquationsMap;

class EqnMHDAD : public EqnBase {

  public:
    
EqnMHDAD(  std::vector<Quantity*> int_map_in,
	           EquationsMap& mg_equations_map_in,
                   std::string eqname_in="Eqn_MHDAD",
                   std::string varname_in="xi");

  ~EqnMHDAD();

 void ic_read(const double * xp, double * u_value, const double * el_xm) const;

 void bc_read(const double * xp,const double * normal, int * bc) const;

 void elem_bc_read(const double xp[],int& surf_id,double normal[],int bc_flag[]) const;

 void GenMatRhsVB(const uint vb, const uint Level);  ///< Volume Assemblying.
  


};



} //end namespace femus


#endif
