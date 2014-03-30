#ifndef __mgsolvermhdadj0_h__
#define __mgsolvermhdadj0_h__


//inherited classes
#include "EqnBase.hpp"

// Forwarded classes
class EquationsMap;

class EqnMHDAD : public EqnBase {

  public:
    
EqnMHDAD(  std::vector<Quantity*> int_map_in,
	           EquationsMap& mg_equations_map_in,
                   std::string eqname_in="Eqn_MHDAD",
                   std::string varname_in="xi");

  ~EqnMHDAD();

  void ic_read(double xp[],double u_value[], double el_xm[]);

  void bc_read(double xp[],double normal[],int u_value[]);

 void elem_bc_read(double xp[],int& surf_id,double normal[],int bc_flag[]);

 void GenMatRhsVB(const uint vb,const double time,const uint Level);  ///< Volume Assemblying.
  


};


#endif
