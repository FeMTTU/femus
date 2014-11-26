#ifndef __mgsolvermhdcont0_h__
#define __mgsolvermhdcont0_h__

//C++ includes
#include <vector>

// Local Includes
#include "EqnBase.hpp"


namespace femus {
  
// Forwarded classes
class NumericVector;

class EqnMHDCONT : public EqnBase {

  public:
    
//====data  
    std::vector<NumericVector *> _x_oldopt;  //old optimization step
  
    
    
EqnMHDCONT(  std::vector<Quantity*> int_map_in,
	     EquationsMap& mg_equations_map_in,
                   std::string eqname_in="Eqn_MHDCONT",
                   std::string varname_in="Becont");

  ~EqnMHDCONT();

  void ic_read(double xp[],double u_value[], double el_xm[]);

  void bc_read(double xp[],double normal[],int u_value[]);

  void elem_bc_read(double el_xm[],int& surf_id, double value[],int el_flag[]);
  
  void GenMatRhsVB(const uint vb,const double time,const uint Level);

  void init_equation_data();

};



} //end namespace femus


#endif
