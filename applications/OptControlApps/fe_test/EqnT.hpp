#ifndef __equationfetest_h__
#define __equationfetest_h__


// Local Includes -----------------
#include "SystemTwo.hpp"


namespace femus {

// Forward declarations ----------
class MultiLevelProblemTwo;


class EqnT : public SystemTwo {

public:

    EqnT( std::vector<Quantity*> int_map_in,
	  MultiLevelProblemTwo& mg_equations_map,
          std::string eqname_in="Eqn_T",
          std::string varname_in="T");

    ~EqnT();  
     

 void GenMatRhs(const uint Level);  ///< Volume Assemblying.
 
 void ic_read(const double * xp, double * u_value, const double * el_xm) const;
 
 void  bc_read(const double * xp,const double * normal, int * bc) const;

 void elem_bc_read(const double */*el_xm*/, int& surf_id, double *value,int* el_flag) const {};
 
};


} //end namespace femus



#endif