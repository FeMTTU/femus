#ifndef __equationt_h__
#define __equationt_h__


// Local Includes -----------------
#include "SystemTwo.hpp"

namespace femus {

// Forward declarations ----------
class MultiLevelProblemTwo;


class EqnT : public SystemTwo {

public:

  
    EqnT( std::vector<Quantity*> int_map_in,
	  MultiLevelProblemTwo& mg_equations_map,
          std::string eqname_in);

    ~EqnT();  
     
 void GenMatRhs(const uint Level);
 
 void elem_bc_read(const double */*el_xm*/, int& surf_id, double *value,int* el_flag) const {};

};


} //end namespace femus



#endif