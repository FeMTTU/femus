#ifndef __equationt_h__
#define __equationt_h__


// Local Includes -----------------
#include "EqnBase.hpp"

// Forward declarations ----------
class EquationsMap;


class EqnT : public EqnBase {

public:

    EqnT( std::vector<Quantity*> int_map_in,
	  EquationsMap& mg_equations_map,
          std::string eqname_in="Eqn_T",
          std::string varname_in="T");

    ~EqnT();  
     

 void GenMatRhsVB(const uint vb,const double time,const uint Level);  ///< Volume Assemblying.
 
 void ic_read(double xp[],double u_value[], double el_xm[]);
 
 void bc_read(double xp[],double normal[],int u_value[]);

 double ComputeIntegral (const uint vb, const uint Level);

 double ComputeNormControl (const uint vb, const uint Level, const uint reg_ord );

};

#endif