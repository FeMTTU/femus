#ifndef __mgsolver_h__
#define __mgsolver_h__


//Inherited classes
#include "EqnBase.hpp"


namespace femus {


// Forwarded classes
class EquationsMap;


class EqnNS : public EqnBase {

  public:

   const uint   _AdvPic_fl;
   const uint   _AdvNew_fl;
   const uint   _Stab_fl;
   const double _Komp_fac;
  
  EqnNS(   std::vector<Quantity*> int_map_in,
	   EquationsMap& mg_equations_map,
           std::string eqname_in="Eqn_NS",/*"Navier-Stokes"*/
           std::string varname_in="u");
   

  ~EqnNS();


 void ic_read(const double * xp, double * u_value, const double * el_xm) const;
 
void  bc_read(const double * xp,const double * normal, int * bc) const;

 void elem_bc_read(const double xp[],int& surf_id,double normal[],int bc_flag[]) const;

 void GenMatRhsVB(const uint vb, const uint Level);  ///< Volume Assemblying.
 
 void ConvertMyselfToChild(EqnBase* mybase);
 
 double ComputeIntegral (const uint vb, const uint Level) ;  //cannot make it const because of set_phi_nds

 
 
};


} //end namespace femus


#endif
