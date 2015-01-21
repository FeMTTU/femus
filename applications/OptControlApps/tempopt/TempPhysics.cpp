#include <cmath>

//library
#include "EquationsMap.hpp"
#include "TimeLoop.hpp"
#include "MeshTwo.hpp"
#include "Box.hpp"

//application
#include "Temp_conf.hpp"
#include "TempPhysics.hpp"
#include "EqnT.hpp"

namespace femus {

TempPhysics::TempPhysics( RunTimeMap<double> & map_in):
  Physics(map_in) { 
  
  
} 

 // ========================================================
 //the Box must already be initialized here
 //it must receive an el_xm already in the REFBox frame

///AAA questa funzione NON lavora se tu fai solo DUE SUDDIVISIONI PER LATO e nolevels=1 !!!

 int TempPhysics::ElFlagControl(const std::vector<double> el_xm)  const {

  Box* box= static_cast<Box*>(_mesh->GetDomain());
   
     int el_flagdom=0;

///optimal control
  #if DIMENSION==2
//============== OUTFLOW     
//        if (   el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
// 	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0])
// 	   && el_xm[1] > 0.75*(box->_le[1] - box->_lb[1]) ) {
//                  el_flagdom=1;
//              }
//============== RIGHT
//      if (   el_xm[0] > 0.85*(box->_le[0] - box->_lb[0])
// 	   && el_xm[1] > 0.25*(box->_le[1] - box->_lb[1])
// 	   && el_xm[1] < 0.75*(box->_le[1] - box->_lb[1]) ) {
//                  el_flagdom=1;
//              }
//============== INFLOW
// facciamo 0.2-0.4; 0.4-0.6; 0.6-0.8     
     if (     el_xm[1] > 0.4*(box->_le[1] - box->_lb[1])
           && el_xm[1] < 0.6*(box->_le[1] - box->_lb[1])
	   && el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0]) ) {
                 el_flagdom=1;
             }
// //============== CENTER
//      if (  /*   el_xm[1] > 0.04*(box->_le[1] - box->_lb[1])*/
//           /* &&*/ el_xm[1] < 0.06*(box->_le[1] - box->_lb[1])
// 	   && el_xm[0] > 0.4*(box->_le[0] - box->_lb[0])
// 	   && el_xm[0] < 0.6*(box->_le[0] - box->_lb[0]) ) {
//                  el_flagdom=1;
//              }
     
 #else
      if ( el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])  
	&& el_xm[0] < 0.75*(box->_le[0] - box->_lb[0]) 
	&& el_xm[1] > 0.75*(box->_le[1] - box->_lb[1])
	&& el_xm[2] > 0.25*(box->_le[2] - box->_lb[2]) 
	&& el_xm[2] < 0.75*(box->_le[2] - box->_lb[2]) ) {
	el_flagdom=1;
        }
 #endif

return el_flagdom; 
}



// ====================================================
/// This function sets the nondimensional groups
void TempPhysics::set_nondimgroups() {

  const double rhof   = _physrtmap.get("rho0");
  const double Uref   = _physrtmap.get("Uref");
  const double Lref   = _physrtmap.get("Lref");
  const double  muf   = _physrtmap.get("mu0");

  _pref = rhof*Uref*Uref;
  _Re  = (rhof*Uref*Lref)/muf;
  _Fr  = (Uref*Uref)/(9.81*Lref);
  _Pr=muf/rhof;    ///< Prandl Number


  return;
}
 


} //end namespace femus


  
  