#include "OptPhysics.hpp"

//C++
#include <cmath>

//library
#include "MeshTwo.hpp"
#include "Box.hpp"
#include "VBTypeEnum.hpp"

//application
#include "Opt_conf.hpp"


namespace femus {


OptPhysics::OptPhysics( RunTimeMap<double> & map_in):
  Physics(map_in) { } 

 // ========================================================
 //the Box must already be initialized here
 //it must receive an el_xm already in the REFBox frame
//  #define LCX1  0.25
// #define LCX2  0.75
// #define LCY1  0.90625
///AAA questa funzione non lavora se tu fai solo DUE SUDDIVISIONI PER LATO e nolevels=1 !!!

 int OptPhysics::ElFlagControl(const std::vector<double> el_xm)  const {

  Box* box= static_cast<Box*>(_mesh->GetDomain());
   
   
     int el_flagdom=0;

///optimal control
  #if DIMENSION==2
   //flag on the controlled region 2D
       if (   el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0])
	   && el_xm[1] > 0.75*(box->_le[1] - box->_lb[1]) ) {
                 el_flagdom=1;
             }
  #else
   //flag on the controlled region 3D
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
void OptPhysics::set_nondimgroups() {

  const double rhof   = _physrtmap.get("rho0");
  const double Uref   = _physrtmap.get("Uref");
  const double Lref   = _physrtmap.get("Lref");
  const double  muf   = _physrtmap.get("mu0");

  _pref = rhof*Uref*Uref;
  _Re  = (rhof*Uref*Lref)/muf;
  _Fr  = (Uref*Uref)/(9.81*Lref);
  _Pr=muf/rhof;    ///< Prandl Number
  const double MUMHD  = _physrtmap.get("MUMHD");
  const double SIGMHD = _physrtmap.get("SIGMHD");
  const double   Bref = _physrtmap.get("Bref");
  _Rem = MUMHD*SIGMHD*Uref*Lref;
  _Hm  = Bref*Lref*sqrt(SIGMHD/muf);
  _S   = _Hm*_Hm/(_Re*_Rem);
  const double sigma  = _physrtmap.get("sigma");
  _We  = (Uref*Uref*Lref*rhof)/sigma;


  return;
}


} //end namespace femus



