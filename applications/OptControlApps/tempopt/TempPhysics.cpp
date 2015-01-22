#include <cmath>

//library
#include "MultiLevelProblemTwo.hpp"
#include "TimeLoop.hpp"
#include "MeshTwo.hpp"
#include "Box.hpp"

//application
#include "Temp_conf.hpp"
#include "TempPhysics.hpp"
#include "EqnT.hpp"

namespace femus {

TempPhysics::TempPhysics( FemusInputParser<double> & map_in):
  Physics(map_in) { 
  
  
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


  
  