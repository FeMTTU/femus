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


OptPhysics::OptPhysics( FemusInputParser<double> & map_in):
  Physics(map_in) { } 




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



