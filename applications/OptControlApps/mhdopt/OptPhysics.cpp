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



  return;
}


} //end namespace femus



