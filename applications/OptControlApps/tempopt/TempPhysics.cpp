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



  return;
}
 


} //end namespace femus


  
  