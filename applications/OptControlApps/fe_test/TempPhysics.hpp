#ifndef __tempphysics__
#define __tempphysics__

#include <vector>
// this function is fully user

#include "Physics.hpp"
#include "RunTimeMap.hpp"

namespace femus {

class EquationsMap;


class Utils;
class MeshTwo;

class TempPhysics : public Physics {

public:
  
//constructor
TempPhysics(  RunTimeMap<double> & map_in); 
  
   //========= CONTROL ====================================
      int ElFlagControl(const std::vector<double> el_xm) const;

};



} //end namespace femus


  
#endif