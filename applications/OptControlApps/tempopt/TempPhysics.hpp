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
  
  
 // =========== Nondimensional groups and ref values ========
  double _pref;  ///< Reference pressure
  double _Re;    ///< Reynolds number
  double _Fr;    ///< Froud Number
  double _Pr;    ///< Prandl Number
   void set_nondimgroups();      /// set nondimensional groups and ref values
  //======================================================

   //========= CONTROL ====================================
      int ElFlagControl(const std::vector<double> el_xm) const;

};



} //end namespace femus


  
#endif