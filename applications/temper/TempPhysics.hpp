#ifndef __tempphysics__
#define __tempphysics__

#include <vector>
// this function is fully user

#include "Physics.hpp"
class EquationsMap;


class Utils;
class Mesh;

class TempPhysics : public Physics {

public:
  
//constructor
TempPhysics(Utils& mgutils_in); 
  
  
 // =========== Nondimensional groups and ref values ========
  double _pref;  ///< Reference pressure
  double _Re;    ///< Reynolds number
  double _Fr;    ///< Froud Number
  double _Pr;    ///< Prandl Number
   void set_nondimgroups();      /// set nondimensional groups and ref values
  //======================================================

   //========= CONTROL ====================================
      int ElFlagControl(const std::vector<double> el_xm) const;

   //================= ALGORITHM ==========================
  void transient_loopPlusJ(EquationsMap & eqmap_in);
  

};

  
#endif