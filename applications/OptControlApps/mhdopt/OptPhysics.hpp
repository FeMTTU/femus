#ifndef __optphysics__
#define __optphysics__

// this function is fully user

#include <vector>
#include "Physics.hpp"


namespace femus {


class Utils;
class MeshTwo;

class OptPhysics : public Physics {

public:
  
//constructor
OptPhysics( RunTimeMap<double> & map_in); 
  
  
 // =========== Nondimensional groups and ref values ========
  double _pref;  ///< Reference pressure
  double _Re;    ///< Reynolds number
  double _Fr;    ///< Froud Number
  double _S;     ///< Coupling Magnetic coeff
  double _Pr;    ///< Prandl Number
  double _Rem;   ///< Magnetic Reynolds number
  double _Hm;    ///< Hartmann Number
  double _We;    ///< Weber Number
   void set_nondimgroups();      /// set nondimensional groups and ref values
  //======================================================

  //========= CONTROL ====================================
      int ElFlagControl(const std::vector<double> el_xm) const;
  //======================================================
 
   
};
  


} //end namespace femus


  
  
#endif