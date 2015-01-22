#ifndef __tempphysics__
#define __tempphysics__

#include <vector>
// this function is fully user

#include "Physics.hpp"
#include "FemusInputParser.hpp"

namespace femus {

class MultiLevelProblemTwo;


class MeshTwo;

class TempPhysics : public Physics {

public:
  
//constructor
TempPhysics(  FemusInputParser<double> & map_in); 
  
  
 // =========== Nondimensional groups and ref values ========
  double _pref;  ///< Reference pressure
  double _Re;    ///< Reynolds number
  double _Fr;    ///< Froud Number
  double _Pr;    ///< Prandl Number
   void set_nondimgroups();      /// set nondimensional groups and ref values
  //======================================================


};



} //end namespace femus


  
#endif