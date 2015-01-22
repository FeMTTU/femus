#ifndef __mg_system2__
#define __mg_system2__

// c/c++
#include <map>
#include <string>

// configure files
#include "Typedefs.hpp"
#include "FemusInputParser.hpp"


namespace femus {



// Forward class
class MeshTwo;

// ========================================
//            Physics class
// ========================================
//this is the base physics class



class Physics {

public:

  FemusInputParser<double>     _physrtmap;

  //mesh
  MeshTwo *               _mesh;   ///< Mesh  pointer
  inline void   set_mesh(MeshTwo * mgmesh_in){    _mesh = mgmesh_in;   return;  }   /// Set mesh pointer

  // Costructor and destructor ----------------------------
   Physics(FemusInputParser<double> & physmap_in); ///< Constructor
  ~Physics(){};                                     ///< Destructor
   
};


} //end namespace femus



#endif
