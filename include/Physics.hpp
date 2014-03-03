#ifndef __mg_system2__
#define __mg_system2__

// std libraries ---------------
#include <map>
#include <string>

// configure files ----------------
#include "Typedefs.hpp"
#include "RunTimeMap.hpp"

// Forward class -------------
class Utils;
class Mesh;

// ========================================
//            Physics class
// ========================================
//this is the base physics class



class Physics {

protected:
  
  Utils&                         _utils;   ///< Utils pointer
  
public:

  RunTimeMap<double>     _physrtmap;

  //mesh
  Mesh*                 _mesh;   ///< Mesh  pointer
  inline void   set_mesh(Mesh * mgmesh_in){    _mesh = mgmesh_in;   return;  }   /// Set mesh pointer

  // Costructor and destructor ----------------------------
   Physics(Utils& mgutils_in); ///< Constructor
  ~Physics(){};                                     ///< Destructor
   
};

#endif
