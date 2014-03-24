#include "Physics.hpp"


#include "MeshTwo.hpp"
#include "Utils.hpp"
#include "Files.hpp"

#include "GeomEl.hpp"


// =================================================
/// This function builds the class
Physics::Physics(Utils& mgutils_in):
    _utils(mgutils_in),
    _physrtmap("Physics",mgutils_in._files.get_basepath())  {  }


