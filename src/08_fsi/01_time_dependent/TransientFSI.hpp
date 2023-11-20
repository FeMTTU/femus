#ifndef __femus_equations_TransientSystem_FSI_hpp__
#define __femus_equations_TransientSystem_FSI_hpp__

#include "TransientSystem.hpp"

namespace femus {
    
class MonolithicFSINonLinearImplicitSystem;

typedef TransientSystem<MonolithicFSINonLinearImplicitSystem> TransientMonolithicFSINonlinearImplicitSystem;

}


#endif
