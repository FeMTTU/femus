#ifndef __femus_mesh_GeomElemHex_hpp__
#define __femus_mesh_GeomElemHex_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemHex : public GeomElemBase  {

public:
         
    unsigned int  get_dimension() const { return 3; };
    unsigned int n_nodes_linear() const { return 8; };
    

};


} //end namespace femus



#endif
