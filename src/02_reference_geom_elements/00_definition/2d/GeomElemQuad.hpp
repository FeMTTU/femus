#ifndef __femus_mesh_GeomElemQuad_hpp__
#define __femus_mesh_GeomElemQuad_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemQuad : public GeomElemBase  {

public:
         
    unsigned int  get_dimension() const { return 2; };
    unsigned int n_nodes_linear() const { return 4; };
    

};


} //end namespace femus



#endif
