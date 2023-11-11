#ifndef __femus_mesh_GeomElemEdge_hpp__
#define __femus_mesh_GeomElemEdge_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemEdge : public GeomElemBase  {

public:
         
    unsigned int  get_dimension() const { return 1; };
    unsigned int n_nodes_linear() const { return 2; };
    

};


} //end namespace femus



#endif
