#ifndef __femus_mesh_GeomElemEdge_hpp__
#define __femus_mesh_GeomElemEdge_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemEdge : public GeomElemBase  {

public:
         
    unsigned int  get_dimension() const { return _dim; };
    unsigned int n_nodes_linear() const { return _n_vertices; };

   static constexpr unsigned int _dim        =  1;
   static constexpr unsigned int _n_vertices =  2;

};


} //end namespace femus



#endif
