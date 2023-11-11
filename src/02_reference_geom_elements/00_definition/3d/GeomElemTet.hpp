#ifndef __femus_mesh_GeomElemTet_hpp__
#define __femus_mesh_GeomElemTet_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemTet : public GeomElemBase  {

public:
         
    unsigned int  get_dimension() const { return _dim; };
    unsigned int n_nodes_linear() const { return _n_vertices; };


   static constexpr unsigned int _dim        =  3;
   static constexpr unsigned int _n_vertices =  4;

};


} //end namespace femus



#endif
