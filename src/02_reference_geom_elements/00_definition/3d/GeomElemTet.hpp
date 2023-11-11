#ifndef __femus_mesh_GeomElemTet_hpp__
#define __femus_mesh_GeomElemTet_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemTet : public GeomElemBase  {

public:
         
    unsigned int  get_dimension() const { return 3; };
    unsigned int n_nodes_linear() const { return 4; };


};


} //end namespace femus



#endif
