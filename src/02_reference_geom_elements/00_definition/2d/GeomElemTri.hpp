#ifndef __femus_mesh_GeomElemTri_hpp__
#define __femus_mesh_GeomElemTri_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemTri : public GeomElemBase  {

public:
         
    unsigned int  get_dimension() const { return 2; };
    unsigned int n_nodes_linear() const { return 3; };
    

};


} //end namespace femus



#endif
