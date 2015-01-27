#ifndef __mggeomel_h__
#define __mggeomel_h__


#include <string>

#include "Typedefs.hpp"
#include "FETypeEnum.hpp"


namespace femus {

// Merge this in Elem... but wait, this is all ABSTRACT, what about Elem?

class GeomEl  {

public:

     GeomEl(const std::string geomel_id, const uint mesh_order);
    ~GeomEl();
    
    std::string _geomel_id;      
    uint      _elnds;     //NVE[6][5]

//===== Multigrid   
    uint n_se;       //NRE[6]             ///< number of subelements 
    
};



} //end namespace femus


#endif