#ifndef __femus_meshGencase_FEElemBase_hpp__
#define __femus_meshGencase_FEElemBase_hpp__

#include <string>
#include <vector>

#include "Typedefs.hpp"
#include "VBTypeEnum.hpp"
#include "GaussPoints.hpp"
#include "ElemType.hpp"

namespace femus {



class FEElemBase  {

public:

    FEElemBase();
    virtual ~FEElemBase();

// FE ==========
    static  FEElemBase* build(const std::string geomel_id_in, const uint fe_family);
    
// Multigrid ===   //only VV
    virtual float get_embedding_matrix(const uint,const uint,const uint) = 0;  //should be geom_el
    virtual double get_prol(const uint) = 0;
    
};


} //end namespace femus



#endif